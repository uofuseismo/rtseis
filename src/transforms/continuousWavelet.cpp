#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <complex> // Put this before fftw
#include <mkl.h>
#include <mkl_vsl.h>
#include <fftw/fftw3.h>
#include "rtseis/transforms/dftRealToComplex.hpp"
#include "rtseis/transforms/continuousWavelet.hpp"
#include "rtseis/transforms/wavelets/morlet.hpp"
#include "rtseis/utilities/math/convolve.hpp"
#include "private/convolve.hpp"
#include "private/pad.hpp"

using namespace RTSeis::Transforms;

namespace Convolve = RTSeis::Utilities::Math::Convolve;

namespace
{

template<class T>
void realToComplex(const int n, const T *x, std::complex<T> *z)
{
    #pragma omp simd
    for (int i=0; i<n; ++i)
    {
        z[i] = std::complex<T> (x[i], 0);
    }
}

/*
int convolve(const int n, const int ldx, const double x[],
             const int nWidths, const double widths[],
             const Wavelets::IContinuousWavelet &wavelet,
             std::complex<double> *cwt,
             const int option,
             const double dt)
{
}
*/

int convolve(const int n, const int ldx, const double x[],
             const int nWidths, const double widths[],
             const Wavelets::IContinuousWavelet &wavelet,
             std::complex<double> *cwt,
             const int option,
             const double dt)
{
    // Get correlation length
    auto indices = computeTrimIndices(Convolve::Mode::SAME, n, n);
    auto nCopy = Convolve::computeConvolutionLength(n, n, Convolve::Mode::SAME); 
    auto i1 = indices.first;
#ifndef NDEBUG
    auto i2 = indices.second;
    assert(i2 - i1 == nCopy);
#endif
    // Figure out the convolution mode
    MKL_INT mode = VSL_CONV_MODE_AUTO;
    if (option == 1){mode = VSL_CONV_MODE_DIRECT;}
    if (option == 2){mode = VSL_CONV_MODE_FFT;}
    // Set space on each thread
    int ierr = 0;
    #pragma omp parallel \
     shared(cwt, x, std::cerr, wavelet, widths) \
     firstprivate(i1, nCopy, mode) \
     reduction(+ : ierr) \
     default(none)
    {
    // Set space for the transform
    auto *xz = static_cast<MKL_Complex16 *>
               (mkl_calloc(static_cast<int> (n), sizeof(MKL_Complex16), 64));
    auto xzPtr = reinterpret_cast<std::complex<double> *> (xz);
    realToComplex(n, x, xzPtr); // Copy real to complex 
    // Since x won't change - set it here.  If using an FFT this will let MKL
    // potentially precompute the FFT once.
    VSLConvTaskPtr task;
    const MKL_INT xStride = 1;
    const MKL_INT yStride = 1;
    const MKL_INT zStride = 1;
    auto xShape = static_cast<MKL_INT> (n);
    auto yShape = static_cast<MKL_INT> (n);
    auto zShape = xShape + yShape - 1;
    auto status = vslzConvNewTaskX1D(&task, mode, xShape, yShape, zShape,
                                     xz, xStride);
    if (status != VSL_STATUS_OK)
    {
        std::cerr << "Failed to initialize convolution engine" << std::endl;
        ierr = 1;
    }
    // Set space for wavelet
    auto *y = static_cast<MKL_Complex16 *> 
              (mkl_calloc(yShape, sizeof(MKL_Complex16), 64));
    auto *z = static_cast<MKL_Complex16 *>
              (mkl_calloc(zShape, sizeof(MKL_Complex16), 64));
    // Have each thread contribute to the correlation
    #pragma omp for
    for (int j=0; j<nWidths; ++j)
    {
        // Evaluate the wavelet transform
        auto wPtr = reinterpret_cast<std::complex<double> *> (y);
        wavelet.evaluate(n, widths[j], &wPtr);
        // Formula is signal convolved with time reversed conjugate of wavelet
        #pragma omp simd
        for (int i=0; i<n/2; ++i)
        {
            auto temp = std::conj(wPtr[n-1-i]);
            wPtr[n-1-i] = std::conj(wPtr[i]);
            wPtr[i] = temp;
        }
/*
// Slow convolution
std::vector<double> wRe(n);
std::vector<double> wIm(n);
std::vector<double> xRe(n);
for (int i=0; i<n; ++i)
{
  wRe[i] = std::real(wPtr[i]);
  wIm[i] = std::imag(wPtr[i]);
  xRe[i] = x[i];
} 
std::cout << n << " " << widths[j] << " " << 1./std::sqrt(widths[j])*wPtr[500]
                                   << " " << 1./std::sqrt(widths[j])*wPtr[510] << std::endl;
auto xwRe = Convolve::convolve(xRe, wRe, Convolve::Mode::SAME);
auto xwIm = Convolve::convolve(xRe, wIm, Convolve::Mode::SAME);
std::vector<std::complex<double>> xw(n);
for (int i=0; i<n; ++i)
{
    xw[i] = std::complex<double> (xwRe[i], xwIm[i]);
}
*/
        // Convolve
        status = vslzConvExecX1D(task, y, yStride, z, zStride); 
        if (status != VSL_STATUS_OK)
        {
             std::cerr << "Convolution failed" << std::endl;
             ierr = 1;
        }
        // Only take the center component of the convolution and normalize
        // by 1/sqrt(|a|) ala the formula.  Note, MKL only implements the
        // convolutional sum; to make it match the convolutional integral
        // we multiply by dt.
        auto xnorm = dt/std::sqrt(std::abs(widths[j]));
        auto cwtPtr = cwt + j*ldx;
        auto zPtr = reinterpret_cast<std::complex<double> *> (z);
        //double emax = 0;
        #pragma omp simd aligned(cwtPtr : 64)
        for (int i=0; i<nCopy; ++i)
        {
             //emax = std::max(emax, std::abs(zPtr[i1 + i] - xw[i]));
             cwtPtr[i] = xnorm*zPtr[i1 + i];
        }
    }
    // Clean up
    vslConvDeleteTask(&task);
    MKL_free(xz);
    MKL_free(y);
    MKL_free(z);
    } // End parallel 
    if (ierr != 0)
    {
        std::cerr << "Errors detected during convolution" << std::endl;
    }
    return ierr;
}

}

///--------------------------------------------------------------------------///
///                            Pointer to Implementation                     ///
///--------------------------------------------------------------------------///
template<class T>
class ContinuousWavelet<T>::ContinuousWaveletImpl
{
public:
    /// Default c'tor
    ContinuousWaveletImpl() = default;
    /// Copy c'tor
    ContinuousWaveletImpl(const ContinuousWaveletImpl &cwt)
    {
        *this = cwt;
    }
    /// Copy assignment
    ContinuousWaveletImpl& operator=(const ContinuousWaveletImpl &cwt)
    {
        if (&cwt == this){return *this;}
        mScales = cwt.mScales;
        if (cwt.mWavelet){mWavelet = cwt.mWavelet->clone();}
        mSamplingRate = cwt.mSamplingRate;
        mSamples = cwt.mSamples;
        mLeadingDimension = cwt.mLeadingDimension;
        mHaveTransform = cwt.mHaveTransform;
        mInitialized = cwt.mInitialized;
        auto nScales = static_cast<int> (mScales.size());
        if (mLeadingDimension > 0 && nScales > 0)
        {
            auto num = static_cast<size_t> (mLeadingDimension)
                      *static_cast<size_t> (nScales);
            if (mDoublePrecision)
            {
                mCWT = mkl_calloc(num, sizeof(MKL_Complex16), 64);
                std::memcpy(mCWT, cwt.mCWT, num*sizeof(MKL_Complex16));
            }
            else
            {
                mCWT = mkl_calloc(num, sizeof(MKL_Complex8), 64);
                std::memcpy(mCWT, cwt.mCWT, num*sizeof(MKL_Complex8)); 
            }
        }
        return *this;
    }
    /// Destructor
    ~ContinuousWaveletImpl()
    {
        clear();
    }
    /// Clear the class
    void clear() noexcept
    {
        if (mCWT){MKL_free(mCWT);}
        mCWT = nullptr;
        mScales.clear();
        mWavelet = nullptr;
        mSamplingRate = 1;
        mSamples = 0;
        mLeadingDimension = 0;
        mHaveTransform = false;
        mInitialized = false;
    } 
    /// CWT.  This is an [nScales x mLeadingDimension] row major matrix.
    //std::vector<std::complex<T>> mCWT;
    void *mCWT = nullptr;
    /// Scales
    std::vector<T> mScales;
    /// Wavelet
    std::unique_ptr<Wavelets::IContinuousWavelet> mWavelet;
    /// Sampling period in seconds
    double mSamplingRate = 1;
    /// Number of samples
    int mSamples = 0;
    /// Leading dimension of mCWT
    int mLeadingDimension = 0;
    /// Have CWT?
    bool mHaveTransform = false;
    /// Initialized?
    bool mInitialized = false;
    /// Double precision?
    const bool mDoublePrecision = (sizeof(T) == 8);
};

///--------------------------------------------------------------------------///
///                           End Implementation                             ///
///--------------------------------------------------------------------------///
 
/// C'tor
template<class T>
ContinuousWavelet<T>::ContinuousWavelet() :
    pImpl(std::make_unique<ContinuousWaveletImpl> ())
{
}

/// Copy c'tor
template<class T>
ContinuousWavelet<T>::ContinuousWavelet(const ContinuousWavelet &cwt)
{
    *this = cwt;
}

/// Move c'tor
template<class T>
ContinuousWavelet<T>::ContinuousWavelet(ContinuousWavelet &&cwt) noexcept
{
    *this = std::move(cwt);
}

/// Copy assignment
template<class T>
ContinuousWavelet<T>& 
ContinuousWavelet<T>::operator=(const ContinuousWavelet &cwt)
{
    if (&cwt == this){return *this;}
    pImpl = std::make_unique<ContinuousWaveletImpl> (*cwt.pImpl);
    return *this;
}

/// Move assignment
template<class T>
ContinuousWavelet<T>& 
ContinuousWavelet<T>::operator=(ContinuousWavelet &&cwt) noexcept
{
    if (&cwt == this){return *this;}
    pImpl = std::move(cwt.pImpl);
    return *this;
}

/// Destructor
template<class T>
ContinuousWavelet<T>::~ContinuousWavelet() = default;

/// Release memory on the class
template<class T>
void ContinuousWavelet<T>::clear() noexcept
{
    pImpl->clear();
}

/// Initialize
template<class T>
void ContinuousWavelet<T>::initialize(
    const int nSamples, const int nScales, const double scales[],
    const Wavelets::IContinuousWavelet &wavelet, const double samplingRate)
{
    clear();
    if (nSamples < 1){throw std::invalid_argument("nSamples must be positive");}
    if (nScales < 1){throw std::invalid_argument("nScales must be positive");}
    if (scales == nullptr){throw std::invalid_argument("scales is NULL");}
    if (samplingRate <= 0)
    {
        throw std::invalid_argument("Sampling rate must be positive");
    }
    int ldx = padLength(nSamples, sizeof(std::complex<T>), 64);
    auto num = static_cast<size_t> (ldx)*static_cast<size_t> (nScales);
    if (pImpl->mDoublePrecision)
    {
        pImpl->mCWT = mkl_calloc(num, sizeof(MKL_Complex16), 64);
    }
    else
    {
        pImpl->mCWT = mkl_calloc(num, sizeof(MKL_Complex8), 64);
    }
    pImpl->mWavelet = wavelet.clone();
    pImpl->mScales.resize(nScales);
    std::copy(scales, scales + nScales, pImpl->mScales.data());
    pImpl->mSamplingRate = samplingRate;
    pImpl->mSamples = nSamples;
    pImpl->mLeadingDimension = ldx;
    pImpl->mHaveTransform = false;
    pImpl->mInitialized = true;
}

/// Transform
template<class T>
void ContinuousWavelet<T>::transform(const int n, const T x[])
{
    pImpl->mHaveTransform = false;
    int nSamples = getNumberOfSamples(); // Throws on initialized
    if (n != nSamples)
    {
        throw std::invalid_argument("Number of samples in x = "
                                  + std::to_string(n) + " must equal "
                                  + std::to_string(nSamples));
    }
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    // Do the work
    auto nScales = getNumberOfScales();
    auto cwt = reinterpret_cast<std::complex<T> *> (pImpl->mCWT);
    auto dt = static_cast<T> (1./getSamplingRate());
    int option = 0; // Let algorithm decide
    auto error = convolve(nSamples, pImpl->mLeadingDimension, x, 
                          nScales, pImpl->mScales.data(),
                          *pImpl->mWavelet, cwt, option, dt);
    if (error != 0)
    {
        throw std::runtime_error("Error computing CWT");
    }
    pImpl->mHaveTransform = true;
}

/// Number of samples
template<class T>
int ContinuousWavelet<T>::getNumberOfSamples() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mSamples;
}

/// Number of scales
template<class T>
int ContinuousWavelet<T>::getNumberOfScales() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return static_cast<int> (pImpl->mScales.size());
}

/// Sampling rate
template<class T>
double ContinuousWavelet<T>::getSamplingRate() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mSamplingRate;
}

/// Initialized?
template<class T>
bool ContinuousWavelet<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Have transform?
template<class T>
bool ContinuousWavelet<T>::haveTransform() const noexcept
{
    return pImpl->mHaveTransform;
}

/// Get the transform
template<class T>
void ContinuousWavelet<T>::getTransform(
    const int nSamples, const int nScales, std::complex<T> *cwtOut[]) const
{
    if (!haveTransform()){throw std::runtime_error("CWT not yet computed");}
    auto n = getNumberOfSamples();
    auto ns = getNumberOfScales();
    auto cwt = *cwtOut;
    if (n != nSamples)
    {
        throw std::invalid_argument("nSamples = " + std::to_string(nSamples)
                                  + " must equal " + std::to_string(n));
    }
    if (ns != nScales)
    {
        throw std::invalid_argument("nScales = " + std::to_string(nScales)
                                  + " must equal " + std::to_string(ns));
    }
    if (cwt == nullptr){throw std::invalid_argument("cwt it NULL");}
    auto ldx = pImpl->mLeadingDimension;
    auto cwtPtr = reinterpret_cast<std::complex<T> *> (pImpl->mCWT);
    for (int i=0; i<ns; ++i)
    {
        auto i1 = i*ldx; 
        auto i2 = i*ldx + nSamples;
        auto j1 = i*nSamples;
        std::copy(cwtPtr + i1, cwtPtr + i2, cwt + j1);
    } 
}

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class RTSeis::Transforms::ContinuousWavelet<double>; 
//template class RTSeis::Utilities::Transforms::ContinuousWavelet<float>;
