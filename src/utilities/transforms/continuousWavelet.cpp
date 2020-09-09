#include <iostream>
#include <string>
#include <vector>
#include <complex> // Put this before fftw
#include <mkl.h>
#include <mkl_vsl.h>
#include <fftw/fftw3.h>
#include "rtseis/utilities/transforms/dftRealToComplex.hpp"
#include "rtseis/utilities/transforms/continuousWavelet.hpp"
#include "rtseis/utilities/transforms/wavelets/morlet.hpp"
#include "rtseis/utilities/filterImplementations/detrend.hpp"
#include "rtseis/utilities/math/convolve.hpp"
#include "private/convolve.hpp"

using namespace RTSeis::Utilities::Transforms;

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
void convolve(const int n, const float x[],
              const int nWidths, const float widths[],
              const Wavelets::IContinuousWavelet &wavelet,
              std::complex<float> *cwt,
              const bool direct = false)
{
}
*/

int convolve(const int n, const double x[],
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
     shared(cwt) \
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
    for (int i=0; i<nWidths; ++i)
    {
        // Evaluate the wavelet transform
        auto wPtr = reinterpret_cast<std::complex<double> *> (y);
        wavelet.evaluate(n, widths[i], &wPtr);
        // Formula is signal times conjugate of wavelet
        #pragma omp simd
        for (int j=0; j<n; ++j){wPtr[j] = std::conj(wPtr[j]);}
        // Convolve
        status = vslzConvExecX1D(task, y, yStride, z, zStride); 
        if (status != VSL_STATUS_OK)
        {
             std::cerr << "Convolution failed" << std::endl;
             ierr = 1;
        }
        // Only take the center component of the convolution and normalize
        // by 1/sqrt(a).  Note, part of dt should be baked into sqrt(a) through
        // the 1/sqrt(a) = sqrt(2*f*pi*dt/w0).  However, we're still off by 
        // a factor of sqrt(dt) which we additionally contribute here.
        auto xnorm = std::sqrt(dt)/std::sqrt(widths[i]);
        auto offset = i*n;
        auto zPtr = reinterpret_cast<std::complex<double> *> (z);
        #pragma omp simd 
        for (int j=0; j<nCopy; ++j)
        {
             cwt[offset+j] = xnorm*zPtr[i1 + j];
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

template<class T>
class ContinuousWavelet<T>::ContinuousWaveletImpl
{
public:
/*
    VSLConvTaskPtr mTask;
    bool mHaveTask = false;
    /// Forward transform signal
    DFTRealToComplex<T> mDFTR2C;
    /// Scaleogram scales (seconds/rad)
    std::vector<T> mScales;
    /// FFTw double precision plan
    fftw_plan  mDoublePlan;
    /// FFTw single precision plan
    fftwf_plan mFloatPlan;
*/
    /// Result
    std::vector<std::complex<T>> mCWT;
    /// Scales
    std::vector<T> mScales;
    /// Wavelet
    std::unique_ptr<Wavelets::IContinuousWavelet> mWavelet;
    /// Sampling period in seconds
    double mSamplingRate = 1;
    /// Number of samples
    int mSamples = 0;
    /// Initialized?
    bool mInitialized = false;
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

/// Destructor
template<class T>
ContinuousWavelet<T>::~ContinuousWavelet() = default;

template<class T>
void ContinuousWavelet<T>::clear() noexcept
{
    pImpl->mCWT.clear();
    pImpl->mScales.clear();
    pImpl->mWavelet = nullptr;
    pImpl->mSamplingRate = 1;
    pImpl->mSamples = 0;
    pImpl->mInitialized = false;
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
    pImpl->mCWT.resize(nSamples*nScales, std::complex<T> (0, 0));
    pImpl->mWavelet = wavelet.clone();
    pImpl->mScales.resize(nScales);
    std::copy(scales, scales + nScales, pImpl->mScales.data());
    pImpl->mSamplingRate = samplingRate;
    pImpl->mSamples = nSamples;
    pImpl->mInitialized = true;
}

/// Transform
template<class T>
void ContinuousWavelet<T>::transform(const int n, const T x[])
{
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
    auto cwt = pImpl->mCWT.data();
    auto dt = static_cast<T> (1./getSamplingRate());
    int option = 0; // Let algorithm decide
    auto error = convolve(nSamples, x, 
                          nScales, pImpl->mScales.data(),
                          *pImpl->mWavelet, cwt, option, dt);
    if (error != 0)
    {
        throw std::runtime_error("Error computing CWT");
    }

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

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class RTSeis::Utilities::Transforms::ContinuousWavelet<double>; 
//template class RTSeis::Utilities::Transforms::ContinuousWavelet<float>;
