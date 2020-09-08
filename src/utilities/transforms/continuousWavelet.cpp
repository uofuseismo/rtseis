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

void convolve(const int n, const double x[],
              const int nWidths, const double widths[],
              const Wavelets::IContinuousWavelet &wavelet,
              std::complex<double> *cwt,
              const bool direct = false)
{
    MKL_INT mode = VSL_CONV_MODE_FFT;
    auto indices = computeTrimIndices(Convolve::Mode::SAME, n, n);
    auto nCopy = Convolve::computeConvolutionLength(n, n, Convolve::Mode::SAME);
#ifndef NDEBUG
    auto i2 = indices.second;
    assert(i2 - i1 == nCopy);
#endif
    if (direct){mode = VSL_CONV_MODE_DIRECT;}
    int ierr = 0;
    // Set space on each thread
    #pragma omp parallel \
     shared(cwt) \
     firstprivate(indices, nCopy, mode) \
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
        //wavelet.evaluate(n, widths[i], &wPtr);
        // Convolve
        vslzConvExecX1D(task, y, yStride, z, zStride); 
        // Only take the center component of the convolution
        auto offset = i*n;
        auto zPtr = reinterpret_cast<std::complex<double> *> (z);
        std::copy(zPtr + indices.first, zPtr + indices.second,
                  cwt + offset);
    }
    // Clean up
    vslConvDeleteTask(&task);
    MKL_free(xz);
    MKL_free(y);
    MKL_free(z);
    } // End parallel 
}

}

template<class T>
class ContinuousWavelet<T>::ContinuousWaveletImpl
{
public:
    ~ContinuousWaveletImpl()
    {
        clear();
    } 
    void clear() noexcept
    {
        if (mHaveTask){vslConvDeleteTask(&mTask);}
        mHaveTask = false;
    }
    void inverseTransform( )
    {

    }
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
    /// Wavelet
    std::unique_ptr<Wavelets::IContinuousWavelet> mWavelet;
    /// Sampling period in seconds
    double mSamplingPeriod = 1;
    /// Number of samples
    int mSamples = 0;
    /// Number of DFT samples
    int mDFTSamples = 0; 
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
    // Demean the signal
    std::vector<T> xDemeaned(n);
    auto xPtr = xDemeaned.data();
    T xmean;
    RTSeis::Utilities::FilterImplementations::removeMean(n, x, &xPtr, &xmean);
    // Fourier transform the signal
    auto lendft = pImpl->mDFTR2C.getTransformLength();
    std::vector<std::complex<T>> xDFT(lendft);
    auto xDFTPtr = xDFT.data();
    pImpl->mDFTR2C.forwardTransform(n, xDemeaned.data(), lendft, &xDFTPtr);
    // Wavenumber array (Eqn 5)
    auto dOmega = 2*M_PI/(pImpl->mDFTSamples*pImpl->mSamplingPeriod);
    std::vector<double> kWave(pImpl->mDFTSamples, 0);
    for (int i=1; i<pImpl->mDFTSamples/2+1; ++i)
    {
        kWave[i] = i*dOmega;
    }
  

    // Main wavelet loop on scales
    std::vector<std::complex<T>> w(kWave.size(), std::complex<T> (0, 0));
    for (const auto &scale : pImpl->mScales)
    {
        auto wPtr = w.data();
        //pImpl->mWavelet->evaluate(kWave.size(), scale, kWave.data(), &wPtr);
    }
}

/// Number of samples
template<class T>
int ContinuousWavelet<T>::getNumberOfSamples() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mSamples;
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
