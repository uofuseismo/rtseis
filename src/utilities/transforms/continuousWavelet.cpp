#include <iostream>
#include <string>
#include <vector>
#include <complex> // Put this before fftw
#include <fftw/fftw3.h>
#include "rtseis/utilities/transforms/dftRealToComplex.hpp"
#include "rtseis/utilities/transforms/continuousWavelet.hpp"
#include "rtseis/utilities/transforms/wavelets/derivativeOfGaussian.hpp"
#include "rtseis/utilities/filterImplementations/detrend.hpp"

using namespace RTSeis::Utilities::Transforms;

template<class T>
class ContinuousWavelet<T>::ContinuousWaveletImpl
{
public:
    void inverseTransform( )
    {

    }
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
        pImpl->mWavelet->evaluate(kWave.size(), scale, kWave.data(), &wPtr);
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
