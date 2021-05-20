#include <cstdio>
#include <cstdlib>
#include <complex>
#include <ipps.h>
#include "rtseis/utilities/transforms/welch.hpp"
#include "rtseis/utilities/transforms/utilities.hpp"
#include "rtseis/utilities/transforms/slidingWindowRealDFT.hpp"
#include "rtseis/utilities/transforms/slidingWindowRealDFTParameters.hpp"

using namespace RTSeis::Utilities::Transforms;

namespace
{

double computeSpectrumScaling(const int npts, const double samplingRate,
                              const double window[])
{
    double wsum2 = 0;
    ippsDotProd_64f(window, window, npts, &wsum2); 
    wsum2 = samplingRate*wsum2; 
    return wsum2;
}

double computeDensityScaling(const int npts, const double window[])
{
    double wsum;
    ippsSum_64f(window, npts, &wsum);
    wsum = wsum*wsum;
    return wsum;
}

IppStatus ippsMulC(const double *pSrc, const double val, double *pDst,
                   const int len)
{
    return ippsMulC_64f(pSrc, val, pDst, len);
}

IppStatus ippsMulC(const float *pSrc, const float val, float *pDst,
                   const int len)
{
    return ippsMulC_32f(pSrc, val, pDst, len);
}


}

template<class T>
class Welch<T>::WelchImpl
{
public:
    class SlidingWindowRealDFT<T> mSlidingWindowRealDFT;
    class SlidingWindowRealDFTParameters mParameters;
    std::vector<T> mSumSpectrum;
    double mSpectrumScaling = 1;
    double mDensityScaling = 1;
    double mSamplingRate = 1;
    bool mInitialized = false;
    bool mHaveTransform = false;
};

/// Constructors
template<class T>
Welch<T>::Welch() :
    pImpl(std::make_unique<WelchImpl> ())
{
}

/// Copy c'tor
template<class T>
Welch<T>::Welch(const Welch &welch)
{
    *this = welch;
}

/// Move c'tor
template<class T>
Welch<T>::Welch(Welch &&welch) noexcept
{
    *this = std::move(welch);
}

/// Copy assignment 
template<class T>
Welch<T>& Welch<T>::operator=(const Welch &welch)
{
    if (&welch == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<WelchImpl> (*welch.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Welch<T>& Welch<T>::operator=(Welch &&welch) noexcept
{
    if (&welch == this){return *this;}
    pImpl = std::move(welch.pImpl);
    return *this;
}

/// Destructor
template<class T>
Welch<T>::~Welch() = default;

/// Clears memory
template<class T>
void Welch<T>::clear() noexcept
{
    pImpl->mSlidingWindowRealDFT.clear();
    pImpl->mParameters.clear();
    pImpl->mSumSpectrum.clear();
    pImpl->mSpectrumScaling = 1;
    pImpl->mDensityScaling = 1;
    pImpl->mSamplingRate = 1;
    pImpl->mInitialized = false;
    pImpl->mHaveTransform = false;
}

/// Initialize
template<class T>
void Welch<T>::initialize(const SlidingWindowRealDFTParameters &parameters,
                          const double samplingRate)
{
    clear();
    if (samplingRate <= 0)
    {
        throw std::invalid_argument("samplingRate = "
                                   + std::to_string(samplingRate)
                                   + " must be positive");
    }
    if (!parameters.isValid())
    {
        throw std::invalid_argument("parameters class is invalid");
    }
    // Initialize the sliding window DFT
    pImpl->mSamplingRate = samplingRate;
    try
    {
        pImpl->mParameters = parameters;
        pImpl->mSlidingWindowRealDFT.initialize(pImpl->mParameters); 
    }
    catch (const std::exception &e)
    {
        clear();
        throw std::runtime_error("Failed to initialize DFT");
    }
    // Compute the scaling
    int nWindow = pImpl->mParameters.getWindowLength();
    auto window = pImpl->mParameters.getWindow();
    pImpl->mSpectrumScaling = computeSpectrumScaling(nWindow,
                                                     pImpl->mSamplingRate,
                                                     window.data());
    pImpl->mDensityScaling = computeDensityScaling(nWindow, window.data()); 
    // Set space for the intermediate output
    int nfreqs = pImpl->mSlidingWindowRealDFT.getNumberOfFrequencies();
    pImpl->mSumSpectrum.resize(nfreqs);
    pImpl->mInitialized = true;
}

/// Initialized?
template<class T>
bool Welch<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Have transform?
template<class T>
bool Welch<T>::haveTransform() const noexcept
{
    return pImpl->mHaveTransform;
}

/// Have number of frequncies?
template<class T>
int Welch<T>::getNumberOfFrequencies() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mSlidingWindowRealDFT.getNumberOfFrequencies();
}

/// Get number of samples
template<class T>
int Welch<T>::getNumberOfSamples() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mSlidingWindowRealDFT.getNumberOfSamples();

}

/// Actually compute welch transform
template<>
void Welch<double>::transform(const int nSamples, const double x[])
{
    pImpl->mHaveTransform = true;
    int nSamplesRef = getNumberOfSamples(); // Throws if not inittialized
    if (nSamples != nSamplesRef)
    {
        throw std::invalid_argument("nSamples = " + std::to_string(nSamples)
                                  + " must equal "
                                  + std::to_string(nSamplesRef));
    }
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    pImpl->mSlidingWindowRealDFT.transform(nSamples, x);
    // Sum over the frequency bins
    auto nFrequencies = getNumberOfFrequencies();
    auto nWindows = pImpl->mSlidingWindowRealDFT.getNumberOfTransformWindows();
    // Initialize the summation
    double *pSumSpectrum = pImpl->mSumSpectrum.data();
    auto cPtr = pImpl->mSlidingWindowRealDFT.getTransform(0);
    auto pDFT = reinterpret_cast<const Ipp64fc *> (cPtr);
    #pragma omp simd
    for (auto k=0; k<nFrequencies; ++k)
    {
        pSumSpectrum[k] = pDFT[k].re*pDFT[k].re
                        + pDFT[k].im*pDFT[k].im;
    }
    // And sum the other windows
    for (auto i=1; i<nWindows; ++i)
    {
        cPtr = pImpl->mSlidingWindowRealDFT.getTransform(i);
        pDFT = reinterpret_cast<const Ipp64fc *> (cPtr);
        #pragma omp simd
        for (auto k=0; k<nFrequencies; ++k)
        {
            pSumSpectrum[k] = pSumSpectrum[k]
                            + pDFT[k].re*pDFT[k].re
                            + pDFT[k].im*pDFT[k].im;
        }
    }
    pImpl->mHaveTransform = true;
}

/// Get the power spectrum
template<class T>
std::vector<T> Welch<T>::getPowerSpectrum() const
{
    auto nFreqs = getNumberOfFrequencies(); // Will throw initializaiton error
    if (!haveTransform())
    {
        throw std::runtime_error("welch transform not yet computed");
    }
    std::vector<T> spectrum(nFreqs);
    auto sPtr = spectrum.data();
    getPowerSpectrum(spectrum.size(), &sPtr);
    return spectrum;
}

/// Get the power spectrum density
template<class T>
std::vector<T> Welch<T>::getPowerSpectralDensity() const
{
    auto nFreqs = getNumberOfFrequencies(); // Will throw initializaiton error
    if (!haveTransform())
    {
        throw std::runtime_error("welch transform not yet computed");
    }
    std::vector<T> spectrum(nFreqs);
    auto sPtr = spectrum.data();
    getPowerSpectralDensity(spectrum.size(), &sPtr);
    return spectrum;
}

/// Get power spectrum
template<class T>
void Welch<T>::getPowerSpectrum(const int nFrequencies,
                                T *powerSpectrum[]) const
{
    auto nFreqs = getNumberOfFrequencies(); // Will throw initializaiton error
    if (nFrequencies != nFreqs)
    {
        throw std::invalid_argument("nFrequencies = "
                                  + std::to_string(nFrequencies)
                                  + " must equal " + std::to_string(nFreqs));
    }
    T *ptr = *powerSpectrum;
    if (ptr == nullptr){throw std::invalid_argument("powerSpectrum is NULL");}
    if (!haveTransform())
    {
        throw std::runtime_error("welch transform not yet computed");
    }
    // Copy and scale
    auto nWindows = pImpl->mSlidingWindowRealDFT.getNumberOfTransformWindows();
    T xscal = 0;
    if (nFreqs > 0)
    {
        xscal = static_cast<T> (2.0/(nWindows*pImpl->mSpectrumScaling));
    }
    ippsMulC(pImpl->mSumSpectrum.data(), xscal, ptr, nFreqs);
}

/// Get power spectral density
template<class T>
void Welch<T>::getPowerSpectralDensity(const int nFrequencies,
                                       T *psd[]) const
{
    auto nFreqs = getNumberOfFrequencies(); // Will throw initializaiton error
    if (nFrequencies != nFreqs)
    {
        throw std::invalid_argument("nFrequencies = "
                                  + std::to_string(nFrequencies)
                                  + " must equal " + std::to_string(nFreqs));
    }
    T *ptr = *psd;
    if (ptr == nullptr){throw std::invalid_argument("psd is NULL");}
    if (!haveTransform())
    {
        throw std::runtime_error("welch transform not yet computed");
    }
    auto nWindows = pImpl->mSlidingWindowRealDFT.getNumberOfTransformWindows();
    // Copy and scale
    T xscal = 0;
    if (nFreqs > 0)
    {
        xscal = static_cast<T> (2.0/(nWindows*pImpl->mDensityScaling));
    }
    ippsMulC(pImpl->mSumSpectrum.data(), xscal, ptr, nFreqs);
}

/// Get frequencies
template<class T>
std::vector<T> Welch<T>::getFrequencies() const
{
    return pImpl->mSlidingWindowRealDFT.getFrequencies(pImpl->mSamplingRate);
}

/// Get frequencies
template<class T>
void Welch<T>::getFrequencies(const int nFrequencies, T *freqsIn[]) const
{
    pImpl->mSlidingWindowRealDFT.getFrequencies(pImpl->mSamplingRate,
                                                nFrequencies, freqsIn);
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class RTSeis::Utilities::Transforms::Welch<double>;
template class RTSeis::Utilities::Transforms::Welch<float>;
