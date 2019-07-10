#include <cstdio>
#include <cstdlib>
#include <complex>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
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

}

class Welch::WelchImpl
{
public:
    class SlidingWindowRealDFT mSlidingWindowRealDFT;
    class SlidingWindowRealDFTParameters mParameters;
    std::vector<double> mSumSpectrum;
    double mSpectrumScaling = 1;
    double mDensityScaling = 1;
    double mSamplingRate = 1;
    bool mInitialized = false;
    bool mHaveTransform = false;
};

/// Constructors
Welch::Welch() :
    pImpl(std::make_unique<WelchImpl> ())
{
}

Welch::Welch(const Welch &welch)
{
    *this = welch;
}

Welch::Welch(Welch &&welch) noexcept
{
    *this = std::move(welch);
}

/// Operators
Welch& Welch::operator=(const Welch &welch)
{
    if (&welch == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<WelchImpl> (*welch.pImpl);
    return *this;
}

Welch& Welch::operator=(Welch &&welch) noexcept
{
    if (&welch == this){return *this;}
    pImpl = std::move(welch.pImpl);
    return *this;
}

/// Destructor
Welch::~Welch() = default;

/// Clears memory
void Welch::clear() noexcept
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
void Welch::initialize(const SlidingWindowRealDFTParameters &parameters,
                       const double samplingRate)
{
    clear();
    if (samplingRate <= 0)
    {
        RTSEIS_THROW_IA("samplingRate = %lf must be posistive", samplingRate);
    }
    if (!parameters.isValid())
    {
        RTSEIS_THROW_IA("%s", "parameters are not valid");
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
        RTSEIS_THROW_RTE("%s", "Failed to initialize sliding DFT");
    }
    // Compute the scaling
    int nWindow = pImpl->mParameters.getWindowLength();
    std::vector<double> window  = pImpl->mParameters.getWindow();
    pImpl->mSpectrumScaling = computeSpectrumScaling(nWindow,
                                                     pImpl->mSamplingRate,
                                                     window.data());
    pImpl->mDensityScaling = computeDensityScaling(nWindow, window.data()); 
    // Set space for the intermediate output
    int nfreqs = pImpl->mSlidingWindowRealDFT.getNumberOfFrequencies();
    pImpl->mSumSpectrum.resize(nfreqs);
    pImpl->mInitialized = true;
}

bool Welch::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

bool Welch::haveTransform() const noexcept
{
    return pImpl->mHaveTransform;
}

int Welch::getNumberOfFrequencies() const
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not initialized");
    }
    return pImpl->mSlidingWindowRealDFT.getNumberOfFrequencies();
}

int Welch::getNumberOfSamples() const
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not initialized");
    }
    return pImpl->mSlidingWindowRealDFT.getNumberOfSamples();

}
void Welch::transform(const int nSamples, const double x[])
{
    pImpl->mHaveTransform = true;
    int nSamplesRef = getNumberOfSamples(); // Throws if not inittialized
    if (nSamples != nSamplesRef)
    {
        RTSEIS_THROW_IA("nSamples = %d must equal %d", nSamples, nSamplesRef);
    }
    if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
    pImpl->mSlidingWindowRealDFT.transform(nSamples, x);
    // Sum over the frequency bins
    auto nFrequencies = getNumberOfFrequencies();
    auto nWindows = pImpl->mSlidingWindowRealDFT.getNumberOfTransformWindows();
    // Initialize the summation
    double *pSumSpectrum = pImpl->mSumSpectrum.data();
    auto cPtr = pImpl->mSlidingWindowRealDFT.getTransform64f(0);
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
        cPtr = pImpl->mSlidingWindowRealDFT.getTransform64f(i);
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

void Welch::getPowerSpectrum(const int nFrequencies,
                             double *powerSpectrum[]) const
{
    auto nFreqs = getNumberOfFrequencies(); // Will throw initializaiton error
    if (nFrequencies != nFreqs)
    {
        RTSEIS_THROW_IA("nFrequencies = %d must equal %d",
                        nFrequencies, nFreqs);
    }
    double *ptr = *powerSpectrum;
    if (ptr == nullptr){RTSEIS_THROW_IA("%s", "powerSpectrum is NULL");}
    if (!haveTransform())
    {
        RTSEIS_THROW_RTE("%s", "welch transform not yet computed");
    }
    // Copy and scale
    auto nWindows = pImpl->mSlidingWindowRealDFT.getNumberOfTransformWindows();
    double xscal = 0;
    if (nFreqs > 0){xscal = 2.0/(nWindows*pImpl->mSpectrumScaling);}
    ippsMulC_64f(pImpl->mSumSpectrum.data(), xscal, ptr, nFreqs);
}

void Welch::getPowerSpectralDensity(const int nFrequencies, double *psd[]) const
{
    auto nFreqs = getNumberOfFrequencies(); // Will throw initializaiton error
    if (nFrequencies != nFreqs)
    {
        RTSEIS_THROW_IA("nFrequencies = %d must equal %d",
                        nFrequencies, nFreqs);
    }
    double *ptr = *psd;
    if (ptr == nullptr){RTSEIS_THROW_IA("%s", "psd is NULL");}
    if (!haveTransform())
    {
        RTSEIS_THROW_RTE("%s", "welch transform not yet computed");
    }
    auto nWindows = pImpl->mSlidingWindowRealDFT.getNumberOfTransformWindows();
    // Copy and scale
    double xscal = 0;
    if (nFreqs > 0){xscal = 2.0/(nWindows*pImpl->mDensityScaling);}
    ippsMulC_64f(pImpl->mSumSpectrum.data(), xscal, ptr, nFreqs);
}

void Welch::getFrequencies(const int nFrequencies, double *freqsIn[]) const
{
    auto nFreqs = getNumberOfFrequencies(); // Will throw initialization error
    if (nFrequencies != nFreqs)
    {
        RTSEIS_THROW_IA("nFrequencies = %d must equal %d",
                        nFrequencies, nFreqs); 
    }
    double *freqs = *freqsIn;
    if (freqs == nullptr)
    {
        RTSEIS_THROW_IA("%s", "frequencies is NULL");
    }
    int nSamples = pImpl->mParameters.getDFTLength();
    try
    {
        DFTUtilities::realToComplexDFTFrequencies(nSamples,
                                                  1.0/pImpl->mSamplingRate,
                                                  nFreqs,
                                                  freqsIn); 
    }
    catch (const std::exception &e)
    {
#ifdef DEBUG
        assert(false);
#else
        RTSEIS_THROW_RTE("%s", "Failed to call realToComplexDFTFrequencies");
#endif
    }
}
