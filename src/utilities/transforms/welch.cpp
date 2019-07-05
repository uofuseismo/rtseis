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
    wsum2 = 1.0/(samplingRate*wsum2); 
    return wsum2;
}

double computeDensityScaling(const int npts, const double window[])
{
    double wsum;
    ippsSum_64f(window, npts, &wsum);
    wsum = 1.0/(wsum*wsum);
    return wsum;
}

}

class Welch::WelchImpl
{
public:
    class SlidingWindowRealDFT mSlidingWindowRealDFT;
    class SlidingWindowRealDFTParameters mParameters;
    double mSamplingRate = 1;
    bool mInitialized = false;
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
    pImpl->mSamplingRate = 1;
    pImpl->mInitialized = false;
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
    int nfreqs = pImpl->mSlidingWindowRealDFT.getNumberOfFrequencies();
    if (nfreqs > 1)
    {
        double fNyq = pImpl->mSamplingRate/2.0;
        double df = fNyq/static_cast<double> (nfreqs - 1);
        std::vector<double> frequencies;
        #pragma omp simd
        for (auto i=0; i<nfreqs; ++i)
        {
            frequencies[i] = static_cast<double> (i)*df;
        }
    }
    pImpl->mInitialized = true;
}

bool Welch::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

int Welch::getNumberOfFrequencies() const
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not initialized");
    }
    return pImpl->mSlidingWindowRealDFT.getNumberOfFrequencies();
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
    int nSamples = 2*(nFreqs - 1); 
    DFTUtilities::realToComplexDFTFrequencies(nSamples, nFreqs,
                                              1.0/pImpl->mSamplingRate,
                                              freqsIn); 
}
