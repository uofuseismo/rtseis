#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <climits>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/filterImplementations/decimate.hpp"
#include "rtseis/utilities/design/fir.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterImplementations/multiRateFIRFilter.hpp"

using namespace RTSeis::Utilities;
using namespace RTSeis::Utilities::FilterImplementations;

class Decimate::DecimateImpl
{
public:
    
    class MultiRateFIRFilter mMRFIRFilter;
    int mDownFactor = 1;
    int mPhase = 0;
    RTSeis::ProcessingMode mMode = RTSeis::ProcessingMode::POST_PROCESSING;
    RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
    bool mInitialized = false;
};

Decimate::Decimate() :
    pImpl(std::make_unique<DecimateImpl> ())
{
}

Decimate::Decimate(const Decimate &decimate)
{
    *this = decimate;
}

/*
Decimate::Decimate(Decimate &&decimate) noexcept
{
   *this = std::move(decimate);
}
*/

Decimate& Decimate::operator=(const Decimate &decimate)
{
    if (&decimate == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<DecimateImpl> (*decimate.pImpl);
    return *this;
}

/*
Decimate::Decimate(Decimate &&decimate) noexcept
{
    if (&decimate == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(decimate.pImpl);
    return *this;
}
*/

Decimate::~Decimate() = default;

void Decimate::clear() noexcept
{
    pImpl->mMRFIRFilter.clear();
    pImpl->mMode = RTSeis::ProcessingMode::POST_PROCESSING;
    pImpl->mPrecision = RTSeis::Precision::DOUBLE;
    pImpl->mDownFactor = 1;
    pImpl->mPhase = 0;
    pImpl->mInitialized = false;
}

void Decimate::initialize(const int downFactor,
                          const int filterLength,
                          const bool lremovePhaseShift,
                          const RTSeis::ProcessingMode mode,
                          const RTSeis::Precision precision)
{
    clear();
    if (downFactor < 1)
    {
        RTSEIS_THROW_IA("Downsampling factor = %d must be positive",
                        downFactor);
    }
    if (filterLength < 5)
    {
        RTSEIS_THROW_IA("Filter length = %d must be greater than 5",
                        filterLength);
    }
    // Set some properties
    pImpl->mDownFactor = downFactor;
    pImpl->mPrecision = precision;
    pImpl->mMode = mode;
    // This is a special case
    if (downFactor == 1)
    {
        pImpl->mInitialized = true;
        return;
    }
    int nfir = filterLength;
    // Postprocessing is a little trickier - may have to extend filter length
    if (mode == RTSeis::ProcessingMode::POST_PROCESSING && lremovePhaseShift)
    {
        bool lfail = true;
        for (auto k=0; k<INT_MAX-1; ++k)
        {
            int groupDelay = nfir/2;
            if ((groupDelay + 1)%downFactor == 0)
            {
                lfail = false;
                break;
            }
            nfir = nfir + 1;
        }
#ifdef DEBUG
        assert(!lfail);
#endif
        if (lfail){RTSEIS_THROW_RTE("%s", "Algorithmic failure");}
    } 
    // Create a hamming filter
    int order = nfir - 1; 
    auto r = 1.0/static_cast<double> (nfir);
    auto fir = FilterDesign::FIR::FIR1Lowpass(order, r,
                                              FilterDesign::FIRWindow::HAMMING);
    // Set the multirate FIR filter
    auto b = fir.getFilterTaps();
    int ntaps = b.size();
    constexpr int upFactor = 1;
    try
    {
        pImpl->mMRFIRFilter.initialize(upFactor, downFactor,
                                       ntaps, b.data(),
                                       mode, precision);
    }
    catch (std::exception &e)
    {
        RTSEIS_THROW_RTE("%sMultirate FIR filter initialization failed",
                         e.what());
    }
    pImpl->mInitialized = true;
}

bool Decimate::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

int Decimate::estimateSpace(const int n) const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (n < 0){RTSEIS_THROW_IA("n=%d cannot be negative", n);}
    /// Estimates the space required to hold the downsampled signal
    int phase = 0;
    if (pImpl->mMode ==  RTSeis::ProcessingMode::REAL_TIME)
    {
        phase = pImpl->mPhase;
    }
    int pDstLen = (n + pImpl->mDownFactor - 1 - phase)/pImpl->mDownFactor;
    return pDstLen;  
}
