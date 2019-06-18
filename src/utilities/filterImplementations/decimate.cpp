#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <climits>
#include <ipps.h>
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
    int mGroupDelay = 0;
    int mFIRLength = 0;
    //int mPhase = 0;
    RTSeis::ProcessingMode mMode = RTSeis::ProcessingMode::POST_PROCESSING;
    RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
    bool mRemovePhaseShift = false;
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
    pImpl->mGroupDelay = 0;
    pImpl->mFIRLength = 0;
    pImpl->mRemovePhaseShift = false;
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
/*
    // This is a special case
    constexpr int upFactor = 1;
    if (downFactor == 1)
    {
        constexpr int downFactor1 = 1;
        constexpr int ntaps1 = 1;
        constexpr double b1[1] = {1};
        pImpl->mMRFIRFilter.initialize(upFactor, downFactor1,
                                       ntaps1, b1,
                                       mode, precision);
        pImpl->mInitialized = true;
        return;
    }
*/
    int nfir = filterLength;
    // Postprocessing is a little trickier - may have to extend filter length
    if (mode == RTSeis::ProcessingMode::POST_PROCESSING && lremovePhaseShift)
    {
        // N.B.  The idea is that we think of decimation as filtering then
        // downsampling.  The filter will introduce a group delay of
        // groupDelay = nfir/2 where nfir is odd (i.e., a unit impulse
        // introduce 1/2 = 0 delay).  The delay is then removed by
        // extracing the filtered signal at yfiltered[groupDelay].
        // However, we have downsampled, so we really want
        // yfiltered[groupDelay/downFactor].  Therefore, our padding
        // has 2 conditions:
        //  1.  nfir be odd.
        //  2.  groupDelay be evenly divisible by downfactor. 
        bool lfail = true;
        for (auto k=0; k<INT_MAX-1; ++k)
        {
            int groupDelay = nfir/2;
            if (groupDelay%downFactor == 0 && nfir%2 == 1)
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
        pImpl->mGroupDelay = nfir/2;
        pImpl->mRemovePhaseShift = true;
    } 
    // Create a hamming filter
    pImpl->mFIRLength = nfir;
    int order = nfir - 1; 
    auto r = 1.0/static_cast<double> (nfir);
    auto fir = FilterDesign::FIR::FIR1Lowpass(order, r,
                                              FilterDesign::FIRWindow::HAMMING);
    // Set the multirate FIR filter
    auto b = fir.getFilterTaps();
    constexpr int upFactor = 1;
    int ntaps = b.size();
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
    if (pImpl->mRemovePhaseShift)
    {
    }
    else
    {
    }
    return pImpl->mMRFIRFilter.estimateSpace(n);
}

int Decimate::getInitialConditionLength() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mMRFIRFilter.getInitialConditionLength(); 
}

void Decimate::setInitialConditions(const int nz,const double zi[])
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    int nzref = getInitialConditionLength();
    if (nz != nzref)
    {
        RTSEIS_THROW_IA("nz = %d must equal %d", nz, nzref);
    }
    if (nz > 0 && zi == nullptr)
    {
        RTSEIS_THROW_IA("%s", "zi cannot be NULL");
    }
    pImpl->mMRFIRFilter.setInitialConditions(nz, zi);
}

void Decimate::resetInitialConditions()
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    pImpl->mMRFIRFilter.resetInitialConditions();
}

void Decimate::apply(const int nx, const double x[],
                     const int ny, int *nyDown, double *yIn[])
{
    *nyDown = 0;
    if (nx <= 0){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
    int nyref = estimateSpace(nx);
    if (ny < nyref){RTSEIS_THROW_IA("ny = %d must be at least %d", ny, nyref);}
    double *y = *yIn;
    if (y == nullptr)
    {
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    // Special case
    if (pImpl->mRemovePhaseShift)
    {
        // Pad a copy of x with group delay samples
        int npad = nx + pImpl->mGroupDelay;
        double *xpad = ippsMalloc_64f(npad);
        ippsCopy_64f(x, xpad, nx); 
        ippsZero_64f(&xpad[nx], pImpl->mGroupDelay);
        // Filter
        int nywork = pImpl->mMRFIRFilter.estimateSpace(npad);
        double *ypad = ippsMalloc_64f(nywork);
        int nyworkDown;
        pImpl->mMRFIRFilter.apply(npad, xpad, nywork, &nyworkDown, &ypad);
        // Think of decimating as filtering then downsampling. In this case
        // we would extract the filtered signal beginning at ypad[groupdelay].
        // However, because we've downsampled we instead extract at 
        // ypad[groupDelay/downFactor].  This is why we took care to pad the
        // number of FIR coefficients during setup.
        int ibeg = pImpl->mGroupDelay/pImpl->mDownFactor; 
        int ncopy = std::max(0, nyworkDown - ibeg);
        *nyDown = ncopy;
        if (ncopy != nyref)
        {
            fprintf(stderr, "ncopy = %d should match %d", ncopy, nyref); 
        }
        ippsCopy_64f(ypad, y, ncopy); 
        ippsFree(xpad);
        ippsFree(ypad);
    }
    else
    {
    } 
}

int Decimate::getDownsamplingFactor() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mDownFactor; 
}

int Decimate::getFIRFilterLength() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mFIRLength;
}
