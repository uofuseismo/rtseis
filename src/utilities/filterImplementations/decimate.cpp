#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <climits>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/filterImplementations/decimate.hpp"
#include "rtseis/utilities/filterDesign/fir.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
//#include "rtseis/utilities/filterImplementations/multiRateFIRFilter.hpp"
#include "rtseis/utilities/filterImplementations/downsample.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"

using namespace RTSeis::Utilities;
using namespace RTSeis::Utilities::FilterImplementations;

class Decimate::DecimateImpl
{
public:
    //class MultiRateFIRFilter mMRFIRFilter; // TODO implementaiton is slow!
    class FIRFilter mFIRFilter;
    class Downsample<double> mDownsampler; 
    int mDownFactor = 1;
    int mGroupDelay = 0;
    int mFIRLength = 0;
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
    //pImpl->mMRFIRFilter.clear();
    pImpl->mFIRFilter.clear();
    pImpl->mDownsampler.clear();
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
            int groupDelay = (nfir-1)/2;
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
    auto r = 1.0/static_cast<double> (downFactor);
    auto fir = FilterDesign::FIR::FIR1Lowpass(order, r,
                                              FilterDesign::FIRWindow::HAMMING);
    // Set the multirate FIR filter
    auto b = fir.getFilterTaps();
    int ntaps = b.size();
    try
    {
        /*
        //constexpr int upFactor = 1;
        pImpl->mMRFIRFilter.initialize(upFactor, downFactor,
                                       ntaps, b.data(),
                                       mode, precision);
        */
        pImpl->mFIRFilter.initialize(ntaps, b.data(),
                                     mode,
                                     precision);
        pImpl->mDownsampler.initialize(downFactor,
                                       mode);
    }
    catch (std::exception &e)
    {
        RTSEIS_THROW_RTE("%sFIR filter initialization failed",
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
    return pImpl->mDownsampler.estimateSpace(n);
}
/* TODO - when lashing in a more performant multirate fir filter use this fn
int Decimate::estimateSpace(const int n) const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (n < 0){RTSEIS_THROW_IA("n=%d cannot be negative", n);}
    if (pImpl->mRemovePhaseShift)
    {
        int npad = n + pImpl->mGroupDelay;
        int len = pImpl->mMRFIRFilter.estimateSpace(npad);
        int i0 = pImpl->mGroupDelay/pImpl->mDownFactor; 
        double d1 = std::ceil(static_cast<double> (npad)
                             /static_cast<double> (pImpl->mDownFactor)); 
        int i1 = static_cast<int> (d1);
        i1 = std::max(i0, std::min(len - 1 - 1, i1 - 1));
        return i1 - i0 + 1;
    }
    else
    {
        return pImpl->mMRFIRFilter.estimateSpace(n);
    }
}
*/

int Decimate::getInitialConditionLength() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mFIRFilter.getInitialConditionLength();
    //return pImpl->mMRFIRFilter.getInitialConditionLength(); 
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
    //pImpl->mMRFIRFilter.setInitialConditions(nz, zi);
    pImpl->mFIRFilter.setInitialConditions(nz, zi);
    pImpl->mDownsampler.resetInitialConditions();
}

void Decimate::resetInitialConditions()
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    //pImpl->mMRFIRFilter.resetInitialConditions();
    pImpl->mFIRFilter.resetInitialConditions();
    pImpl->mDownsampler.resetInitialConditions();
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
#ifdef DEBUG
        assert(pImpl->mGroupDelay%pImpl->mDownFactor == 0); 
#endif
        // Postpend zeros so that the filter delay pushes our desired
        // output to the end of the temporary array
        int npad = nx + pImpl->mGroupDelay; 
        double *xpad = ippsMalloc_64f(npad);
        ippsCopy_64f(x, xpad, nx);
        ippsZero_64f(&xpad[nx], pImpl->mGroupDelay);
        // The FIR filter seems to perform better than the multirate FIR
        // filter.  So don't fight it.
        double *yfilt = ippsMalloc_64f(npad);
        pImpl->mFIRFilter.apply(npad, xpad, &yfilt);
        pImpl->mDownsampler.apply(nx, &yfilt[pImpl->mGroupDelay],  
                                  ny, nyDown, &y);
        ippsFree(xpad);
        ippsFree(yfilt);
    }
    else
    {
        // Just apply the downsampler
        //pImpl->mMRFIRFilter.apply(nx, x, ny, nyDown, &y);
        double *yfilt = ippsMalloc_64f(nx);
        pImpl->mFIRFilter.apply(nx, x, &yfilt);
        pImpl->mDownsampler.apply(nx, yfilt, ny, nyDown, &y); 
        ippsFree(yfilt);
    }
}

/*
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
#ifdef DEBUG
        assert(pImpl->mGroupDelay%pImpl->mDownFactor == 0);
#endif
        // Postpend zeros so that the filter delay pushes our desired
        // output to the end of the temporary array
        int npad = nx + pImpl->mGroupDelay; 
        double *xpad = ippsMalloc_64f(npad);
        ippsCopy_64f(x, xpad, nx);
        ippsZero_64f(&xpad[nx], pImpl->mGroupDelay); 

        int nywork = pImpl->mMRFIRFilter.estimateSpace(npad);
        int nyworkDown;
        double *ypad = ippsMalloc_64f(nywork);
        pImpl->mMRFIRFilter.apply(npad, xpad, nywork, &nyworkDown, &ypad); 
        // Recall, that to remove the phase shift I would nominally start at
        // the groupDelay'th index.  However, I've downsampled the signal so
        // I instead start at the groupDelay/downFactor'th sample.
        int i0 = pImpl->mGroupDelay/pImpl->mDownFactor; 
        // The filter will delay nx + groupDelay samples.  However, I've only
        // retained every (nx + groupDelay)/downFactor samples so this should
        // correspond to the last `useful' or `reliable' sample.  I've noticed
        // that sometimes IPP will compute the i1'th sample however,
        // sometimes it doesn't, so this is relatively safe for safety 
        // I use nyworkDown - 1.  The additional -1 is to make this inclusive. 
        double d1 = std::ceil(static_cast<double> (npad)
                             /static_cast<double> (pImpl->mDownFactor)); 
        int i1 = static_cast<int> (d1);
        i1 = std::max(i0, std::min(nyworkDown - 1 - 1, i1 - 1));
        // Copy
        *nyDown = i1 - i0 + 1;
        ippsCopy_64f(&ypad[i0], y, *nyDown); 
        ippsFree(xpad);
        ippsFree(ypad);
    }
    else
    {
        // Just apply the downsampler
        pImpl->mMRFIRFilter.apply(nx, x, ny, nyDown, &y);
    } 
}
*/

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
