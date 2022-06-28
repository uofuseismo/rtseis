#include <stdexcept>
#include <string>
#include <climits>
#include <string>
#ifndef NDEBBUG
#include <cassert>
#endif
#include <ipps.h>
#include "rtseis/enums.hpp"
#include "rtseis/filterImplementations/decimate.hpp"
#include "rtseis/filterDesign/fir.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include "rtseis/filterImplementations/downsample.hpp"
#include "rtseis/filterImplementations/firFilter.hpp"

using namespace RTSeis::FilterImplementations;

template<RTSeis::ProcessingMode E, class T>
class Decimate<E, T>::DecimateImpl
{
public:
    void initialize(const int downFactor,
                    const int filterLength,
                    const bool lRemovePhaseShift)
    {
        mDownFactor = downFactor;
        int nfir = filterLength;
        // Postprocessing is a little trickier - may have to extend filter length
        if (mMode == RTSeis::ProcessingMode::POST_PROCESSING &&
            lRemovePhaseShift)
        {
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
#ifndef NDEBUG
            assert(!lfail);
#endif
            if (lfail)
            {
                throw std::runtime_error("Padding algorithmic failure");
            }
            mGroupDelay = nfir/2;
            mRemovePhaseShift = true;
        }
        // Create a hamming filter.
        mFIRLength = nfir;
        int order = nfir - 1;
        auto r = 1.0/static_cast<double> (downFactor);
        auto fir = FilterDesign::FIR::FIR1Lowpass(order, r,
                                                  FilterDesign::FIRWindow::HAMMING);
        // Set the multirate FIR filter
        auto b = fir.getFilterTaps();
        int ntaps = b.size();
        try
        {
            mFIRFilter.initialize(ntaps, b.data());
            mDownsampler.initialize(downFactor);
        }
        catch (std::exception &e)
        {
            auto errmsg = std::string("FIR filter initialization failed with: ")
                        + e.what();
            throw std::runtime_error(errmsg);
        }
        mInitialized = true;
    }
    void apply(const int nx, const T x[],
               const int ny, int *nyDown, T y[])
    {
        if (mRemovePhaseShift)
        {
#ifndef NDEBUG
            assert(mGroupDelay%mDownFactor == 0);
#endif
            // Postpend zeros so that the filter delay pushes our desired
            // output to the end of the temporary array
            int npad = nx + mGroupDelay;
            auto xpad = reinterpret_cast<T *> (ippsMalloc_8u(npad*sizeof(T)));
            std::copy(x, x + nx, xpad);
            std::fill(xpad + nx, xpad + nx + mGroupDelay, 0);
            //ippsCopy_64f(x, xpad, nx);
            //ippsZero_64f(&xpad[nx], mGroupDelay);
            // The FIR filter seems to perform better than the multirate FIR
            // filter.  So don't fight it.
            //double *yfilt = ippsMalloc_64f(npad);
            auto yfilt = reinterpret_cast<T *> (ippsMalloc_8u(npad*sizeof(T)));
            mFIRFilter.apply(npad, xpad, &yfilt);
            mDownsampler.apply(nx, &yfilt[mGroupDelay],
                               ny, nyDown, &y);
            ippsFree(xpad);
            ippsFree(yfilt);
        }
        else
        {
            // Just apply the downsampler
            //pImpl->mMRFIRFilter.apply(nx, x, ny, nyDown, &y);
            //double *yfilt = ippsMalloc_64f(nx);
            auto yfilt = reinterpret_cast<T *> (ippsMalloc_8u(nx*sizeof(T)));
            mFIRFilter.apply(nx, x, &yfilt);
            mDownsampler.apply(nx, yfilt, ny, nyDown, &y);
            ippsFree(yfilt);
        }
    }
    //class MultiRateFIRFilter mMRFIRFilter; // TODO implementation is slow!
    class FIRFilter<E, T> mFIRFilter;
    class Downsample<E, T> mDownsampler; 
    int mDownFactor = 1;
    int mGroupDelay = 0;
    int mFIRLength = 0;
    const RTSeis::ProcessingMode mMode = E;
    //RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
    bool mRemovePhaseShift = false;
    bool mInitialized = false;
};

/// C'tor
template<RTSeis::ProcessingMode E, class T>
Decimate<E, T>::Decimate() :
    pImpl(std::make_unique<DecimateImpl> ())
{
}

/// Copy c'tor
template<RTSeis::ProcessingMode E, class T>
Decimate<E, T>::Decimate(const Decimate &decimate)
{
    *this = decimate;
}

/// Move c'tor
template<RTSeis::ProcessingMode E, class T>
Decimate<E, T>::Decimate(Decimate &&decimate) noexcept
{
   *this = std::move(decimate);
}

/// Copy assignment
template<RTSeis::ProcessingMode E, class T>
Decimate<E, T>& Decimate<E, T>::operator=(const Decimate &decimate)
{
    if (&decimate == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<DecimateImpl> (*decimate.pImpl);
    return *this;
}

/// Move assignemnt
template<RTSeis::ProcessingMode E, class T>
Decimate<E, T>& Decimate<E, T>::operator=(Decimate &&decimate) noexcept
{
    if (&decimate == this){return *this;}
    pImpl = std::move(decimate.pImpl);
    return *this;
}

/// Destructor
template<RTSeis::ProcessingMode E, class T>
Decimate<E, T>::~Decimate() = default;

/// Clears the class
template<RTSeis::ProcessingMode E, class T>
void Decimate<E, T>::clear() noexcept
{
    //pImpl->mMRFIRFilter.clear();
    pImpl->mFIRFilter.clear();
    pImpl->mDownsampler.clear();
    //pImpl->mMode = RTSeis::ProcessingMode::POST_PROCESSING;
    //pImpl->mPrecision = RTSeis::Precision::DOUBLE;
    pImpl->mDownFactor = 1;
    pImpl->mGroupDelay = 0;
    pImpl->mFIRLength = 0;
    pImpl->mRemovePhaseShift = false;
    pImpl->mInitialized = false;
}

/// Initilalization
template<RTSeis::ProcessingMode E, class T>
void Decimate<E, T>::initialize(const int downFactor,
                                const int filterLength,
                                const bool lRemovePhaseShift)
{
    clear();
    if (downFactor < 2)
    {
        auto errmsg = "Downsampling factor = " + std::to_string(downFactor)
                    + " must be at least 2";
        throw std::invalid_argument(errmsg);
    }
    if (filterLength < 5)
    {
        auto errmsg = "Filter length = " + std::to_string(filterLength)
            + " must be at least 5";
    }
    pImpl->initialize(downFactor, filterLength, lRemovePhaseShift);
    /*
    // Set some properties
    pImpl->mDownFactor = downFactor;
    //pImpl->mPrecision = precision;
    //pImpl->mMode = mode;
    int nfir = filterLength;
    // Postprocessing is a little trickier - may have to extend filter length
    if (pImpl->mode == RTSeis::ProcessingMode::POST_PROCESSING &&
        lremovePhaseShift)
    {
        // N.B.  The idea is that we think of decimation as filtering then
        // downsampling.  The filter will introduce a group delay of
        // groupDelay = nfir/2 where nfir is odd (i.e., a unit impulse
        // introduce 1/2 = 0 delay).  The delay is then removed by
        // extracting the filtered signal at yfiltered[groupDelay].
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
#ifndef NDEBUG
        assert(!lfail);
#endif
        if (lfail){throw std::runtime_error("Algorithmic failure");}
        pImpl->mGroupDelay = nfir/2;
        pImpl->mRemovePhaseShift = true;
    } 
    // Create a hamming filter.  Note, downFactor = 1 can trigger an error.
    pImpl->mFIRLength = nfir;
    int order = nfir - 1; 
    auto r = 1.0/static_cast<double> (downFactor);
    auto fir = FilterDesign::FIR::FIR1Lowpass(order, r,
                                              FilterDesign::FIRWindow::HAMMING);
    // Set the FIR filter
    auto b = fir.getFilterTaps();
    int ntaps = b.size();
    try
    {
        //constexpr int upFactor = 1;
        ////pImpl->mMRFIRFilter.initialize(upFactor, downFactor,
        ////                               ntaps, b.data(),
        ////                               mode, precision);
        pImpl->mFIRFilter.initialize(ntaps, b.data());
        pImpl->mDownsampler.initialize(downFactor,
                                       mode);
    }
    catch (std::exception &e)
    {
        auto errmsg = "FIR filter initialization failed with: "
                    + std::to_string(e.what());
        throw std::runtime_error(errmsg);
    }
    pImpl->mInitialized = true;
     */
}

/*
template<RTSeis::ProcessingMode E, class T>
void Decimate<E, T>::initialize(const int downFactor,
                                const int filterLength,
                                const bool lRemovePhaseShift)
{
    //constexpr RTSeis::Precision precision = RTSeis::Precision::FLOAT;
    clear();
    if (downFactor < 2)
    {
        auto errmsg = "Downsampling factor = " + std::to_string(downFactor)
                    + " must be at least 2";
        throw std::invalid_argument(errmsg);
    }
    if (filterLength < 5)
    {
        auto errmsg = "Filter length = " + std::to_string(filterLength)
                    + " must be at least 5";
    }
    // Set some properties
    pImpl->mDownFactor = downFactor;
    //pImpl->mPrecision = precision;
    //pImpl->mMode = mode;
    int nfir = filterLength;
    // Postprocessing is a little trickier - may have to extend filter length
    if (pImpl->mMode == RTSeis::ProcessingMode::POST_PROCESSING &&
        lRemovePhaseShift)
    {
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
#ifndef NDEBUG
        assert(!lfail);
#endif
        if (lfail){throw std::runtime_error("Algorithmic failure");}
        pImpl->mGroupDelay = nfir/2;
        pImpl->mRemovePhaseShift = true;
    }
    // Create a hamming filter.
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
        //constexpr int upFactor = 1;
        ////pImpl->mMRFIRFilter.initialize(upFactor, downFactor,
        ////                               ntaps, b.data());
        pImpl->mFIRFilter.initialize(ntaps, b.data());
        pImpl->mDownsampler.initialize(downFactor,
                                       mode);
    }
    catch (std::exception &e)
    {
        auto errmsg = "FIR filter initialization failed with: "
            + std::to_string(e.what());
        throw std::runtime_error(errmsg);
    }
    pImpl->mInitialized = true;
}
*/

/// Initialized?
template<RTSeis::ProcessingMode E, class T>
bool Decimate<E, T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Estimate output space
template<RTSeis::ProcessingMode E, class T>
int Decimate<E, T>::estimateSpace(const int n) const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (n < 0){throw std::invalid_argument("n cannot be negative");}
    return pImpl->mDownsampler.estimateSpace(n);
}
/* TODO - when lashing in a more performant multirate fir filter use this fn
int Decimate::estimateSpace(const int n) const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (n < 0){throw std::invalid_argument("n cannot be negative");}
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

/// Get initial conditions length
template<RTSeis::ProcessingMode E, class T>
int Decimate<E, T>::getInitialConditionLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mFIRFilter.getInitialConditionLength();
}

/// Set initial conditions
template<RTSeis::ProcessingMode E, class T>
void Decimate<E, T>::setInitialConditions(const int nz,const double zi[])
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int nzref = getInitialConditionLength();
    if (nz != nzref)
    {
        throw std::invalid_argument("nz = " + std::to_string(nz)
                                  + " must be equal to "
                                  + std::to_string(nzref));
    }
    if (nz > 0 && zi == nullptr){throw std::invalid_argument("zi is NULL");}
    //pImpl->mMRFIRFilter.setInitialConditions(nz, zi);
    pImpl->mFIRFilter.setInitialConditions(nz, zi);
    pImpl->mDownsampler.resetInitialConditions();
}

/// Reset initial conditions
template<RTSeis::ProcessingMode E, class T>
void Decimate<E, T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    //pImpl->mMRFIRFilter.resetInitialConditions();
    pImpl->mFIRFilter.resetInitialConditions();
    pImpl->mDownsampler.resetInitialConditions();
}

/// Apply decimator (double)
template<RTSeis::ProcessingMode E, class T>
void Decimate<E, T>::apply(const int nx, const T x[],
                           const int ny, int *nyDown, T *yIn[])
{
    *nyDown = 0;
    if (nx <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    int nyref = estimateSpace(nx);
    if (ny < nyref)
    {
        auto errmsg = "ny = " + std::to_string(ny) + " must be at least "
                    + std::to_string(nyref);
        throw std::invalid_argument(errmsg);
    }
    auto y = *yIn;
    if (y == nullptr){throw std::invalid_argument("y is NULL");}
    pImpl->apply(nx, x, ny, nyDown, y);
    /*
    // Special case
    if (pImpl->mRemovePhaseShift)
    {   
#ifndef NDEBUG
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
    */
}

/*
/// Apply decimator (float)
template<RTSeis::ProcessingMode E>
void Decimate<E, float>::apply(const int nx, const float x[],
                               const int ny, int *nyDown, float *yIn[])
{
    *nyDown = 0;
    if (nx <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (x == nullptr){throw std::invalid_argumnet("x is NULL");}
    int nyref = estimateSpace(nx);
    if (ny < nyref){throw std::invalid_argument("ny = " + std:::to_string(ny) + " must be at least " + std::to_string(nyref));}
    float *y = *yIn;
    if (y == nullptr)
    {
        throw std::invalid_argument("y is NULL");
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
        float *xpad = ippsMalloc_32f(npad);
        ippsCopy_32f(x, xpad, nx);
        ippsZero_32f(&xpad[nx], pImpl->mGroupDelay);
        // The FIR filter seems to perform better than the multirate FIR
        // filter.  So don't fight it.
        float *yfilt = ippsMalloc_32f(npad);
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
        float *yfilt = ippsMalloc_32f(nx);
        pImpl->mFIRFilter.apply(nx, x, &yfilt);
        pImpl->mDownsampler.apply(nx, yfilt, ny, nyDown, &y);
        ippsFree(yfilt);
    }
}
*/

/*
void Decimate::apply(const int nx, const double x[],
                     const int ny, int *nyDown, double *yIn[])
{
    *nyDown = 0;
    if (nx <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    int nyref = estimateSpace(nx);
    if (ny < nyref){throw std::invalid_argument("ny = " + std::to_string(ny) + " must be at least " + std::to_string(nyref));}
    double *y = *yIn;
    if (y == nullptr)
    {
        throw std::invalid_argument("y is NULL");
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

/// Get downsampling factor
template<RTSeis::ProcessingMode E, class T>
int Decimate<E, T>::getDownsamplingFactor() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mDownFactor; 
}

/// Get FIR filter length
template<RTSeis::ProcessingMode E, class T>
int Decimate<E, T>::getFIRFilterLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mFIRLength;
}

///--------------------------------------------------------------------------///
///                         Template instantiation                           ///
///--------------------------------------------------------------------------///
template class RTSeis::FilterImplementations::Decimate<RTSeis::ProcessingMode::POST, double>;
template class RTSeis::FilterImplementations::Decimate<RTSeis::ProcessingMode::REAL_TIME, double>;
template class RTSeis::FilterImplementations::Decimate<RTSeis::ProcessingMode::POST, float>;
template class RTSeis::FilterImplementations::Decimate<RTSeis::ProcessingMode::REAL_TIME, float>;
