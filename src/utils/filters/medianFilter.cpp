#include <stdio.h>
#include <stdlib.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#define IPPS_CORE_SRC 1
#include "rtseis/utils/filters.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utils::Filters;

/*!
 * @defgroup rtseis_utils_filters_median Median Filter
 * @brief This is the core implementation for median filtering.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
/*!
 * @brief Default constructor.
 * @ingroup rtseis_utils_filters_median
 */
MedianFilter::MedianFilter(void)
{
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_utils_filters_median
 */
MedianFilter::~MedianFilter(void)
{
    clear();
    return;
}
/*!
 * @brief Copy constructor.
 * @param[in] median  Class from which to initialize.
 * @ingroup rtseis_utils_filters_median
 */
MedianFilter::MedianFilter(const MedianFilter &median)
{
    *this = median;
    return;
}
/*!
 * @brief Copy operator.
 * @param[in] median   Median filter class to copy.
 * @result A deep copy of the input class.
 * @ingroup rtseis_utils_filters_median
 */
MedianFilter& MedianFilter::operator=(const MedianFilter &median)
{
    if (&median == this){return *this;}
    clear();
    if (!median.linit_){return *this;}
    // Reinitialize the filter
    initialize(median.maskSize_, median.isRealTime(), median.getPrecision());
    // Now copy the filter states
    if (bufferSize_ > 0)
    {
        Ipp8u *pBufIn = static_cast<Ipp8u *> (median.pBuf_);
        Ipp8u *pBufOut = static_cast<Ipp8u *> (pBuf_);
        ippsCopy_8u(pBufIn, pBufOut, bufferSize_);
    }
    if (nwork_ > 0)
    {
        ippsCopy_64f(median.zi_, zi_, nwork_);
        if (isDoublePrecision())
        {
            Ipp64f *dlysrcIn  = static_cast<Ipp64f *> (median.dlysrc_);
            Ipp64f *dlysrcOut = static_cast<Ipp64f *> (dlysrc_);
            ippsCopy_64f(dlysrcIn, dlysrcOut, nwork_);
            Ipp64f *dlydstIn  = static_cast<Ipp64f *> (median.dlydst_);
            Ipp64f *dlydstOut = static_cast<Ipp64f *> (dlydst_);
            ippsCopy_64f(dlydstIn, dlydstOut, nwork_);
        }
        else
        {
            Ipp32f *dlysrcIn  = static_cast<Ipp32f *> (median.dlysrc_);
            Ipp32f *dlysrcOut = static_cast<Ipp32f *> (dlysrc_);
            ippsCopy_32f(dlysrcIn, dlysrcOut, nwork_); 
            Ipp32f *dlydstIn  = static_cast<Ipp32f *> (median.dlydst_);
            Ipp32f *dlydstOut = static_cast<Ipp32f *> (dlydst_);
            ippsCopy_32f(dlydstIn, dlydstOut, nwork_);
        }
    }
    return *this;
}
/*!
 * @brief Clears the module and resets all parameters.
 * @ingroup rtseis_utils_filters_median
 */
void MedianFilter::clear(void)
{
    if (pBuf_ != nullptr){ippsFree(pBuf_);}
    if (dlysrc_ != nullptr){ippsFree(dlysrc_);}
    if (dlydst_ != nullptr){ippsFree(dlydst_);}
    if (zi_ != nullptr){ippsFree(zi_);}
    setPrecision(RTSEIS_DOUBLE);
    toggleRealTime(false);
    pBuf_ = nullptr;
    dlysrc_ = nullptr;
    dlydst_ = nullptr;
    zi_ = nullptr;
    maskSize_ = 0;
    bufferSize_ = 0;
    linit_ = false;
    return;
}
/*!
 * @brief Initializes the median filter.
 * @param[in] n            The window size of the median filter.  This must 
 *                         be a positive and odd number.  If n is not add then
 *                         it's length will be increased by 1.
 * @param[in] lisRealTime  Flag indicating that this is for real-time
 *                         application.
 * @param[in] precision    Determines the precision of the underlying
 *                         median calculation. 
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_median
 */
int MedianFilter::initialize(const int n,
                             const bool lisRealTime,
                             const enum rtseisPrecision_enum precision)
{
    clear();
    // Set the mask size
    if (n < 1)
    {
        RTSEIS_ERRMSG("Mask size=%d must be postive", n);
        return -1;
    }
    maskSize_ = n;
    if (maskSize_%2 == 0)
    {
        maskSize_ = maskSize_ + 1;
        RTSEIS_WARNMSG("n=%d should be odd; setting to maskSize=%d",
                       n, maskSize_);
    }
    // Set the space
    nwork_ = std::max(8, maskSize_ - 1);
    zi_ = ippsMalloc_64f(nwork_);
    ippsZero_64f(zi_, nwork_);
    if (precision == RTSEIS_DOUBLE)
    {
        IppStatus status = ippsFilterMedianGetBufferSize(maskSize_, ipp64f,
                                                         &bufferSize_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting buffer size");
            return -1; 
        }
        Ipp64f *dlysrc = ippsMalloc_64f(nwork_);
        Ipp64f *dlydst = ippsMalloc_64f(nwork_);
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        ippsZero_64f(dlysrc, nwork_);
        ippsZero_64f(dlydst, nwork_);
        pBuf_   = pBuf;
        dlysrc_ = dlysrc;
        dlydst_ = dlydst;
    }
    else
    {
        IppStatus status = ippsFilterMedianGetBufferSize(maskSize_, ipp32f,
                                                         &bufferSize_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting buffer size");
            return -1;
        }
        Ipp32f *dlysrc = ippsMalloc_32f(nwork_);
        Ipp32f *dlydst = ippsMalloc_32f(nwork_);
        ippsZero_32f(dlysrc, nwork_);
        ippsZero_32f(dlydst, nwork_);
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        pBuf_   = pBuf;
        dlysrc_ = dlysrc;
        dlydst_ = dlydst;
    }
    setPrecision(precision);
    toggleRealTime(lisRealTime);
    linit_ = true;
    return 0;
}
/*!
 * @brief Sets the initial conditions for the filter.  This should be called
 *        prior to filter application as it will reset the filter.
 * @param[in] nz   The median filter initial conditions.  This should be
 *                 equal to getInitialConditionLength().
 * @param[in] zi   The initial conditions.  This has dimension [nz].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_median
 */
int MedianFilter::setInitialConditions(const int nz, const double zi[])
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    resetInitialConditions();
    int nzRef = getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef){RTSEIS_ERRMSG("nz=%d should equal %d", nz, nzRef);}
        if (zi == nullptr){RTSEIS_ERRMSG("%s", "zi is NULL");}
        return -1;
    }
    ippsCopy_64f(zi, zi_, nzRef);
    if (isDoublePrecision())
    {
        Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
        ippsCopy_64f(zi, dlysrc, nzRef);
    }
    else
    {
        Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
        ippsConvert_64f32f(zi, dlysrc, nzRef);
    }
    return 0;
}
/*!
 * @brief Resets the initial conditions on the source delay line to the
 *        default initial conditions or the initial conditions set 
 *        when MedianFilter::setInitialConditions() was called.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_median
 */
int MedianFilter::resetInitialConditions(void)
{
    if (!linit_)
    {   
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1; 
    }
    if (isDoublePrecision())
    {
        Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
        ippsCopy_64f(zi_, dlysrc, nwork_);
    }
    else
    {
        Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
        ippsConvert_64f32f(zi_, dlysrc, nwork_);
    }
    return 0;
}
/*!
 * @brief Utility routine to determine the initial condition length.
 * @retval A non-negative number is the length of the initial condition
 *         vector.
 * @retval 1 Indicates failure.
 * rtseis_utils_filters_median
 */
int MedianFilter::getInitialConditionLength(void) const
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    int len = maskSize_ - 1;
    return len;
}
/*!
 * @brief Appplies the median filter to the array x.
 * @param[in] n   Number of points in x.
 * @param[in] x   The signal to filter.  This has dimension [n].
 * @param[out] y  The filtered signal.  This has dimension [n].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_median
 */
int MedianFilter::apply(const int n, const double x[], double y[])
{
    if (n <= 0){return 0;} // Nothing to do
    if (!linit_)
    {   
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        RTSEIS_ERRMSG("%s", "x is NULL");
        RTSEIS_ERRMSG("%s", "y is NULL");
    }
    if (isFloatPrecision())
    {
        Ipp32f *x32 = ippsMalloc_32f(n);
        Ipp32f *y32 = ippsMalloc_32f(n);
        ippsConvert_64f32f(x, x32, n);
        int ierr = apply(n, x32, y32);
        ippsFree(x32);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to apply filter");
            ippsFree(y32);
            return -1;
        }
        ippsConvert_32f64f(y32, y, n);
        ippsFree(y32);
        return 0;
    }
    IppStatus status;
    Ipp8u *pBuf = static_cast<Ipp8u *> (pBuf_);
    Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
    Ipp64f *dlydst = nullptr;
    if (isRealTime()){dlydst = static_cast<Ipp64f *> (dlydst_);}
    status = ippsFilterMedian_64f(x, y, n, maskSize_,
                                  dlysrc, dlydst, pBuf);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Error applying filter");
        return -1;
    }
    if (isRealTime() && maskSize_ > 1)
    {
        ippsCopy_64f(dlydst, dlysrc, maskSize_-1);
    }
    return 0;
}
/*!
 * @brief Appplies the median filter to the array x.
 * @param[in] n   Number of points in x.
 * @param[in] x   The signal to filter.  This has dimension [n].
 * @param[out] y  The filtered signal.  This has dimension [n].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_median
 */
int MedianFilter::apply(const int n, const float x[], float y[])
{
    if (n <= 0){return 0;} // Nothing to do
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        RTSEIS_ERRMSG("%s", "x is NULL");
        RTSEIS_ERRMSG("%s", "y is NULL");
    }
    if (isDoublePrecision())
    {
        Ipp64f *x64 = ippsMalloc_64f(n);
        Ipp64f *y64 = ippsMalloc_64f(n);
        ippsConvert_32f64f(x, x64, n); 
        int ierr = apply(n, x64, y64);
        ippsFree(x64);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to apply filter");
            ippsFree(y64);
            return -1; 
        }
        ippsConvert_64f32f(y64, y, n); 
        ippsFree(y64);
        return 0;
    }
    IppStatus status;
    Ipp8u *pBuf = static_cast<Ipp8u *> (pBuf_);
    Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
    Ipp32f *dlydst = nullptr;
    if (isRealTime()){dlydst = static_cast<Ipp32f *> (dlydst_);}
    status = ippsFilterMedian_32f(x, y, n, maskSize_,
                                  dlysrc, dlydst, pBuf);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Error applying filter");
        return -1; 
    }
    if (isRealTime() && maskSize_ > 1)
    {
        ippsCopy_32f(dlydst, dlysrc, maskSize_-1);
    }
    return 0;
}
/*!
 * @brief Returns the group delay of the filter.  Note, that this shift is
 *        required to get a correspondence to Matlab.
 * @result The group delay.
 * @ingroup rtseis_utils_filters_median
 */
int MedianFilter::getGroupDelay(void) const
{
    int grpDelay = maskSize_/2;
    return grpDelay;
}
