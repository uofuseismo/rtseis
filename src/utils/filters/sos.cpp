#include <stdio.h>
#include <stdlib.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#define IPPS_CORE_SRC 1
#include "rtseis/utils/filters.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utils::Filters;

/*!
 * @defgroup rtseis_utils_filters_sos Second Order Sections
 * @brief This is the core implementation for second order section (biquad)
 *        infinite impulse response filtering.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
/*!
 * @brief Default constructor.
 * @ingroup rtseis_utils_filters_sos
 */
SOSFilter::SOSFilter(void)
{
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_utils_filters_sos
 */
SOSFilter::~SOSFilter(void)
{
    clear();
    return;
}
/*!
 * @brief Releases memory and resets the filter.
 * @ingroup rtseis_utils_filters_sos
 */
void SOSFilter::clear(void)
{
    if (pBuf_ != nullptr){ippsFree(pBuf_);}
    if (pTaps_ != nullptr){ippsFree(pTaps_);}
    if (dlysrc_ != nullptr){ippsFree(dlysrc_);}
    if (dlydst_ != nullptr){ippsFree(dlydst_);}
    if (bsRef_ != nullptr){ippsFree(bsRef_);}
    if (asRef_ != nullptr){ippsFree(asRef_);}
    if (zi_ != nullptr){ippsFree(zi_);}
    pState_ = nullptr;
    pBuf_ = nullptr;
    pTaps_ = nullptr;
    dlysrc_ = nullptr;
    dlydst_ = nullptr;
    bsRef_ = nullptr;
    asRef_ = nullptr;
    zi_ = nullptr;
    setPrecision(RTSeis::Precision::DOUBLE);
    toggleRealTime(false);
    nsections_ = 0;
    tapsLen_ = 0;
    nwork_ = 0;
    bufferSize_ = 0;
    linit_ = false;
    return;
}
/*!
 * @brief Copy constructor.
 * @param[in] sos  Class from which to initialize.
 * @ingroup rtseis_utils_filters_sos
 */
SOSFilter::SOSFilter(const SOSFilter &sos)
{
    *this = sos;
    return;
}
/*!
 * @brief Copy operator.
 * @param[in] sos   SOS filter class to copy.
 * @result A deep copy of the input class.
 * @ingroup rtseis_utils_filters_sos
 */
SOSFilter& SOSFilter::operator=(const SOSFilter &sos)
{
    if (&sos == this){return *this;}
    clear();
    if (!sos.linit_){return *this;}
    // Reinitialize the filter
    initialize(sos.nsections_, sos.bsRef_, sos.asRef_, 
               sos.isRealTime(), sos.getPrecision());
    // Now copy the filter states
    if (bufferSize_ > 0)
    {
        Ipp8u *pBufIn = static_cast<Ipp8u *> (sos.pBuf_);
        Ipp8u *pBufOut = static_cast<Ipp8u *> (pBuf_);
        ippsCopy_8u(pBufIn, pBufOut, bufferSize_);
    }
    // Copy the initial conditions
    if (nsections_ > 0){ippsCopy_64f(sos.zi_, zi_, 2*nsections_);}
    // And the delay lines
    if (nwork_ > 0)
    {
        if (isDoublePrecision())
        {
            Ipp64f *dlysrcIn  = static_cast<Ipp64f *> (sos.dlysrc_);
            Ipp64f *dlysrcOut = static_cast<Ipp64f *> (dlysrc_);
            ippsCopy_64f(dlysrcIn, dlysrcOut, nwork_); 
            Ipp64f *dlydstIn  = static_cast<Ipp64f *> (sos.dlydst_);
            Ipp64f *dlydstOut = static_cast<Ipp64f *> (dlydst_);
            ippsCopy_64f(dlydstIn, dlydstOut, nwork_);
        }
        else
        {
            Ipp32f *dlysrcIn  = static_cast<Ipp32f *> (sos.dlysrc_);
            Ipp32f *dlysrcOut = static_cast<Ipp32f *> (dlysrc_);
            ippsCopy_32f(dlysrcIn, dlysrcOut, nwork_); 
            Ipp32f *dlydstIn  = static_cast<Ipp32f *> (sos.dlydst_);
            Ipp32f *dlydstOut = static_cast<Ipp32f *> (dlydst_);
            ippsCopy_32f(dlydstIn, dlydstOut, nwork_);
        }
    }
    return *this;
}
/*!
 * @brief Initializes the second order section filter.
 * @param[in] ns           The number of second order sections.
 * @param[in] bs           Numerator coefficients.  This is an array of
 *                         dimension [3 x ns] with leading dimension 3.
 *                         There is a further requirement that b[3*is]
 *                         for \f$ i_s=0,1,\cdots,n_s-1 \f$ not be zero.
 * @param[in] as           Denominator coefficients.  This is an array of
 *                         dimension [3 x ns] with leading dimension 3. 
 *                         There is a further requirement that a[3*is]
 *                         for \f$ i_s=0,1,\cdots,n_s-1 \f$ not be zero.
 * @param[in] lisRealTime  Flag indicating that this is for real-time
 *                         application.
 * @param[in] precision    Determines the precision of the underlying
 *                         filter application.
 * @ingroup rtseis_utils_filters_sos
 */
int SOSFilter::initialize(const int ns,
                          const double bs[],
                          const double as[],
                          const bool lisRealTime,
                          const RTSeis::Precision precision)
{
    clear();
    // Checks
    if (ns < 1 || bs == NULL || as == NULL)
    {
        if (ns < 1){RTSEIS_ERRMSG("%s", "No sections\n");}
        if (bs == NULL){RTSEIS_ERRMSG("%s", "bs is NULL\n");}
        if (as == NULL){RTSEIS_ERRMSG("%s", "as is NULL\n");}
        return -1;
    }
    // Verify the highest order coefficients make sense
    for (int i=0; i<ns; i++)
    {
        if (bs[3*i] == 0.0)
        {
            RTSEIS_ERRMSG("Leading bs coefficient of section %d is zero", i);
            return -1;
        }
        if (as[3*i] == 0.0)
        {
            RTSEIS_ERRMSG("Leading as coefficient of section %d is zero", i);
            return -1;
        }
    }
    // Figure out sizes and copy the inputs
    nsections_ = ns;
    tapsLen_ = 6*nsections_;
    nwork_ = std::max(128, 2*nsections_);
    bsRef_ = ippsMalloc_64f(3*nsections_);
    ippsCopy_64f(bs, bsRef_, 3*nsections_);
    asRef_ = ippsMalloc_64f(3*nsections_);
    ippsCopy_64f(as, asRef_, 3*nsections_);
    zi_ = ippsMalloc_64f(2*nsections_);
    ippsZero_64f(zi_, 2*nsections_);
    IppStatus status;
    if (precision == RTSeis::Precision::DOUBLE)
    {
        status = ippsIIRGetStateSize_BiQuad_64f(nsections_, &bufferSize_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to get state size");
            clear();
            return -1;
        }
        IppsIIRState_64f *pState = nullptr; 
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        Ipp64f *pTaps = ippsMalloc_64f(tapsLen_);
        for (int i=0; i<nsections_; i++)
        {
            pTaps[6*i+0] = bs[3*i+0];
            pTaps[6*i+1] = bs[3*i+1];
            pTaps[6*i+2] = bs[3*i+2];
            pTaps[6*i+3] = as[3*i+0];
            pTaps[6*i+4] = as[3*i+1];
            pTaps[6*i+5] = as[3*i+2];
        }
        Ipp64f *dlysrc = ippsMalloc_64f(nwork_);
        ippsZero_64f(dlysrc, nwork_);
        Ipp64f *dlydst = ippsMalloc_64f(nwork_);
        ippsZero_64f(dlydst, nwork_);
        status = ippsIIRInit_BiQuad_64f(&pState, pTaps, nsections_,
                                        dlysrc, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to initialized biquad filter");
            clear();
            return -1;
        }
        pState_ = pState;
        pBuf_ = pBuf;
        pTaps_ = pTaps;
        dlysrc_ = dlysrc;
        dlydst_ = dlydst;
    }
    else
    {
        status = ippsIIRGetStateSize_BiQuad_32f(nsections_, &bufferSize_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to get state size");
            clear();
            return -1;
        }
        IppsIIRState_32f *pState = nullptr;
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        Ipp32f *pTaps = ippsMalloc_32f(tapsLen_);
        for (int i=0; i<nsections_; i++)
        {
            pTaps[6*i+0] = static_cast<Ipp32f> (bs[3*i+0]);
            pTaps[6*i+1] = static_cast<Ipp32f> (bs[3*i+1]);
            pTaps[6*i+2] = static_cast<Ipp32f> (bs[3*i+2]);
            pTaps[6*i+3] = static_cast<Ipp32f> (as[3*i+0]);
            pTaps[6*i+4] = static_cast<Ipp32f> (as[3*i+1]);
            pTaps[6*i+5] = static_cast<Ipp32f> (as[3*i+2]);
        }
        Ipp32f *dlysrc = ippsMalloc_32f(nwork_);
        ippsZero_32f(dlysrc, nwork_);
        Ipp32f *dlydst = ippsMalloc_32f(nwork_);
        ippsZero_32f(dlydst, nwork_);
        status = ippsIIRInit_BiQuad_32f(&pState, pTaps, nsections_,
                                        dlysrc, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to initialized biquad filter");
            clear();
            return -1;
        }
        pState_ = pState;
        pBuf_ = pBuf;
        pTaps_ = pTaps;
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
 * @param[in] nz   The second order section filter initial conditions.
 *                 This should be equal to getInitialConditionLength().
 * @param[in] zi   The initial conditions.  This has dimension [nz].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_sos
 */
int SOSFilter::setInitialConditions(const int nz, const double zi[])
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
 *        when SOSFilter::setInitialConditions() was called.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_sos
 */
int SOSFilter::resetInitialConditions(void)
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    if (isDoublePrecision())
    {
        Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
        ippsCopy_64f(zi_, dlysrc, 2*nsections_);
    }
    else
    {
        Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
        ippsConvert_64f32f(zi_, dlysrc, 2*nsections_);
    }
    return 0;
}
/*!
 * @brief Applies the second order section filter to the data.
 * @param[in] n   Number of points in signals.
 * @param[in] x   The signal to filter.  This has dimension [n].
 * @param[out] y  The filtered signal.  This has dimension [n].
 * @ingroup rtseis_utils_filters_sos
 */
int SOSFilter::apply(const int n, const double x[], double y[]) 
{
    if (n <= 0){return 0;}
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1;
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
    // Get a pointer to the filter state and set the initial conditions
    IppsIIRState_64f *pState = static_cast<IppsIIRState_64f *> (pState_);
    Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
    IppStatus status = ippsIIRSetDlyLine_64f(pState, dlysrc);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to set delay line");
        return -1;
    }
    // Apply the filters
    status = ippsIIR_64f(x, y, n, pState); 
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
    }
    if (isRealTime())
    {
        Ipp64f *dlydst = static_cast<Ipp64f *> (dlydst_);
        status = ippsIIRGetDlyLine_64f(pState, dlydst);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to get delay line");
            return -1;
        }
        ippsCopy_64f(dlydst, dlysrc, 2*nsections_);
    }
    return 0;
}
/*!
 * @brief Applies the second order section filter to the data.
 * @param[in] n   Number of points in signals.
 * @param[in] x   The signal to filter.  This has dimension [n].
 * @param[out] y  The filtered signal.  This has dimension [n].
 * @ingroup rtseis_utils_filters_sos
 */
int SOSFilter::apply(const int n, const float x[], float y[])
{
    if (n <= 0){return 0;}
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1;
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
    // Get a pointer to the filter state and set the initial conditions
    IppsIIRState_32f *pState = static_cast<IppsIIRState_32f *> (pState_);
    Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
    IppStatus status = ippsIIRSetDlyLine_32f(pState, dlysrc);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to set delay line");
        return -1;
    }
    // Apply the filters
    status = ippsIIR_32f(x, y, n, pState);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
    }
    if (isRealTime())
    {
        Ipp32f *dlydst = static_cast<Ipp32f *> (dlydst_);
        status = ippsIIRGetDlyLine_32f(pState, dlydst);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to get delay line");
            return -1;
        }
        ippsCopy_32f(dlydst, dlysrc, 2*nsections_);
    }
    return 0;
}
/*!
 * @brief Returns the length of the initial conditions.
 * @result The length of the initial condtions array.
 * @ingroup rtseis_utils_filters_sos
 */
int SOSFilter::getInitialConditionLength(void) const
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    return 2*nsections_;
}
/*!
 * @brief Returns the number of sections in the filter.
 * @result The number of cascaded second order sections.
 * @ingroup rtseis_utils_filters_sos
 */
int SOSFilter::getNumberOfSections(void) const
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    return nsections_;
}
