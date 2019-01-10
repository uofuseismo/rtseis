#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/ippsHelper.h"
#include <ipps.h>

/*!
 * @defgroup rtseis_ipps_median Median Filtering
 * @brief Interface to the IPP median filter.
 * @ingroup rtseis_ipps
 * @author Ben Baker
 * @copyright ISTI distributed under the Apache 2 license.
 */
/*!
 * @brief Initializes the median filter for IPP.
 * @param[in] n            The size of the window when computing the
 *                         median.  This must be odd and positive.
 * @param[in] lisRealTime  If true then the filter is for real-time
 *                         applications.
 * @param[in] lisRealTime  Otherwise, it is for post-processing.
 * @param[in] precision    Precision of the median filter. 
 * @param[out] ippsMedian  Median filter structure to be used by IPP.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_median
 */
int rtseis_ippsMedianFilter_initialize(
    const int n,
    const bool lisRealTime,
    const enum rtseisPrecision_enum precision,
    struct ippsMedianFilter_struct *ippsMedian)
{
    IppStatus status;
    int maskSize;
    memset(ippsMedian, 0, sizeof(struct ippsMedianFilter_struct));
    maskSize = n;
    if (maskSize < 1)
    {
        RTSEIS_ERRMSG("n=%d must be positive", maskSize);
        return -1;
    }
    if (maskSize%2 == 0)
    {
        maskSize = maskSize + 1;
        RTSEIS_WARNMSG("n=%d should be odd; setting to maskSize=%d", n, maskSize);
    }
    if (precision == RTSEIS_DOUBLE)
    {
        status = ippsFilterMedianGetBufferSize(maskSize, ipp64f,
                                               &ippsMedian->bufferSize);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting buffer size");
            return -1;
        }
        ippsMedian->pBuf = ippsMalloc_8u(ippsMedian->bufferSize);
        ippsMedian->zi     = ippsMalloc_64f(MAX(1, maskSize-1));
        ippsMedian->dlysrc = ippsMalloc_64f(MAX(1, maskSize-1));
        ippsMedian->dlydst = ippsMalloc_64f(MAX(1, maskSize-1));
        ippsZero_64f(ippsMedian->zi,     MAX(1, maskSize-1));
        ippsZero_64f(ippsMedian->dlysrc, MAX(1, maskSize-1));
        ippsZero_64f(ippsMedian->dlydst, MAX(1, maskSize-1));
    }
    else if (precision == RTSEIS_DOUBLE)
    {
        status = ippsFilterMedianGetBufferSize(maskSize, ipp32f,
                                               &ippsMedian->bufferSize);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting buffer size");
            return -1;
        }
        ippsMedian->pBuf = (void *) ippsMalloc_8u(ippsMedian->bufferSize);
        ippsMedian->zi     = ippsMalloc_32f(MAX(1, maskSize-1));
        ippsMedian->dlysrc = ippsMalloc_32f(MAX(1, maskSize-1));
        ippsMedian->dlydst = ippsMalloc_32f(MAX(1, maskSize-1));
        ippsZero_32f(ippsMedian->zi,     MAX(1, maskSize-1));
        ippsZero_32f(ippsMedian->dlysrc, MAX(1, maskSize-1));
        ippsZero_32f(ippsMedian->dlydst, MAX(1, maskSize-1));
    }
    else
    {
        RTSEIS_ERRMSG("Invalid precision %d", (int) precision);
        return -1;
    }
    ippsMedian->precision = precision;
    ippsMedian->maskSize = maskSize;
    rtseis_ippsMedianFilter_toggleRealTime(lisRealTime, ippsMedian);
    ippsMedian->linit = true;
    return 0;
}
//============================================================================//
/*!
 * @brief Releases the memory on the median filter and sets all variables to
 *        zero or NULL.
 * @param[in,out] ippsMedian   On input contains the initialized median
 *                             filtering structure.
 * @param[in,out] ippsMedian   On exit all memory has been released and
 *                             variables set to 0 or NULL.
 * @ingroup rtseis_ipps_median
 */
int rtseis_ippsMedianFilter_finalize(
    struct ippsMedianFilter_struct *ippsMedian)
{
    if (!ippsMedian->linit){return 0;}
    ippsFree(ippsMedian->pBuf);
    ippsFree(ippsMedian->zi);
    ippsFree(ippsMedian->dlysrc);
    ippsFree(ippsMedian->dlydst);
    memset(ippsMedian, 0, sizeof(struct ippsMedianFilter_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Copies an initialized median filter structure, ippsMedianIn, to
 *        ippsMedianOut.  Note ,that this will also copy the initial conditions
 *        and delay lines.
 */
int rtseis_ippsMedianFilter_copy(
    const struct ippsMedianFilter_struct ippsMedianIn,
    struct ippsMedianFilter_struct *ippsMedianOut)
{
    int ierr;
    memset(ippsMedianOut, 0, sizeof(struct ippsMedianFilter_struct));
    if (!ippsMedianIn.linit)
    {
        RTSEIS_ERRMSG("%s", "Input median filter not yet set");
        return -1;
    }
    ierr = rtseis_ippsMedianFilter_initialize(
              ippsMedianIn.maskSize,  ippsMedianIn.lrt,
              ippsMedianIn.precision, ippsMedianOut);
    if (ierr != 0) 
    {
        RTSEIS_ERRMSG("%s", "ERror initializing ippsMedianOut");
        return -1;
    }
    if (ippsMedianIn.precision == RTSEIS_DOUBLE)
    {
        if (ippsMedianIn.maskSize > 1)
        {
            ippsCopy_64f(ippsMedianIn.dlysrc, ippsMedianOut->dlysrc,
                         ippsMedianIn.maskSize-1);    
            ippsCopy_64f(ippsMedianIn.dlydst, ippsMedianOut->dlydst,
                         ippsMedianIn.maskSize-1);
            ippsCopy_64f(ippsMedianIn.zi,     ippsMedianOut->zi,
                         ippsMedianIn.maskSize-1);
         }
    }
    else if (ippsMedianIn.precision == RTSEIS_FLOAT)
    {
        if (ippsMedianIn.maskSize > 1)
        {
            ippsCopy_32f(ippsMedianIn.dlysrc, ippsMedianOut->dlysrc,
                         ippsMedianIn.maskSize-1);    
            ippsCopy_32f(ippsMedianIn.dlydst, ippsMedianOut->dlydst,
                         ippsMedianIn.maskSize-1);
            ippsCopy_32f(ippsMedianIn.zi,     ippsMedianOut->zi,
                         ippsMedianIn.maskSize-1);
         }
    }
    else 
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the median filter for real-time or post-processing.
 * @param[in] lisRealTime  If true then the filter is for real-time
 *                         applications. \n
 *                         Otherwise, it is for post-processing.
 * @param[out] ippsMedian  On exit the real-time flag has been set.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_median
 */
int rtseis_ippsMedianFilter_toggleRealTime(
    const bool lisRealTime,
    struct ippsMedianFilter_struct *ippsMedian)
{
    ippsMedian->lrt = lisRealTime;
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the initial conditions for the median filter.  Note, that this
 *        will reset the filtering and therefore should be called prior to the
 *        first filter application.
 * @param[in] nz   Number of coefficients in zi.  This should equal the
 *                 maskSize - 1 where maskSize is the number of points
 *                 in the median.
 * @param[in] zi   The initial conditions.  This is an array of dimension [nz].
 * @param[in,out] ippsMedian  On input contains the initialized median filter.
 * @param[in,out] ippsMedian  On exit the initial conditions have been set to 
 *                            zi.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_median
 */
int rtseis_ippsMedianFilter_setInitialConditions(
    const int nz, const double zi[],
    struct ippsMedianFilter_struct *ippsMedian)
{
    double *dlysrc64, *zi64;
    float *dlysrc32, *zi32;
    int ierr;
    ierr = rtseis_ippsMedianFilter_resetInitialConditions(ippsMedian);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Error resetting filter");
        return -1; 
    }
    if (nz != ippsMedian->maskSize-1)
    {
        RTSEIS_ERRMSG("Error nz=%d not equal to maskSize-1=%d",
                    nz, ippsMedian->maskSize-1);
        return -1; 
    }
    if (nz == 0){return 0;}
    if (zi == NULL)
    {
        RTSEIS_ERRMSG("%s", "Error zi is NULL");
        return -1; 
    }
    if (ippsMedian->precision == RTSEIS_DOUBLE)
    {
        zi64 = (double *) ippsMedian->zi;
        dlysrc64 = (double *) ippsMedian->dlysrc;
        ippsCopy_64f(zi, zi64, nz-1);
        ippsCopy_64f(zi, dlysrc64, nz-1);
    }
    else if (ippsMedian->precision == RTSEIS_FLOAT)
    {
        zi32 = (float *) ippsMedian->zi;
        dlysrc32 = (float *) ippsMedian->dlysrc;
        ippsConvert_64f32f(zi, zi32, nz-1);
        ippsConvert_64f32f(zi, dlysrc32, nz-1);
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Resets the initial conditoins on the IPP median filter structure.
 *        The initial conditions for this filter were initially set with
 *        setInitialConditions.  If new initial conditions are desired then
 *        one should call that function.
 * @param[in,out] ippsMedian  On input contains the initial ippsMedian
 *                            structure. 
 * @param[in,out] ippsMedian  On exit the the initial conditions have been
 *                            reset.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_median
 */
int rtseis_ippsMedianFilter_resetInitialConditions(
    struct ippsMedianFilter_struct *ippsMedian)
{
    double *dlysrc64;
    float *dlysrc32;
    Ipp64f *zi64;
    Ipp32f *zi32;
    if (!ippsMedian->linit)
    {
        RTSEIS_ERRMSG("%s", "Filter never initialized");
        return -1;
    }
    if (ippsMedian->maskSize == 1){return 0;} // Nothing to do
    if (ippsMedian->precision == RTSEIS_DOUBLE)
    {
        dlysrc64 = (double *) ippsMedian->dlysrc;
        zi64 = (Ipp64f *) ippsMedian->zi;
        ippsCopy_64f(zi64, dlysrc64, ippsMedian->maskSize-1);
    }
    else if (ippsMedian->precision == RTSEIS_FLOAT)
    {
        dlysrc32 = (float *) ippsMedian->dlysrc;
        zi32 = (Ipp32f *) ippsMedian->zi;
        ippsCopy_32f(zi32, dlysrc32, ippsMedian->maskSize-1);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Applies the median filter to the array x.
 * @param[in] n      Number of points in x.
 * @param[in] x      Array to filter.  This has dimension [n].
 * @param[out] y     Median filtered version of y.  This has dimension [n].
 * @param[in,out] ippsMedian  On input contains the initialized median filter.
 * @param[in,out] ippsMedian  On exit contains the updated initial conditions.
 * @ingroup rtseis_ipps_median
 */
int rtseis_ippsMedianFilter_apply32f(
    const int n,
    const float x[],
    float y[],
    struct ippsMedianFilter_struct *ippsMedian)
{
    Ipp8u *pBuf;
    IppStatus status;
    float *dlysrc, *dlydst;
    double *x64, *y64;
    int ierr;
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    if (!ippsMedian->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsMedian not initialized");
        return -1;
    }
    if (ippsMedian->precision == RTSEIS_DOUBLE)
    {   
        x64 = ippsMalloc_64f(n);
        y64 = ippsMalloc_64f(n);
        ippsConvert_32f64f(x, x64, n); 
        ierr = rtseis_ippsMedianFilter_apply64f(n, x64, y64, ippsMedian);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Error calling ippsMedianFilter_apply64f");
            ierr = 1;
        }
        ippsConvert_64f32f(y64, y, n); 
        ippsFree(x64);
        ippsFree(y64);
        return ierr;
    }   
    pBuf = (Ipp8u *) ippsMedian->pBuf;
    dlysrc = (float *) ippsMedian->dlysrc;
    dlydst = NULL;
    if (ippsMedian->lrt)
    {
        dlydst = (float *) ippsMedian->dlydst;
    }
    status = ippsFilterMedian_32f(x, y, n, ippsMedian->maskSize,
                                  dlydst, dlysrc, pBuf);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Error applying filter");
        return -1; 
    }
    if (ippsMedian->lrt)
    {
        ippsCopy_32f(dlydst, dlysrc, ippsMedian->maskSize);
    }
    // Dererence pointers
    pBuf = NULL;
    dlysrc = NULL;
    dlydst = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Applies the median filter to the array x.
 * @param[in] n     Number of points in x.
 * @param[in] x     Array to filter.  This has dimension [n].
 * @param[out] y    Median filtered version of y.  This has dimension [n].
 * @param[in,out] ippsMedian  On input contains the initialized median filter.
 * @param[in,out] ippsMedian  On exit contains the updated initial conditions.
 * @ingroup rtseis_ipps_median
 */
int rtseis_ippsMedianFilter_apply64f(
    const int n,
    const double x[],
    double y[],
    struct ippsMedianFilter_struct *ippsMedian)
{
    Ipp8u *pBuf;
    IppStatus status;
    double *dlysrc, *dlydst;
    float *x32, *y32;
    int ierr;
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    if (!ippsMedian->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsMedian not initialized");
        return -1;
    }
    if (ippsMedian->precision == RTSEIS_FLOAT)
    {
        x32 = ippsMalloc_32f(n);
        y32 = ippsMalloc_32f(n);
        ippsConvert_64f32f(x, x32, n);
        ierr = rtseis_ippsMedianFilter_apply32f(n, x32, y32, ippsMedian);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Error calling ippsMedianFilter_apply32f");
            ierr = 1;
        }
        ippsConvert_32f64f(y32, y, n);
        ippsFree(x32);
        ippsFree(y32);
        return ierr;
    }
    pBuf = (Ipp8u *) ippsMedian->pBuf;
    dlysrc = (double *) ippsMedian->dlysrc;
    dlydst = NULL;
    if (ippsMedian->lrt)
    {
        dlydst = (double *) ippsMedian->dlydst;
    }
    status = ippsFilterMedian_64f(x, y, n, ippsMedian->maskSize,
                                  dlysrc, dlydst, pBuf);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Error applying filter");
        return -1;
    }
    if (ippsMedian->lrt && ippsMedian->maskSize > 2)
    {
        ippsCopy_64f(dlydst, dlysrc, ippsMedian->maskSize-1);
    }
    // Dererence pointers
    pBuf = NULL;
    dlysrc = NULL;
    dlydst = NULL;
    return 0;
}
                                    
