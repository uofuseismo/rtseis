#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/ippsHelper.h"
#include <ipps.h>

/*!
 * @defgroup rtseis_ipps_downsample Downsample
 * @brief Interface to the IPP downsampler.
 * @ingroup rtseis_ipps
 * @author Ben Baker
 * @copyright ISTI distributed under the Apache 2 license.
 */
/*!
 * @brief Initializes the downsampler. 
 * @param[in] downFactor   The downsampling factor.  This must be positive.
 * @param[in] lisRealTime  If true then the downsampler is for real-time
 *                         applications.
 * @param[in] lisRealTime  Otherwise, it is for post-processing.
 * @param[in] precision    Precision of the downsampler.
 *
 * @param[out] ippsDS      On successful exit this is the downsampling
 *                         structure.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_downsample
 */
int rtseis_ippsDownsample_initialize(const int downFactor,
                                     const bool lisRealTime,
                                     const enum rtseisPrecision_enum precision,
                                     struct ippsDownsample_struct *ippsDS)
{
    memset(ippsDS, 0, sizeof(struct ippsDownsample_struct));
    if (downFactor < 1)
    {
        RTSEIS_ERRMSG("Downsampling factor=%d must be positive", downFactor);
        return -1;
    }
    if (precision != RTSEIS_DOUBLE && precision != RTSEIS_FLOAT)
    {
        RTSEIS_ERRMSG("%s", "Only float64 and float32 supported");
        return -1;
    }
    ippsDS->precision = precision;
    ippsDS->downFactor = downFactor;
    rtseis_ippsDownsample_toggleRealTime(lisRealTime, ippsDS); 
    ippsDS->linit = true;
    return 0;
}
//============================================================================//
/*!
 * @brief Resets the downsampler phase to ippsDR->phase0.  To set phase0
 *        one must call setInitialConditions.
 * @param[in,out] ippsDS  On input contains the initialized downsampling
 *                        structure. \n
 *                        On exit the phase has been reset to 0.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_downsample
 */
int rtseis_ippsDownsample_resetInitialConditions(
    struct ippsDownsample_struct *ippsDS)
{
    if (!ippsDS->linit)
    { 
        RTSEIS_ERRMSG("%s", "downsampler not yet initialized");
        return -1;
    }
    ippsDS->phase = ippsDS->phase0;
    return 0;
}
//============================================================================//
/*!
 * @brief The downsampler does not have initial conditions.  This simply
 *        sets the phase.
 *  @param[in] phase       Phase of downsampler.  This must be in the range
 *                         [0, ippsDS->downFactor - 1].
 * @param[in,out] ippsDS   On in put this is the initialized downsampling
 *                         structure. \n
 *                         On exit the phase has been set.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_downsample
 */
int rtseis_ippsDownsample_setInitialConditions(
    const int phase, struct ippsDownsample_struct *ippsDS)
{
    int ierr;
    ierr = rtseis_ippsDownsample_resetInitialConditions(ippsDS);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Error resetting ippsDS");
        return -1;
    }
    if (phase < 0 || phase > ippsDS->downFactor - 1)
    {
        RTSEIS_ERRMSG("phase=%d must be in range[0,%d]",
                      phase, ippsDS->downFactor-1);
        return -1;
    }
    ippsDS->phase0 = phase;
    ippsDS->phase = phase;
    return 0; 
}
//============================================================================//
/*!
 * @brief Sets the downsampler for real-time or post-processing.
 * @param[in] lisRealTime  If true then the filter is for real-time
 *                         applications. \n
 *                         Otherwise, it is for post-processing.
 * @param[out] ippsDS      On exit the real-time flag has been set.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_downsample
 */
int rtseis_ippsDownsample_toggleRealTime(const bool lisRealTime,
                                         struct ippsDownsample_struct *ippsDS)
{
    ippsDS->lrt = lisRealTime;
    return 0;
}
//============================================================================//
/*!
 * @brief Applies the downsampler to the data.
 * @param[in] n           Number of points in input signal.
 * @param[in] x           Array to downsample.  This has dimension [n].
 * @param[in] ny          Max space allocated to y.  This must be at least
 *                        (n + downFactor - 1 - phase)/downFactor.
 * @param[out] len        The number of elements accessed in y.
 * @param[out] y          The downsampled version of x.  This is an array
 *                        of dimension [ny] however only the first len 
 *                        elements have been accessed.
 * @param[in,out] ippsDS  On input contains the initialized downsampling
 *                        structure.
 * @param[in,out] ippsDS  When in real-time mode the phase has been updated.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_downsample
 */
int rtseis_ippsDownsample_apply64f(const int n,
                                   const double x[],
                                   const int ny, int *len,
                                   double y[],
                                   struct ippsDownsample_struct *ippsDS)
{
    IppStatus status;
    int ierr, pDstLen, phase;
    if (x == NULL)
    {
        RTSEIS_ERRMSG("%s", "Error x is NULL");
        return -1;
    }
    if (!ippsDS->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsDS structure not intitialized");
        return -1; 
    }
    pDstLen = rtseis_ippsDownsample_estimateSpace(n, *ippsDS);
    if (ny < pDstLen || pDstLen < 0)
    {
        if (pDstLen < 0)
        {
            RTSEIS_ERRMSG("%s", "Space estimate error");
            return -1;
        }
        RTSEIS_ERRMSG("ny=%d must be at least length=%d", ny, pDstLen);
        return -1;
    }
    if (y == NULL)
    {
        RTSEIS_ERRMSG("%s", "Error y is NULL");
        return -1;
    }
    // Handle the wrong precision
    if (ippsDS->precision != RTSEIS_DOUBLE)
    {
        Ipp32f *x32 = ippsMalloc_32f(n);
        Ipp32f *y32 = ippsMalloc_32f(n);
        ippsConvert_64f32f(x, x32, n);
        ierr = rtseis_ippsDownsample_apply32f(n, x32, n, len, y32, ippsDS);
        if (ierr == 0)
        {
            ippsConvert_32f64f(y32, y, *len);
        }
        else
        {
            RTSEIS_ERRMSG("%s", "Error in apply32f");
        }
        ippsFree(x32);
        ippsFree(y32);
        return ierr;
    }
    // Set the phase and apply the downsampler
    phase = 0;
    if (ippsDS->lrt){phase = ippsDS->phase;}
    status = ippsSampleDown_64f(x, n, y, len,
                                ippsDS->downFactor, &phase);
    if (status != ippStsNoErr)      
    {
        RTSEIS_ERRMSG("%s", "Failed to downsample signal");
        return -1;
    }
    if (ippsDS->lrt){ippsDS->phase = phase;}
    return 0;
}
//============================================================================//
/*!
 * @brief Applies the downsampler to the data.
 * @param[in] n           Number of points in input signal.
 * @param[in] x           Array to downsample.  This has dimension [n].
 * @param[in] ny          Max space allocated to y.  This must be at least
 *                        (n + downFactor - 1 - phase)/downFactor.
 * @param[out] len        The number of elements accessed in y.
 * @param[out] y          The downsampled version of x.  This is an array
 *                        of dimension [ny] however only the first len 
 *                        elements have been accessed.
 * @param[in,out] ippsDS  On input contains the initialized downsampling
 *                        structure.
 * @param[in,out] ippsDS  When in real-time mode the phase has been updated.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_downsample
 */
int rtseis_ippsDownsample_apply32f(const int n,
                                   const float x[],
                                   const int ny, int *len,
                                   float y[],
                                   struct ippsDownsample_struct *ippsDS)
{
    IppStatus status;
    int ierr, pDstLen, phase;
    if (x == NULL)
    {
        RTSEIS_ERRMSG("%s", "Error x is NULL");
        return -1;
    }
    if (!ippsDS->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsDS structure not intitialized");
        return -1;
    }
    pDstLen = rtseis_ippsDownsample_estimateSpace(n, *ippsDS);
    if (ny < pDstLen || pDstLen < 0)
    {   
        if (pDstLen < 0)
        {   
            RTSEIS_ERRMSG("%s", "Space estimate error");
            return -1;
        }   
        RTSEIS_ERRMSG("ny=%d must be at least length=%d", ny, pDstLen);
        return -1; 
    }
    if (y == NULL)
    {
        RTSEIS_ERRMSG("%s", "Error y is NULL");
        return -1;
    }
    // Handle the wrong precision
    if (ippsDS->precision != RTSEIS_FLOAT)
    {
        Ipp64f *x64 = ippsMalloc_64f(n);
        Ipp64f *y64 = ippsMalloc_64f(n);
        ippsConvert_32f64f(x, x64, n);
        ierr = rtseis_ippsDownsample_apply64f(n, x64, n, len, y64, ippsDS);
        if (ierr == 0)
        {
            ippsConvert_64f32f(y64, y, *len);
        }
        else
        {
            RTSEIS_ERRMSG("%s", "Error in apply32f");
        }
        ippsFree(x64);
        ippsFree(y64);
        return ierr;
    }
    phase = 0;
    if (ippsDS->lrt){phase = ippsDS->phase;}
    status = ippsSampleDown_32f(x, n, y, len,
                                ippsDS->downFactor, &phase);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to downsample signal");
        return -1;
    }
    if (ippsDS->lrt){ippsDS->phase = phase;}
    return 0;
}
//============================================================================//
/*!
 * @brief Finalizes the IPP downsampling structure.
 * @param[in,out] ippsDS   On input this is the initialized downsampling
 *                         structure. \n
 *                         On exit all variables have been reset to 0.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_downsample
 */
int rtseis_ippsDownsample_finalize(struct ippsDownsample_struct *ippsDS)
{
    memset(ippsDS, 0, sizeof(struct ippsDownsample_struct));
    return 0;
} 
//============================================================================//
/*!
 * @brief Estimates the length of the downsampled signal.
 * @param[in] n       Number of points in the signal to downsample.
 * @param[in] ippsDS  The initialized downsampler.
 * @retval The estimated length of the downsampled array.
 * @retval If negative than error has occurred. 
 * @ingroup rtseis_ipps_downsample
 */
int rtseis_ippsDownsample_estimateSpace(
    const int n,
    const struct ippsDownsample_struct ippsDS)
{
    int pDstLen, phase;
    if (!ippsDS.linit)
    {
        RTSEIS_ERRMSG("%s", "ippsDS structure not intitialized");
        return 0;
    }
    phase = 0;
    if (ippsDS.lrt){phase = ippsDS.phase;} 
    pDstLen = (n + ippsDS.downFactor - 1 - phase)/ippsDS.downFactor;
    return pDstLen;
}
