#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/ippsHelper.h"
#include <ipps.h>

/*!
 * @defgroup rtseis_ipps_fir FIR Filtering
 * @brief Interface to the IPP finite impulse response filter.
 * @ingroup rtseis_ipps
 * @author Ben Baker
 * @copyright ISTI distributed under the Apache 2 license.
 */
/*!
 * @brief Initializes the FIR filter for IPP.
 *
 * @param[in] nb           Number of numerator coefficients.
 * @param[in] b            Numerator coefficients.  This is an array of
 *                         dimension [nb].
 * @param[in] lisRealTime  If true then the filter is for real-time
 *                         processing.
 * @param[in] lisRealTime  Otherwise, the filter is for post-processing.
 * @param[in] precision    Precision of FIR filter. 
 * 
 * @param[out] ippsFIR     FIR filter structure appropriate for IPP.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_fir
 */
int rtseis_ippsFIRFilter_initialize(const int nb, const double b[],
                                    const bool lisRealTime,
                                    const enum rtseisPrecision_enum precision,
                                    struct ippsFIRFilter_struct *ippsFIR)
{
    IppsFIRSpec_64f *pFIRspec64 = NULL;
    IppsFIRSpec_32f *pFIRspec32 = NULL;
    Ipp64f *pTaps64 = NULL;
    Ipp32f *pTaps32 = NULL;
    Ipp8u *pBuf = NULL;
    int tapsLen;
    IppStatus status;
    //------------------------------------------------------------------------//
    //
    // Initialize and error check
    memset(ippsFIR, 0, sizeof(struct ippsFIRFilter_struct));
    // Checks
    if (nb < 1 || b == NULL)
    {
        if (nb < 1){RTSEIS_ERRMSG("%s", "No b coefficients\n");}
        if (b == NULL){RTSEIS_ERRMSG("%s", "b is NULL\n");}
        return -1;
    }
    tapsLen = nb;
    ippsFIR->order = nb - 1;
    ippsFIR->nbDly = MAX(128, ippsFIR->order+1);
    ippsFIR->tapsRef = ippsMalloc_64f(nb);
    ippsCopy_64f(b, ippsFIR->tapsRef, nb);
    ippsFIR->zi = ippsMalloc_64f(MAX(1, ippsFIR->order));
    ippsZero_64f(ippsFIR->zi, MAX(1, ippsFIR->order));
    if (precision == RTSEIS_DOUBLE)
    {
        ippsFIR->dlysrc = ippsMalloc_64f(ippsFIR->nbDly);
        ippsFIR->dlydst = ippsMalloc_64f(ippsFIR->nbDly);
        ippsZero_64f(ippsFIR->dlysrc, ippsFIR->nbDly);
        ippsZero_64f(ippsFIR->dlydst, ippsFIR->nbDly);
        pTaps64 = ippsMalloc_64f(tapsLen);
        ippsZero_64f(pTaps64, tapsLen);
        ippsCopy_64f(b, pTaps64, nb);
        status = ippsFIRSRGetSize(tapsLen, ipp64f,
                                  &ippsFIR->specSize, &ippsFIR->bufferSize);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting state size");
            return -1;
        }
        pFIRspec64 = (IppsFIRSpec_64f *) ippsMalloc_8u(ippsFIR->specSize);
        pBuf = ippsMalloc_8u(ippsFIR->bufferSize);
        status = ippsFIRSRInit_64f(pTaps64, tapsLen, ippAlgDirect,
                                   pFIRspec64);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error initializing state structure");
            return -1;
        }
        ippsFIR->pTaps = (void *) pTaps64;
        ippsFIR->spec = (void *) pFIRspec64;
        ippsFIR->pBuf = (void *) pBuf;
    }
    else if (precision == RTSEIS_FLOAT) 
    {
        ippsFIR->dlysrc = ippsMalloc_32f(ippsFIR->nbDly);
        ippsFIR->dlydst = ippsMalloc_32f(ippsFIR->nbDly);
        ippsZero_32f(ippsFIR->dlysrc, ippsFIR->nbDly);
        ippsZero_32f(ippsFIR->dlydst, ippsFIR->nbDly);
        pTaps32 = ippsMalloc_32f(tapsLen);
        ippsZero_32f(pTaps32, tapsLen);
        ippsConvert_64f32f(b, pTaps32, nb);
        status = ippsFIRSRGetSize(tapsLen, ipp32f,
                                  &ippsFIR->specSize, &ippsFIR->bufferSize);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting state size");
            return -1; 
        }
        pFIRspec32 = (IppsFIRSpec_32f *) ippsMalloc_8u(ippsFIR->specSize);
        pBuf = ippsMalloc_8u(ippsFIR->bufferSize);
        status = ippsFIRSRInit_32f(pTaps32, tapsLen, ippAlgDirect,
                                   pFIRspec32);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error initializing state structure");
            return -1; 
        }
        ippsFIR->pTaps = (void *) pTaps32;
        ippsFIR->spec = (void *) pFIRspec32;
        ippsFIR->pBuf = (void *) pBuf;
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    rtseis_ippsFIRFilter_toggleRealTime(lisRealTime, ippsFIR);
    ippsFIR->precision = precision;
    ippsFIR->linit = true;
    return 0;
}
//============================================================================//
/*!
 * @brief Copies the initialized FIR structure ippsFIRIn to ippsFIROut.
 *        Note, that this will also copy the initial conditions and delay
 *        lines.
 *
 * @param[in] ippsFIRIn    The ippsFIR structure to copy.
 * @param[out] ippsFIROut  The copied version of ippsFIRIn. 
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_fir
 */
int rtseis_ippsFIRFilter_copy(const struct ippsFIRFilter_struct ippsFIRIn,
                              struct ippsFIRFilter_struct *ippsFIROut)
{
    int ierr, n;
    memset(ippsFIROut, 0, sizeof(struct ippsFIRFilter_struct));
    if (!ippsFIRIn.linit)
    {
        RTSEIS_ERRMSG("%s", "Input FIR filter not yet set");
        return -1;
    }
    n = ippsFIRIn.order + 1; // Number of taps
    ierr =  rtseis_ippsFIRFilter_initialize(n, ippsFIRIn.tapsRef,
                                            ippsFIRIn.lrt, ippsFIRIn.precision,
                                            ippsFIROut);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Error initializing ippsFIROut");
        return -1;
    }
    if (ippsFIRIn.precision == RTSEIS_DOUBLE)
    {
        ippsCopy_64f(ippsFIRIn.dlysrc, ippsFIROut->dlysrc, ippsFIRIn.nbDly);
        ippsCopy_64f(ippsFIRIn.dlydst, ippsFIROut->dlydst, ippsFIRIn.nbDly);
        if (ippsFIRIn.order > 0)
        {
            ippsCopy_64f(ippsFIRIn.zi, ippsFIROut->zi, ippsFIRIn.order); 
        }
    }
    else if (ippsFIRIn.precision == RTSEIS_FLOAT)
    {
        ippsCopy_32f(ippsFIRIn.dlysrc, ippsFIROut->dlysrc, ippsFIRIn.nbDly);
        ippsCopy_32f(ippsFIRIn.dlydst, ippsFIROut->dlydst, ippsFIRIn.nbDly);
        if (ippsFIRIn.order > 0)
        {
            ippsCopy_64f(ippsFIRIn.zi, ippsFIROut->zi, ippsFIRIn.order); 
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
 * @brief Resets the initial conditions on the IPP FIR filter structure.
 *        The initial conditions for this filter were initially set with
 *        setInitialConditions.  If new initial conditions are desired
 *        then one should call that function.
 *
 * @param[in,out] ippsFIR  On input contains the initialized ippsFIR filter
 *                         structure. \n
 *                         On exit dlysrc has been reset to the initial
 *                         conditions.
 *
 * @result 0 indicates success.
 *
 * @ingroup rtseis_ipps_fir
 *
 */
int rtseis_ippsFIRFilter_resetInitialConditions(
    struct ippsFIRFilter_struct *ippsFIR)
{
    double *dlysrc64;
    float *dlysrc32;
    Ipp64f *zi64;
    if (!ippsFIR->linit)
    {
        RTSEIS_ERRMSG("%s", "Filter never initialized");
        return -1;
    }
    zi64 = (Ipp64f *) ippsFIR->zi;
    if (ippsFIR->precision == RTSEIS_DOUBLE)
    {
        dlysrc64 = (double *) ippsFIR->dlysrc;
        ippsZero_64f(dlysrc64, ippsFIR->nbDly);
        if (ippsFIR->order > 0)
        {
            ippsCopy_64f(zi64, dlysrc64, ippsFIR->order);
        }
    }
    else if (ippsFIR->precision == RTSEIS_FLOAT)
    {
        dlysrc32 = (float *) ippsFIR->dlysrc;
        ippsZero_32f(dlysrc32, ippsFIR->nbDly);
        if (ippsFIR->order > 0)
        {
            ippsConvert_64f32f(zi64, dlysrc32, ippsFIR->order);
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the initial conditions on the FIR filter.  Note that this
 *        will reset the filtering and should therefore be called prior
 *        to the first filter application.
 *
 * @param[in] nz     Number of coefficients in zi.  This must be equal
 *                   to ippsFIR->order.
 * @param[in] zi     Initial conditions.  This is an array of dimension [nz].
 * @param[in,out] ippsFIR  On input contains the initialized FIR filter 
 *                         structure. \n
 *                         On exit the initial conditions have been set.
 *
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_fir
 */
int rtseis_ippsFIRFilter_setInitialConditions(
    const int nz, const double zi[],
    struct ippsFIRFilter_struct *ippsFIR)
{
    double *dlysrc64;
    float *dlysrc32;
    Ipp64f *zi64;
    int ierr;
    ierr = rtseis_ippsFIRFilter_resetInitialConditions(ippsFIR);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Error resetting filter");
        return -1;
    }
    if (nz != ippsFIR->order)
    {
        RTSEIS_ERRMSG("Error nz=%d not equal to order=%d",
                    nz, ippsFIR->order);
        return -1;
    }   
    if (nz == 0){return 0;}
    if (zi == NULL)
    {
        RTSEIS_ERRMSG("%s", "Error zi is NULL");
        return -1; 
    }
    zi64 = (Ipp64f *) ippsFIR->zi;
    ippsCopy_64f(zi, zi64, nz);
    if (ippsFIR->precision == RTSEIS_DOUBLE)
    {
        dlysrc64 = (double *) ippsFIR->dlysrc;
        ippsCopy_64f(zi, dlysrc64, nz);
    }   
    else if (ippsFIR->precision == RTSEIS_FLOAT)
    {
        dlysrc32 = (float *) ippsFIR->dlysrc;
        ippsConvert_64f32f(zi, dlysrc32, nz);
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
 * @brief Applies the FIR filtering to the double array x.
 *
 * @param[in] n            Number of points in x.
 * @param[in] x            Array to filter.  This has dimension [n].
 *
 * @param[out] y           Filtered version of x.  This has dimension [n].
 *
 * @param[in,out] ippsFIR  Holds the FIR filtering states.
 *
 * @result 0 indicates success.
 *
 * @ingroup rtseis_ipps_fir
 *
 */
int rtseis_ippsFIRFilter_apply64f(const int n,
                                  const double x[],
                                  double y[],
                                  struct ippsFIRFilter_struct *ippsFIR)
{
    IppsFIRSpec_64f *pFIRspec = NULL;
    Ipp8u *pBuf;
    double *dlysrc, *dlydst;
    float *x32, *y32;
    int ierr;
    IppStatus status;
    ierr = 0;
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {  
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    } 
    if (!ippsFIR->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsFIR not initialized");
        return -1;
    }
    if (ippsFIR->precision == RTSEIS_FLOAT)
    {
        x32 = ippsMalloc_32f(n);
        y32 = ippsMalloc_32f(n);
        ippsConvert_64f32f(x, x32, n);
        ierr = rtseis_ippsFIRFilter_apply32f(n, x32, y32, ippsFIR);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Error calling rtseis_ippsFIRFilter_apply32f");
            ierr = 1;
        } 
        ippsConvert_32f64f(y32, y, n);
        ippsFree(x32);
        ippsFree(y32);
        return ierr;
    }
    // Can't deal with any other precision
    if (ippsFIR->precision != RTSEIS_DOUBLE)
    {
        RTSEIS_ERRMSG("%s", "Error must convert data");
        return -1;
    }
    // Apply the filter
    pFIRspec = (IppsFIRSpec_64f *) ippsFIR->spec;
    pBuf = (Ipp8u *) ippsFIR->pBuf;
    dlysrc = (double *) ippsFIR->dlysrc; 
    dlydst = NULL; 
    if (ippsFIR->order > 0) //ippsFIR->lrt && ippsFIR->order > 0)
    {
        dlydst = (double *) ippsFIR->dlydst;
    }
    status = ippsFIRSR_64f(x, y, n, pFIRspec, dlysrc, dlydst, pBuf);
    if (status != ippStsNoErr)
    { 
        RTSEIS_ERRMSG("%s", "Error applying filter");
        return -1;
    }
    // Update the initial conditions
    if (ippsFIR->lrt && ippsFIR->order > 0)
    {
        ippsCopy_64f(dlydst, dlysrc, ippsFIR->order);
    }
    // Dereference pointers
    //pFIRspec = NULL;
    //pBuf = NULL;
    return 0;
} 
//============================================================================//
/*!
 * @brief Applies the FIR filtering to the float array x.
 * @param[in] n            Number of points in x.
 * @param[in] x            Array to filter.  This has dimension [n].
 * @param[out] y           Filtered version of x.  This has dimension [n].
 * @param[in,out] ippsFIR  Holds the FIR filtering states.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_fir
 */
int rtseis_ippsFIRFilter_apply32f(const int n,
                                  const float x[],
                                  float y[],
                                  struct ippsFIRFilter_struct *ippsFIR)
{
    IppsFIRSpec_32f *pFIRspec = NULL;
    Ipp8u *pBuf;
    float *dlysrc, *dlydst;
    double *x64, *y64;
    int ierr;
    IppStatus status;
    ierr = 0;
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {   
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1; 
    }   
    if (!ippsFIR->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsFIR not initialized");
        return -1; 
    }   
    if (ippsFIR->precision == RTSEIS_FLOAT)
    {
        x64 = ippsMalloc_64f(n);
        y64 = ippsMalloc_64f(n);
        ippsConvert_32f64f(x, x64, n);
        ierr = rtseis_ippsFIRFilter_apply64f(n, x64, y64, ippsFIR);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Error calling rtseis_ippsFIRFilter_apply32f");
            ierr = 1;
        }
        ippsConvert_64f32f(y64, y, n);
        ippsFree(x64);
        ippsFree(y64);
        return ierr;
    }
    // Can't deal with any other precision
    if (ippsFIR->precision != RTSEIS_FLOAT)
    {
        RTSEIS_ERRMSG("%s", "Error must convert data");
        return -1;
    }
    // Apply the filter
    pFIRspec = (IppsFIRSpec_32f *) ippsFIR->spec;
    pBuf = (Ipp8u *) ippsFIR->pBuf;
    dlysrc = (float *) ippsFIR->dlysrc;
    dlydst = NULL; 
    if (ippsFIR->order > 0) //ippsFIR->lrt && ippsFIR->order > 0)
    {
        dlydst = (float *) ippsFIR->dlydst;
    }
    status = ippsFIRSR_32f(x, y, n, pFIRspec, dlysrc, dlydst, pBuf);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Error applying filter");
        return -1;
    }
    // Update the initial conditions
    if (ippsFIR->lrt && ippsFIR->order > 0)
    {
        ippsCopy_32f(dlydst, dlysrc, ippsFIR->order);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the FIR filter for real-time or post-processing.
 * @param[in] lisRealTime  If true then the filter is for real-time
 *                         applications. \n
 *                         Otherwise, it is for post-processing.
 *
 * @param[out] ippsFIR     On exit the real-time flag has been set.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_fir
 */
int rtseis_ippsFIRFilter_toggleRealTime(const bool lisRealTime,
                                        struct ippsFIRFilter_struct *ippsFIR)
{
    ippsFIR->lrt = lisRealTime;
    return 0;
}
//============================================================================//
/*!
 * @brief Determines if the filter is real-time or for post-processing.
 * @param[in] ippsFIR   Initialized filter structure.
 * @retval True if the filter is for real-time.  \n
 * @retval False if the filter is for post-processing. 
 * @ingroup rtseis_ipps_fir
 */
bool rtseis_ippsFIRFilter_isRealTime(const struct ippsFIRFilter_struct ippsFIR)
{
    if (!ippsFIR.linit)
    {
        RTSEIS_WARNMSG("%s", "ippsFIR not initialized");
    }
    return ippsFIR.lrt;
} 
//============================================================================//
/*!
 * @brief Frees all memory on the FIR filter structure and sets all variables
 *        to 0.
 * @param[in,out] ippsFIR  On input this is an initialized FIR filter 
 *                         structure. \n
 *                         On exit all memory has been freed and variables set
 *                         to 0.
 * @result 0 indicates success.
 * @ingroup rtseis_ipp_fir
 */
int rtseis_ippsFIRFilter_finalize(struct ippsFIRFilter_struct *ippsFIR)
{
    Ipp64f *pTaps64;
    Ipp32f *pTaps32;
    IppsFIRSpec_64f *pFIRspec64 = NULL;
    IppsFIRSpec_32f *pFIRspec32 = NULL;
    Ipp8u *pBuf;
    if (!ippsFIR->linit){return 0;} // Never initialized 
    pBuf = (Ipp8u *) ippsFIR->pBuf;
    ippsFree(pBuf);
    ippsFree(ippsFIR->zi);
    ippsFree(ippsFIR->dlysrc);
    ippsFree(ippsFIR->dlydst);
    if (ippsFIR->precision == RTSEIS_DOUBLE)
    {
        pFIRspec64 = (IppsFIRSpec_64f *) ippsFIR->spec;
        pTaps64 = (Ipp64f *) ippsFIR->pTaps;
        ippsFree(pTaps64);
        ippsFree(pFIRspec64);
    }
    else if (ippsFIR->precision == RTSEIS_FLOAT)
    {
        pFIRspec32 = (IppsFIRSpec_32f *) ippsFIR->spec;
        pTaps32 = (Ipp32f *) ippsFIR->pTaps;
        ippsFree(pTaps32);
        ippsFree(pFIRspec32);
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    ippsFree(ippsFIR->tapsRef);
    memset(ippsFIR, 0, sizeof(struct ippsFIRFilter_struct));
    return 0;
}
