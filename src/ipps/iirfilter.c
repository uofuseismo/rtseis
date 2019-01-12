#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/ippsHelper.h"
#include <ippversion.h>
#include <ippcore.h>
#include <ipps.h>

#if IPP_VERSION_MAJOR < 2019
static int iirdf2_64f(const int n, const double x[],
                      double y[],
                      struct ippsIIRFilter_struct *iir);

static int iirdf2_32f(const int n, const float x[],
                      float y[],
                      struct ippsIIRFilter_struct *iir);
#endif
int rtseis_ippsIIRFilter_apply64f_ori(const int n, const double x[], 
                                      double y[], 
                                      struct ippsIIRFilter_struct *ippsIIR);
int rtseis_ippsIIRFilter_apply32f_ori(const int n, const float x[],
                                      float y[],
                                      struct ippsIIRFilter_struct *ippsIIR);

/*!
 * @defgroup rtseis_ipps_iir IIR Filtering
 * @brief Interface to the IPP infinite impulse response filter.
 * @ingroup rtseis_ipps
 * @author Ben Baker
 * @copyright ISTI distributed under the Apache 2 license.
 */
/*!
 * @brief Initializes the IIR filter for IPP.
 * @param[in] nb           Number of numerator coefficients.
 * @param[in] b            Numerator coefficients.  This is an array of
 *                         dimension [nb].
 * @param[in] na           Number of denominator coefficients.
 * @param[in] a            Denominator coefficients.  This is an array of
 *                         dimension [na].
 * @param[in] lisRealTime  If true then the filter is for real-time
 *                         applications.
 * @param[in] lisRealTime  Otherwise, it is for post-processing.
 * @param[in] precision    Precision of IIR filter. 
 * 
 * @param[out] ippsIIR     IIR filter structure appropriate for IPP.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFilter_initialize(const int nb, const double b[],
                                    const int na, const double a[],
                                    const bool lisRealTime,
                                    const enum rtseisPrecision_enum precision,
                                    struct ippsIIRFilter_struct *ippsIIR)
{
    Ipp64u featureMask, enabledMask;
    //Ipp32u pCpuidInfoRegs[4];
    IppsIIRState_64f *pIIRState64 = NULL;
    IppsIIRState_32f *pIIRState32 = NULL;
    Ipp64f *zi64 = NULL;
    Ipp64f *pTaps64 = NULL;
    Ipp32f *zi32 = NULL;
    Ipp32f *pTaps32 = NULL;
    Ipp8u *pBuf = NULL;
    IppStatus status;
    double a0;
    int bufSrcLen, order;
    bool useLegacy;
    //------------------------------------------------------------------------//
    //
    // Initialize and error check
    memset(ippsIIR, 0, sizeof(struct ippsIIRFilter_struct));
    // Checks
    if (nb < 1 || na < 1 || b == NULL || a == NULL)
    {
        if (nb < 1){RTSEIS_ERRMSG("%s", "No b coefficients\n");}
        if (na < 1){RTSEIS_ERRMSG("%s", "No a coefficients\n");}
        if (b == NULL){RTSEIS_ERRMSG("%s", "a is NULL\n");}
        if (a == NULL){RTSEIS_ERRMSG("%s", "b is NULL\n");}
        return -1;
    }
    a0 = a[0];
    if (a0 == 0.0)
    {
        RTSEIS_ERRMSG("%s", "Error a[0] cannot be zero\n");
        return -1;
    }
    // Get processor info
    order = MAX(nb, na) - 1;
    useLegacy = false;
#if IPP_VERSION_MAJOR > 2018
    if (order > 8)
    {
        status = ippGetCpuFeatures(&featureMask, NULL);
        if (status == ippStsNoErr)
        {
            enabledMask = ippGetEnabledCpuFeatures();
            if ((featureMask & ippCPUID_AVX2) &&
                (enabledMask & ippCPUID_AVX2))
            {
                useLegacy = true;
            }
            if ((featureMask & ippCPUID_AVX512F) &&
                (enabledMask & ippCPUID_AVX512F))
            {
                useLegacy = true;
            }
        }
        else
        {
            RTSEIS_ERRMSG("%s", "Failed to get processor info");
            return -1;
        }
    }
#endif
    if (useLegacy){RTSEIS_WARNMSG("%s", "Use legacy code for IIR order > 8");}

    ippsIIR->useLegacy = useLegacy;
    ippsIIR->nbRef = nb;
    ippsIIR->naRef = na;
    ippsIIR->bRef = ippsMalloc_64f(MAX(nb, na)); //ippsIIR->nbRef);
    ippsIIR->aRef = ippsMalloc_64f(MAX(nb, na)); //ippsIIR->naRef);
    // Save the normalized coefficients
    ippsZero_64f(ippsIIR->bRef, MAX(nb, na));
    ippsZero_64f(ippsIIR->aRef, MAX(nb, na));
    ippsDivC_64f(b, a0, ippsIIR->bRef, nb);
    ippsDivC_64f(a, a0, ippsIIR->aRef, na);
    ippsIIR->order = order; //MAX(nb, na) - 1;
    ippsIIR->nbDly = MAX(128, ippsIIR->order+1);
    bufSrcLen = 2*MAX(ippsIIR->order+1, 1024);
    if (precision == RTSEIS_DOUBLE)
    {
        ippsIIR->bufsrc = ippsMalloc_64f(bufSrcLen);
        ippsIIR->bufdst = ippsMalloc_64f(bufSrcLen);
        ippsIIR->bufdly = ippsMalloc_64f(ippsIIR->nbDly);
        ippsZero_64f(ippsIIR->bufsrc, bufSrcLen);
        ippsZero_64f(ippsIIR->bufdst, bufSrcLen);
        ippsZero_64f(ippsIIR->bufdly, ippsIIR->nbDly);
        pTaps64 = ippsMalloc_64f(2*(ippsIIR->order+1));
        zi64 = ippsMalloc_64f(MAX(1, ippsIIR->order));
        ippsZero_64f(pTaps64, 2*(ippsIIR->order+1));
        ippsZero_64f(zi64, MAX(1, ippsIIR->order));
        ippsCopy_64f(b, &pTaps64[0],                nb);
        ippsCopy_64f(a, &pTaps64[ippsIIR->order+1], na);
        ippsIIR->pTaps = (void *) pTaps64;
        status = ippsIIRGetStateSize_64f(ippsIIR->order, &ippsIIR->bufferSize);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting state size\n");
            return -1;
        }
        pBuf = ippsMalloc_8u(ippsIIR->bufferSize);
        status = ippsIIRInit_64f(&pIIRState64, pTaps64, ippsIIR->order,
                                 zi64, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error initializing state structure\n");
            return -1;
        }
        ippsIIR->state = (void *) pIIRState64;
        ippsIIR->pBuf = (void *) pBuf;
        ippsIIR->zi = (void *) zi64;
    }
    else if (precision == RTSEIS_FLOAT)
    {
        ippsIIR->bufsrc = ippsMalloc_32f(bufSrcLen);
        ippsIIR->bufdst = ippsMalloc_32f(bufSrcLen);
        ippsIIR->bufdly = ippsMalloc_32f(ippsIIR->nbDly);
        ippsZero_32f(ippsIIR->bufsrc, bufSrcLen);
        ippsZero_32f(ippsIIR->bufdst, bufSrcLen);
        ippsZero_32f(ippsIIR->bufdly, ippsIIR->nbDly);
        pTaps32 = ippsMalloc_32f(2*(ippsIIR->order+1));
        zi32 = ippsMalloc_32f(MAX(ippsIIR->order,1));
        ippsZero_32f(pTaps32, 2*(ippsIIR->order+1));
        ippsZero_32f(zi32, MAX(1,ippsIIR->order));
        ippsConvert_64f32f(b, &pTaps32[0],                nb);
        ippsConvert_64f32f(a, &pTaps32[ippsIIR->order+1], na);
        ippsIIR->pTaps = (void *) pTaps32;
        status = ippsIIRGetStateSize_32f(ippsIIR->order, &ippsIIR->bufferSize);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting state size\n");
            return -1;
        }
        pBuf = ippsMalloc_8u(ippsIIR->bufferSize);
        status = ippsIIRInit_32f(&pIIRState32, pTaps32, ippsIIR->order,
                                 zi32, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error initializing state structure\n");
            return -1;
        }
        ippsIIR->state = (void *) pIIRState32;
        ippsIIR->pBuf = (void *) pBuf;
        ippsIIR->zi = (void *) zi32; 
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    rtseis_ippsIIRFilter_toggleRealTime(lisRealTime, ippsIIR);
    ippsIIR->precision = precision;
    ippsIIR->linit = true;
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the IIR filter for real-time or post-processing.
 * @param[in] lisRealTime  If true then the filter is for real-time
 *                         applications.
 * @param[in] lisRealTime  Otherwise, it is for post-processing.
 * @param[out] ippsIIR     On exit the real-time flag has been set.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFilter_toggleRealTime(const bool lisRealTime,
                                        struct ippsIIRFilter_struct *ippsIIR)
{
    ippsIIR->lrt = lisRealTime;
    return 0;
}
//============================================================================//
/*!
 * @brief Copies an ipps IIR filter structure.  Note, that this will also
 *        copy the initial conditions and delay lines.
 * @param[in] ippsIIRIn    The ippsIIR filter structure to copy.
 * @param[out] ippsIIROut  The copied version of ippsIIRIn. 
 * @result 0 indicate success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFIlter_copy(const struct ippsIIRFilter_struct ippsIIRIn,
                            struct ippsIIRFilter_struct *ippsIIROut)
{
    int bufSrcLen, ierr;
    memset(ippsIIROut, 0, sizeof(struct ippsIIRFilter_struct));
    if (!ippsIIRIn.linit)
    {
        RTSEIS_ERRMSG("%s", "Input IIR filter not yet set");
        return -1;
    }
    ierr = rtseis_ippsIIRFilter_initialize(ippsIIRIn.nbRef, ippsIIRIn.bRef,
                                           ippsIIRIn.naRef, ippsIIRIn.aRef,
                                           ippsIIRIn.lrt,
                                         ippsIIRIn.precision,
                                         ippsIIROut);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Error initializing ippsIIROut");
        return -1;
    }
    bufSrcLen = 2*MAX(ippsIIROut->order+1, 1024);
    ippsIIROut->lnextBlock = false;
    if (ippsIIRIn.precision == RTSEIS_DOUBLE)
    {
        ippsCopy_64f(ippsIIRIn.bufsrc, ippsIIROut->bufsrc, bufSrcLen);
        ippsCopy_64f(ippsIIRIn.bufdst, ippsIIROut->bufdst, bufSrcLen);
        ippsCopy_64f(ippsIIRIn.bufdly, ippsIIROut->bufdly, ippsIIRIn.nbDly);
        if (ippsIIRIn.order > 0)
        {
            ippsCopy_64f(ippsIIRIn.zi, ippsIIROut->zi, ippsIIRIn.order);
        }
    }
    else if (ippsIIRIn.precision == RTSEIS_FLOAT)
    {
        ippsCopy_32f(ippsIIRIn.bufsrc, ippsIIROut->bufsrc, bufSrcLen);
        ippsCopy_32f(ippsIIRIn.bufdst, ippsIIROut->bufdst, bufSrcLen);
        ippsCopy_32f(ippsIIRIn.bufdly, ippsIIROut->bufdly, ippsIIRIn.nbDly);
        if (ippsIIRIn.order > 0)
        {
            ippsCopy_32f(ippsIIRIn.zi, ippsIIROut->zi, ippsIIRIn.order);
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
 * @brief Resets the initial conditions on the IPP IIR filter structure.
 *        The initial conditions for this filter were initially set with
 *        setInitialConditions.  If new initial conditions are desired 
 *        then one should call that function.
 * @param[in,out] ippsIIR   On input contains the initialized ippsIIR filter
 *                          structure.
 * @param[in,out] ippsIIR   On exit the accumLen and workspace buffers have
 *                          been set to zero, and bufsrc reset to the initial
 *                          conditions.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFilter_resetInitialConditions(
    struct ippsIIRFilter_struct *ippsIIR)
{
    double *bufdly64, *bufsrc64, *bufdst64;
    float *bufdly32, *bufsrc32, *bufdst32;
    int bufSrcLen;
    Ipp64f *zi64;
    Ipp32f *zi32;
    if (!ippsIIR->linit)
    {
        RTSEIS_ERRMSG("%s", "Error ippsIIR not initialized");
        return -1;
    }
    ippsIIR->lnextBlock = false;
    ippsIIR->accumLen = 0;
    bufSrcLen = 2*MAX(ippsIIR->order+1, 1024);
    if (ippsIIR->precision == RTSEIS_DOUBLE)
    {
        bufsrc64 = (double *) ippsIIR->bufsrc;
        bufdst64 = (double *) ippsIIR->bufdst;
        bufdly64 = (double *) ippsIIR->bufdly;
        zi64 = (Ipp64f *) ippsIIR->zi;
        ippsZero_64f(bufsrc64, bufSrcLen);
        ippsZero_64f(bufdst64, bufSrcLen);
        ippsZero_64f(bufdly64, ippsIIR->nbDly); 
        if (ippsIIR->order > 0)
        {
            ippsCopy_64f(zi64, bufsrc64, ippsIIR->order);
        }
    }
    else if (ippsIIR->precision == RTSEIS_FLOAT) 
    {
        bufsrc32 = (float *) ippsIIR->bufsrc;
        bufdst32 = (float *) ippsIIR->bufdst;
        bufdly32 = (float *) ippsIIR->bufdly;
        zi32 = (Ipp32f *) ippsIIR->zi;
        ippsZero_32f(bufsrc32, bufSrcLen);
        ippsZero_32f(bufdst32, bufSrcLen);
        ippsZero_32f(bufdly32, ippsIIR->nbDly);
        if (ippsIIR->order > 0)
        { 
            ippsCopy_32f(zi32, bufsrc32, ippsIIR->order);
        }
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the initial conditions for the filter.  Note, that this will
 *        reset the filtering and should therefore be called prior to the
 *        first filter application.
 *
 * @param[in] nz           Number of coefficients in zi.  This must be equal
 *                         to the filter order which is ippsIIR->order; 
 * @param[in] zi           Initial conditions.  This is an array of dimension
 *                         [nz].
 * @param[in,out] ippsIIR  On input contains the initialized IIR filter
 *                         structure. \n
 *                         On output the filter buffers have been reset and
 *                         the initial conditions set.
 * 
 * @result 0 indicates success.
 *
 * @ingroup rtseis_ipps_iir
 *
 */
int rtseis_ippsIIRFilter_setInitialConditions(
    const int nz, const double zi[],
    struct ippsIIRFilter_struct *ippsIIR)
{
    double *bufdly64;
    float *bufdly32;
    Ipp64f *zi64;
    Ipp32f *zi32;
    int ierr;
    ierr = rtseis_ippsIIRFilter_resetInitialConditions(ippsIIR);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Error resetting initial conditions");
        return -1;
    }
    if (nz != ippsIIR->order)
    {
        RTSEIS_ERRMSG("Error nz=%d not equal to order=%d",
                      nz, ippsIIR->order);
        return -1;
    }
    if (nz == 0){return 0;}
    if (zi == NULL)
    {
        RTSEIS_ERRMSG("%s", "Error zi is NULL");
        return -1;
    }
    if (ippsIIR->precision == RTSEIS_DOUBLE)
    {
        zi64 = (Ipp64f *) ippsIIR->zi;
        bufdly64 = (double *) ippsIIR->bufdly;
        ippsCopy_64f(zi, zi64, nz);
        ippsCopy_64f(zi, bufdly64, nz);
    }
    else if (ippsIIR->precision == RTSEIS_FLOAT)
    {
        zi32 = (Ipp32f *) ippsIIR->zi;
        bufdly32 = (float *) ippsIIR->bufdly;
        ippsConvert_64f32f(zi, zi32, nz);
        ippsConvert_64f32f(zi, bufdly32, nz);
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
 * @brief Resets the final conditions which aren't applicable.  This function
 *        exists for completeness.
 * @param[in] ippsIIR   This is an initialized IPP IIR filter structure. 
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFilter_resetFinalConditions(
    struct ippsIIRFilter_struct *ippsIIR)
{
    if (!ippsIIR->linit)
    {
        RTSEIS_ERRMSG("%s", "Error struct not initialized");
        return -1;
    }
    return 0;
}
/*!
 * @brief Applies the IIR filtering to the array x.
 * @param[in] n            Number of points in x.
 * @param[in] x            Array to filter.  This has dimension [n].
 * @param[out] y           Filtered version of x.  This has dimension [n].
 * @param[in,out] ippsIIR  Holds the IIR filtering states.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFilter_apply64f(const int n, const double x[],
                                  double y[],
                                  struct ippsIIRFilter_struct *ippsIIR)
{
    int ierr;
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    // Verify the structure was initialized.
    if (!ippsIIR->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsIIR not initialized");
        return -1;
    }
    if (ippsIIR->precision != RTSEIS_DOUBLE)
    {
        RTSEIS_ERRMSG("%s", "Error must convert data");
        return -1; 
    }
    ierr = 0;
    if (!ippsIIR->useLegacy)
    {
        IppsIIRState_64f *pIIRState = (IppsIIRState_64f *) ippsIIR->state;
        if (!ippsIIR->lnextBlock)
        {
            Ipp64f *bufdly = (Ipp64f *) ippsIIR->bufdly;
            IppStatus status = ippsIIRSetDlyLine_64f(pIIRState, bufdly);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error filtering");
                return -1;
            }
            if (ippsIIR->lrt){ippsIIR->lnextBlock = true;}
        }
        ippsIIR_64f(x, y, n, pIIRState);
        if (ippsIIR->lrt)
        {
            ippsIIRGetDlyLine_64f(pIIRState, ippsIIR->bufdly);
        }
    }
    else
    {
        //iirdf2_64f(n, x, y, ippsIIR);
        ierr = rtseis_ippsIIRFilter_apply64f_ori(n, x, y, ippsIIR);
    }
    return ierr;
}
/*!
 * @brief Applies the IIR filtering to the array x.
 *
 * @param[in] n            Number of points in x.
 * @param[in] x            Array to filter.  This has dimension [n].
 * @param[out] y           Filtered version of x.  This has dimension [n].
 * @param[in,out] ippsIIR  Holds the IIR filtering states.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFilter_apply32f(const int n, const float x[],
                                  float y[],
                                  struct ippsIIRFilter_struct *ippsIIR)
{
    int ierr;
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    // Verify the structure was initialized.
    if (!ippsIIR->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsIIR not initialized");
        return -1;
    }
    if (ippsIIR->precision != RTSEIS_FLOAT)
    {
        RTSEIS_ERRMSG("%s", "Error must convert data");
        return -1;
    }
    ierr = 0; 
    if (!ippsIIR->useLegacy)
    {
        IppsIIRState_32f *pIIRState = (IppsIIRState_32f *) ippsIIR->state;
        if (!ippsIIR->lnextBlock)
        {
            Ipp32f *bufdly = (Ipp32f *) ippsIIR->bufdly;
            IppStatus status = ippsIIRSetDlyLine_32f(pIIRState, bufdly);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error filtering");
                return -1;
            }
            if (ippsIIR->lrt){ippsIIR->lnextBlock = true;}
        }
        ippsIIR_32f(x, y, n, pIIRState);
        if (ippsIIR->lrt)
        {
            ippsIIRGetDlyLine_32f(pIIRState, ippsIIR->bufdly);
        }
    }
    else
    {
        //iirdf2_32f(n, x, y, ippsIIR);
        ierr = rtseis_ippsIIRFilter_apply32f_ori(n, x, y, ippsIIR);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Applies the IIR filtering to the array x.
 * @param[in] n            Number of points in x.
 * @param[in] x            Array to filter.  This has dimension [n].
 * @param[out] y           Filtered version of x.  This has dimension [n].
 * @param[in,out] ippsIIR  Holds the IIR filtering states.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFilter_apply64f_ori(const int n, const double x[],
                                      double y[],
                                      struct ippsIIRFilter_struct *ippsIIR)
{
    IppsIIRState_64f *pIIRState = NULL;
    double *bufdly = NULL;
    double *bufdst = NULL;
    double *bufsrc = NULL;
    Ipp64f *pSrc, *pDst;
    IppStatus status;
    int accumLen, ierr, k, lento1024, nptsLoc, nptsFilter, offset;
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    // Verify the structure was initialized.
    if (!ippsIIR->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsIIR not initialized");
        return -1;
    }
    if (ippsIIR->precision != RTSEIS_DOUBLE)
    {
        RTSEIS_ERRMSG("%s", "Error must convert data");
        return -1; 
    }
 //iirdf2(n, x, y, ippsIIR);
//return 0;
    // Get handle to the state pointer
    pIIRState = (IppsIIRState_64f *) ippsIIR->state;
    bufdly = (double *) ippsIIR->bufdly;
    // This is post-processing so just apply the filter with the I.C.'s
    if (!ippsIIR->lrt) //lhaveZI && !lwantZF)
    {
        status = ippsIIRSetDlyLine_64f(pIIRState, bufdly);
        status = ippsIIR_64f(x, y, n, pIIRState);
        if (status != ippStsNoErr)
        {   
            RTSEIS_ERRMSG("%s", "Error filtering");
            return -1; 
        }
        pIIRState = NULL;
        bufdly = NULL;
        return 0;
    }
    // Real-time filtering is more complicated.  This doesn't seem to work
    // anymore
    bufsrc = (double *) ippsIIR->bufsrc;
    bufdst = (double *) ippsIIR->bufdst;
/*
status = ippsIIRSetDlyLine_64f(pIIRState, bufdly);
status = ippsIIR_64f(x, y, n, pIIRState);
status = ippsIIRGetDlyLine_64f(pIIRState, bufdly);
return 0;
*/
    // Filter in chunks of 1024
    if (n <= 1024)
    {
        nptsLoc = n; //MIN(1024, n - k);
        pSrc = (Ipp64f *) x;
        pDst = (Ipp64f *) y;
        // Figure out the filter length and offset (npts%1024)
        nptsFilter = nptsLoc;
        offset = ippsIIR->accumLen%1024;
        if (offset < 0 || offset > 1023)
        {
            RTSEIS_ERRMSG("offset=%d out of bounds [0,1023]", offset);
            return -1;
        }
        // Set the signal into the source array w/ appropriate offset
        status = ippsCopy_64f(x, &bufsrc[offset], nptsFilter);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error setting signal");
            return -1;
        }
        // Set the delay line from the buffer
        status = ippsIIRSetDlyLine_64f(pIIRState, bufdly);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error setting delay line");
            return -1;
        }
        accumLen = offset + nptsFilter; // Number of points accumulated
        lento1024 = 1024 - offset; // This can't exceed 1024
        //printf("(accumlen, offset, n): (%d,%d,%d,%d)\n",
        //       ippsIIR->accumLen, accumLen, offset, nptsFilter);
        // Trigger filtering and update
        if (accumLen >= 1024)
        {
            // Apply filter
            status = ippsIIR_64f(bufsrc, bufdst, 1024, pIIRState);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying filter in update");
                return -1;
            }
            ippsCopy_64f(&bufdst[offset], pDst, lento1024);
            // Get the delay line for the 1024 sample packet
            status = ippsIIRGetDlyLine_64f(pIIRState, bufdly);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error getting 1024 packet delay line");
                return -1;
            }
            // Delay line update has occurred - update counters
            nptsFilter = nptsFilter - lento1024; // Remaining points
            accumLen = accumLen - 1024; // Now this many accumulated points
            offset = 0; // There's no more offset
            pDst = pDst + (size_t) lento1024;
            // Update the source
            if (nptsFilter > 0)
            {
                ippsCopy_64f(&pSrc[lento1024], bufsrc, nptsFilter);
            }
        }
        // If there are remaining points then filter them
        if (accumLen > 0)
        {
            status = ippsIIR_64f(bufsrc, bufdst, accumLen, pIIRState);
        }
        // If there are
        if (nptsFilter > 0)
        {
            ippsCopy_64f(&bufdst[offset], pDst, nptsFilter);
        }
        // Update the number of filtered points
        ippsIIR->accumLen = (ippsIIR->accumLen + nptsLoc)%1024;
    } // Loop on chunks of 1024
    else
    {
        // Apply this function on packets of length 1024 
        for (k=0; k<n; k=k+1024)
        {
            nptsLoc = MIN(1024, n - k);
            ierr = rtseis_ippsIIRFilter_apply64f(nptsLoc, &x[k], &y[k],
                                                 ippsIIR); 
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Error applying filter in chunks");
                return -1;
            }
        }
    }
    // Dereference pointers
    pIIRState = NULL;
    bufdly = NULL; 
    bufdst = NULL;
    bufsrc = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Applies the IIR filtering to the array x.
 * @param[in] n            Number of points in x.
 * @param[in] x            Array to filter.  This has dimension [n].
 * @param[out] y           Filtered version of x.  This has dimension [n].
 * @param[in,out] ippsIIR  Holds the IIR filtering states.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 * @author Ben Baker with substantial help from IPP Dev team.
 */
int rtseis_ippsIIRFilter_apply32f_ori(const int n, const float x[],
                                      float y[],
                                      struct ippsIIRFilter_struct *ippsIIR)
{
    IppsIIRState_32f *pIIRState = NULL;
    float *bufdly = NULL;
    float *bufdst = NULL;
    float *bufsrc = NULL;
    Ipp32f *pSrc, *pDst;
    IppStatus status;
    int accumLen, ierr, k, lento1024, nptsLoc, nptsFilter, offset;
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    // Verify the structure was initialized.
    if (!ippsIIR->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsIIR not initialized");
        return -1;
    }
    if (ippsIIR->precision != RTSEIS_FLOAT)
    {
        RTSEIS_ERRMSG("%s", "Error must convert data");
        return -1;
    }
    // Get handle to the state pointer
    pIIRState = (IppsIIRState_32f *) ippsIIR->state;
    bufdly = (float *) ippsIIR->bufdly;
    // This is easy to do - just apply the filter
    if (!ippsIIR->lrt) //lhaveZI && !lwantZF)
    {
        status = ippsIIRSetDlyLine_32f(pIIRState, bufdly);
        status = ippsIIR_32f(x, y, n, pIIRState);
        if (status != ippStsNoErr)
        {   
            RTSEIS_ERRMSG("%s", "Error filtinering\n");
            return -1; 
        }
        pIIRState = NULL;
        bufdly = NULL;
        return 0;
    }
    // Real-time filtering is more complicated - get handles to pointers
    bufsrc = (float *) ippsIIR->bufsrc;
    bufdst = (float *) ippsIIR->bufdst;
    // Filter in chunks of 1024
    if (n <= 1024)
    {
        nptsLoc = n; //MIN(1024, n - k);
        pSrc = (Ipp32f *) x;
        pDst = (Ipp32f *) y;
        // Figure out the filter length and offset (npts%1024)
        nptsFilter = nptsLoc;
        offset = ippsIIR->accumLen%1024;
        if (offset < 0 || offset > 1023)
        {
            RTSEIS_ERRMSG("offset=%d out of bounds [0,1023]", offset);
            return -1;
        }
        // Set the signal into the source array w/ appropriate offset
        status = ippsCopy_32f(x, &bufsrc[offset], nptsFilter);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error setting signal");
            return -1;
        }
        // Set the delay line from the buffer
        status = ippsIIRSetDlyLine_32f(pIIRState, bufdly);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error setting delay line");
            return -1;
        }
        accumLen = offset + nptsFilter; // Number of points accumulated
        lento1024 = 1024 - offset; // This can't exceed 1024
        //printf("(accumlen, offset, n): (%d,%d,%d,%d)\n",
        //       ippsIIR->accumLen, accumLen, offset, nptsFilter);
        // Trigger filtering and update
        if (accumLen >= 1024)
        {
            // Apply filter
            status = ippsIIR_32f(bufsrc, bufdst, 1024, pIIRState);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying filter in update");
                return -1;
            }
            ippsCopy_32f(&bufdst[offset], pDst, lento1024);
            // Get the delay line for the 1024 sample packet
            status = ippsIIRGetDlyLine_32f(pIIRState, bufdly);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error getting 1024 packet delay line");
                return -1;
            }
            // Delay line update has occurred - update counters
            nptsFilter = nptsFilter - lento1024; // Remaining points
            accumLen = accumLen - 1024; // Now this many accumulated points
            offset = 0; // There's no more offset
            pDst = pDst + (size_t) lento1024;
            // Update the source
            if (nptsFilter > 0)
            {
                ippsCopy_32f(&pSrc[lento1024], bufsrc, nptsFilter);
            }
        }
        // If there are remaining points then filter them
        if (accumLen > 0)
        {
            status = ippsIIR_32f(bufsrc, bufdst, accumLen, pIIRState);
        }
        // If there are
        if (nptsFilter > 0)
        {
            ippsCopy_32f(&bufdst[offset], pDst, nptsFilter);
        }
        // Update the number of filtered points
        ippsIIR->accumLen = (ippsIIR->accumLen + nptsLoc)%1024;
    } // Loop on chunks of 1024
    else
    {
        // Apply this function on packets of length 1024 
        for (k=0; k<n; k=k+1024)
        {
            nptsLoc = MIN(1024, n - k);
            ierr = rtseis_ippsIIRFilter_apply32f(nptsLoc, &x[k], &y[k],
                                                 ippsIIR); 
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Error applying filter in chunks");
                return -1;
            }
        }
    }
    // Dereference pointers
    pIIRState = NULL;
    bufdly = NULL; 
    bufdst = NULL;
    bufsrc = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Releases memory on the IIR filter structure and sets all variables
 *        to 0.
 * @param[out] ippsIIR   On exit all memory has been released and variables
 *                       set to 0.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_iir
 */
int rtseis_ippsIIRFilter_finalize(struct ippsIIRFilter_struct *ippsIIR)
{
    Ipp64f *pTaps64, *zi64;
    Ipp32f *pTaps32, *zi32;
    Ipp8u *pBuf;
    if (!ippsIIR->linit){return 0;} // Never initialized
    ippsFree(ippsIIR->bRef);
    ippsFree(ippsIIR->aRef);
    ippsFree(ippsIIR->bufdst);
    ippsFree(ippsIIR->bufsrc);
    ippsFree(ippsIIR->bufdly);
    if (ippsIIR->precision == RTSEIS_DOUBLE)
    {
        zi64 = (Ipp64f *) ippsIIR->zi;
        pTaps64 = (Ipp64f *) ippsIIR->pTaps;
        pBuf  = (Ipp8u *) ippsIIR->pBuf;
        ippsFree(zi64);
        ippsFree(pTaps64);
        ippsFree(pBuf);
    }
    else if (ippsIIR->precision == RTSEIS_FLOAT)
    {
        zi32 = (Ipp32f *) ippsIIR->zi;
        pTaps32 = (Ipp32f *) ippsIIR->pTaps;
        pBuf  = (Ipp8u *) ippsIIR->pBuf;
        ippsFree(zi32);
        ippsFree(pTaps32);
        ippsFree(pBuf);
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    memset(ippsIIR, 0, sizeof(struct ippsIIRFilter_struct));
    return 0;
}

#if IPP_VERSION_MAJOR < 2019
static int iirdf2_64f(const int n, const double x[],
                      double y[],
                      struct ippsIIRFilter_struct *iir)
{
    double Xi, Yi;
    int i, j;
    int order = iir->order; 
    Ipp64f *b  = iir->bRef; //&pTaps[0];
    Ipp64f *a  = iir->aRef; //&pTaps[order+1];
    Ipp64f *vi = iir->bufdly;
    Ipp64f *v  = iir->bufdst;
    // Loop on samples
    for (i=0; i<n; i++)
    {
        #pragma omp simd
        for (j=order; j>=1; j--){v[j] = vi[j-1];}
        Xi = 0;
        Yi = 0;
        #pragma omp simd reduction(+:Xi,Yi)
        for (j=1; j<=order; j++)
        {
            Xi = Xi - a[j]*v[j];
            Yi = Yi + b[j]*v[j];
        }
        Xi   = Xi + x[i];
        v[0] = Xi;
        y[i] = Yi + b[0]*v[0];
        #pragma omp simd
        for (j=0; j<=order; j++){vi[j] = v[j];}
    }
    if (!iir->lrt)
    {   
        #pragma omp simd
        for (j=0; j<=order; j++){vi[j] = 0;} 
    }
    return 0;
}

static int iirdf2_32f(const int n, const float x[],
                      float y[],
                      struct ippsIIRFilter_struct *iir)
{
    double Xi, Yi;
    int i, j; 
    int order = iir->order; 
    Ipp64f *b  = iir->bRef; //&pTaps[0];
    Ipp64f *a  = iir->aRef; //&pTaps[order+1];
    Ipp64f *vi = iir->bufdly;
    Ipp64f *v  = iir->bufdst;
    // Loop on samples
    for (i=0; i<n; i++)
    {
        #pragma omp simd
        for (j=order; j>=1; j--){v[j] = vi[j-1];}
        Xi = 0;
        #pragma omp simd reduction(+:Xi)
        for (j=1; j<=order; j++)
        {   
            Xi = Xi - a[j]*v[j];
        }
        Xi   = Xi + (double) x[i];
        v[0] = Xi;
        Yi = 0;
        #pragma omp simd reduction(+:Yi)
        for (j=0; j<=order; j++)
        {   
            Yi = Yi + b[j]*v[j];
        }
        y[i] = (float) Yi;
        v[0] = Xi;
        #pragma omp simd
        for (j=0; j<=order; j++){vi[j] = v[j];}
    }
    if (!iir->lrt)
    {   
        #pragma omp simd
        for (j=0; j<=order; j++){vi[j] = 0;}
    }
    return 0;
}
#endif
/*
  // This is the working code from ispl
 for (i=0; i<npts; i++){
           // This is the working code from ispl
            for (j=order; j>=1; j--)
            {
                v[j] = v[j-1];
            }
            Xi = 0.0;
            Yi = 0.0;
            #pragma omp simd reduction(+:Xi, Yi)
            for (j=order; j>=1; j--)
            {
                Xi = Xi - a[j]*v[j];
                Yi = Yi + b[j]*v[j];
            }
            Xi   = Xi + x[i];
            y[i] = Yi + b[0]*Xi;
            v[0] = Xi;
      }
*/
