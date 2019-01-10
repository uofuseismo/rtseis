#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/ippsHelper.h"
#include <ipps.h>

/*!
 * @defgroup rtseis_ipps_dft DFT
 * @brief Interface to the IPP discrete Fourier transform.
 *        This is from ISTI's ISPL but function names and enums
 *        have names have changed to conform to RTSeis.
 * @ingroup rtseis_ipps
 * @author Ben Baker
 * @copyright ISTI distributed under the Apache 2 license.
 */
/*!
 * @brief Initializes the real to complex discrete Fourier transform structure.
 * @param[in] luseFFT      If true then this will force the base 2 FFT to be 
 *                         used.
 * @param[in] luseFFT      If false then the DFT or FFT will be selected 
 *                         according to the the length of the input signal.
 * @param[in] length       Length of signal to transform.  If this is a power of
 *                         2 then the FFT will be used instead of the DFT.
 * @param[in] precision    Precision of the transform.  This an be RTSEIS_Float
 *                         or RTSEIS_Double
 *
 * @param[out] ippsDFTR2C  The initialized Fourier transform structure.
 * @result 0 indicates success.
 * @ingroup rtseis_ipp_dft
 */
int rtseis_ippsDFTR2C_initialize(const bool luseFFT,
                                 const int length,
                                 const enum rtseisPrecision_enum precision,
                                 struct ippsDFTR2C_struct *ippsDFTR2C)
{
    IppStatus status;
    Ipp8u *pDFTInitBuf, *pFFTSpecBuf;
    IppsFFTSpec_R_64f *pFFTSpec64 = NULL;
    IppsDFTSpec_R_64f *pDFTSpec64 = NULL;
    IppsFFTSpec_R_32f *pFFTSpec32 = NULL;
    IppsDFTSpec_R_32f *pDFTSpec32 = NULL;
    int n2, order, sizeInit;
    if (length < 2)
    {
        RTSEIS_ERRMSG("Length=%d must be at least 2", length);
        return -1;
    }
    memset(ippsDFTR2C, 0, sizeof(struct ippsDFTR2C_struct));
    // Decide between FFT and DFT 
    ippsDFTR2C->order =-1;
    ippsDFTR2C->length = length;
    order = (int) (round(log2((double) length)));
    n2 = (int) pow(2, order);
    if (luseFFT || n2 == length)
    {
        // Forcing the FFT
        if (n2 != length)
        {
            // Pad 
            if (n2 < length)
            {
                order = order + 1;
                n2 = (int) pow(2, order);
            }
            if (n2 < length)
            {
                RTSEIS_ERRMSG("%s", "Algorithmic failure");
                return -1;
            }
        }
        ippsDFTR2C->luseFFT = true;
        ippsDFTR2C->order = order;
        ippsDFTR2C->length = n2;
    }
    ippsDFTR2C->lenft = ippsDFTR2C->length/2 + 1;
    // Initialize the structures 
    if (precision == RTSEIS_DOUBLE)
    {
        if (ippsDFTR2C->luseFFT)
        {
            status = ippsFFTGetSize_R_64f(ippsDFTR2C->order,
                                          IPP_FFT_DIV_INV_BY_N,
                                          ippAlgHintNone,
                                          &ippsDFTR2C->specSize,
                                          &sizeInit,
                                          &ippsDFTR2C->bufferSize);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                return -1; 
            }
            pFFTSpecBuf = ippsMalloc_8u(ippsDFTR2C->specSize); 
            pDFTInitBuf = ippsMalloc_8u(sizeInit);
            status = ippsFFTInit_R_64f(&pFFTSpec64, 
                                       ippsDFTR2C->order, 
                                       IPP_FFT_DIV_INV_BY_N, 
                                       ippAlgHintNone,
                                       pFFTSpecBuf,
                                       pDFTInitBuf);
            if (pDFTInitBuf){ippsFree(pDFTInitBuf);}
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                return -1;
            }
            ippsDFTR2C->fftBuffer = (void *) pFFTSpecBuf;
            ippsDFTR2C->ftHandle = (void *) pFFTSpec64;
        }
        else
        {
            status = ippsDFTGetSize_R_64f(ippsDFTR2C->length, 
                                          IPP_FFT_DIV_INV_BY_N,
                                          ippAlgHintNone,
                                          &ippsDFTR2C->specSize,
                                          &sizeInit,
                                          &ippsDFTR2C->bufferSize);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                return -1;
            }
            pDFTSpec64
                = (IppsDFTSpec_R_64f *) ippsMalloc_8u(ippsDFTR2C->specSize);
            pDFTInitBuf = ippsMalloc_8u(sizeInit);
            status = ippsDFTInit_R_64f(ippsDFTR2C->length,
                                       IPP_FFT_DIV_INV_BY_N,
                                       ippAlgHintNone,
                                       pDFTSpec64,
                                       pDFTInitBuf);
            if (pDFTInitBuf){ippsFree(pDFTInitBuf);}
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                return -1; 
            }
            ippsDFTR2C->ftHandle = (void *) pDFTSpec64;
        }
        ippsDFTR2C->pBuf  = (void *) ippsMalloc_8u(ippsDFTR2C->bufferSize);
        ippsDFTR2C->xwork = (void *) ippsMalloc_64f(ippsDFTR2C->length);
        ippsDFTR2C->ywork = (void *) ippsMalloc_64f(2*ippsDFTR2C->length);
        ippsZero_64f(ippsDFTR2C->xwork,   ippsDFTR2C->length);
        ippsZero_64f(ippsDFTR2C->ywork, 2*ippsDFTR2C->length);
    }
    else if (precision == RTSEIS_FLOAT)
    {
        if (ippsDFTR2C->luseFFT)
        {
            status = ippsFFTGetSize_R_32f(ippsDFTR2C->order,
                                          IPP_FFT_DIV_INV_BY_N,
                                          ippAlgHintNone,
                                          &ippsDFTR2C->specSize,
                                          &sizeInit,
                                          &ippsDFTR2C->bufferSize);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                return -1; 
            }
            pFFTSpecBuf = ippsMalloc_8u(ippsDFTR2C->specSize); 
            pDFTInitBuf = ippsMalloc_8u(sizeInit);
            status = ippsFFTInit_R_32f(&pFFTSpec32, 
                                       ippsDFTR2C->order, 
                                       IPP_FFT_DIV_INV_BY_N, 
                                       ippAlgHintNone,
                                       pFFTSpecBuf,
                                       pDFTInitBuf);
            if (pDFTInitBuf){ippsFree(pDFTInitBuf);}
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                return -1; 
            }
            ippsDFTR2C->fftBuffer = (void *) pFFTSpecBuf;
            ippsDFTR2C->ftHandle = (void *) pFFTSpec32;
        }
        else
        {
            status = ippsDFTGetSize_R_32f(ippsDFTR2C->length, 
                                          IPP_FFT_DIV_INV_BY_N,
                                          ippAlgHintNone,
                                          &ippsDFTR2C->specSize,
                                          &sizeInit,
                                          &ippsDFTR2C->bufferSize);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                return -1; 
            }
            pDFTSpec32
                = (IppsDFTSpec_R_32f *) ippsMalloc_8u(ippsDFTR2C->specSize);
            pDFTInitBuf = ippsMalloc_8u(sizeInit);
            status = ippsDFTInit_R_32f(ippsDFTR2C->length,
                                       IPP_FFT_DIV_INV_BY_N,
                                       ippAlgHintNone,
                                       pDFTSpec32,
                                       pDFTInitBuf);
            if (pDFTInitBuf){ippsFree(pDFTInitBuf);}
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                return -1; 
            }
            ippsDFTR2C->ftHandle = (void *) pDFTSpec32;
        }
        ippsDFTR2C->pBuf  = (void *) ippsMalloc_8u(ippsDFTR2C->bufferSize);
        ippsDFTR2C->xwork = (void *) ippsMalloc_32f(ippsDFTR2C->length);
        ippsDFTR2C->ywork = (void *) ippsMalloc_32f(2*ippsDFTR2C->length);
        ippsZero_32f(ippsDFTR2C->xwork,   ippsDFTR2C->length);
        ippsZero_32f(ippsDFTR2C->ywork, 2*ippsDFTR2C->length);
    }
    else
    {
        RTSEIS_ERRMSG("Invalid precision %d", (int) precision);
        return -1;
    }
    ippsDFTR2C->precision = precision;
    ippsDFTR2C->linit = true;
    return 0;
}
//============================================================================//
/*!
 * @brief Releases the memory on the ippsDFT real-to-complex Fourier transform
 *        structure.
 * @param[in,out] ippsDFTR2C  On input this is the initialized structure.
 * @param[in,out] ippsDFTR2C  On exit all memory has been released.
 * @result 0 indicates success.
 * @ingroup rtseis_ipps_dft
 */
int rtseis_ippsDFTR2C_finalize(struct ippsDFTR2C_struct *ippsDFTR2C)
{
    if (!ippsDFTR2C->linit){return 0;}
    if (ippsDFTR2C->pBuf != NULL){ippsFree(ippsDFTR2C->pBuf);}
    if (ippsDFTR2C->xwork != NULL){ippsFree(ippsDFTR2C->xwork);}
    if (ippsDFTR2C->ywork != NULL){ippsFree(ippsDFTR2C->ywork);}
    if (ippsDFTR2C->luseFFT)
    {
        if (ippsDFTR2C->fftBuffer != NULL){ippsFree(ippsDFTR2C->fftBuffer);}
    }
    else
    {
        if (ippsDFTR2C->ftHandle != NULL){ippsFree(ippsDFTR2C->ftHandle);}
    }
    memset(ippsDFTR2C, 0, sizeof(struct ippsDFTR2C_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Fourier transforms a real signal. 
 *
 * @param[in] n      The length of the signal to transform.  This must be 
 *                   0 < n <= ippsDFTR2C->length.  If n < ippsDFTR2C->length
 *                   then x will be 0 padded.
 * @param[in] x      Signal to transform.  This is an array of dimension [n].
 * @param[in] maxy   Max space allocated to y.  If -1 then this is a space 
 *                   inquiry.
 *
 * @param[out] y     The Fourier transformed version of x.  This is an array
 *                   of dimension [maxy].
 * @param[out] ny    Length of signal y.
 *
 * @param[in,out] ippsDFTR2C  On input this is the initialized DFT structure.
 * @param[in,out] ippsDFTR2C  On exit the workspace arrays have been updated.
 *
 * @ingroup rtseis_ipps_dft
 *
 */
int rtseis_ippsDFTR2C_applyForwardDFT64f(const int n,
                                         const double x[],
                                         const int maxy,
                                         int *ny,
                                         void *yptr,
                                         struct ippsDFTR2C_struct *ippsDFTR2C)
{
    double complex *y = (double complex *) yptr;
    double *x64;
    float *x32;
    float complex *y32;
    IppsFFTSpec_R_64f *pFFTSpec;
    IppsDFTSpec_R_64f *pDFTSpec;
    const Ipp64f *pSrc;
    Ipp64f *pDst;
    Ipp8u *pBuf;
    IppStatus status;
    int i, ierr;
    *ny = 0;
    if (n < 0)
    {
        RTSEIS_ERRMSG("No points to transform: %d", n);
        return -1;
    }
    if (!ippsDFTR2C->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsDFTR2C not initialized");
        return -1;
    }
    *ny = ippsDFTR2C->lenft; //ippsDFTR2C->length/2 + 1;
    // Space inquiry
    if (maxy < 0){return 0;}
    // Handle the other precision
    if (ippsDFTR2C->precision == RTSEIS_FLOAT)
    {
        x32 = ippsMalloc_32f(n);
        ippsConvert_64f32f(x, x32, n);
        y32 = (float complex *) ippsMalloc_32f(*ny*2);
        ierr = rtseis_ippsDFTR2C_applyForwardDFT32f(n, x32, maxy, ny,
                                                    y32, ippsDFTR2C);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Error applying forward DFT");
        }
        else
        {
            #pragma omp simd aligned(y32:64)
            for (i=0; i<*ny; i++)
            {
                y[i] = (double complex) y32[i];
            }
        }
        ippsFree(x32);
        ippsFree(y32); 
        return ierr;
    }
    if (ippsDFTR2C->precision != RTSEIS_DOUBLE)
    {
        RTSEIS_ERRMSG("%s", "Invalid precision"); 
        return -1;
    }
    // Need to pad
    if (n < ippsDFTR2C->length)
    {
        x64 = (double *) ippsDFTR2C->xwork;
        ippsCopy_64f(x, x64, n);
        status = ippsZero_64f(&x64[n], ippsDFTR2C->length - n);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error padding signal");
            return -1;
        }
        pSrc = x64;
        pDst = (Ipp64f *) y;
    }
    // Otherwise just set the pointer
    else
    {
        pSrc = x;
        pDst = (Ipp64f *) y;
    }
    pBuf = (Ipp8u *) ippsDFTR2C->pBuf;
    // Apply the transform 
    if (ippsDFTR2C->luseFFT)
    {
        pFFTSpec = (IppsFFTSpec_R_64f *) ippsDFTR2C->ftHandle;
        status = ippsFFTFwd_RToCCS_64f(pSrc, pDst, pFFTSpec, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error applying FFT");
            return -1;
        } 
        pFFTSpec= NULL;
    }
    else
    {
        pDFTSpec = (IppsDFTSpec_R_64f *) ippsDFTR2C->ftHandle;
        status = ippsDFTFwd_RToCCS_64f(pSrc, pDst, pDFTSpec, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error applying DFT");
            return -1;
        }
        pDFTSpec = NULL;
    }
    pBuf = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Inverse transforms a signal back to the time domain.
 * @param[in] lenft  Length of the signal x.
 * @param[in] x      Signal to inverse transform.  This has array of
 *                   dimension [n].
 * @param[in] ny     Number of points in output signal.  This should equal
 *                   ippsDFTR2C.length.
 * @param[out] y     The inverse Fourier transformed signal.  This has dimension
 *                   [ny].
 *
 * @param[in,out] ippsDFTR2C  On input contains the initialized DFT structure.
 * @param[in,out] ippsDFTR2C  On exit workspace arrays have been modified.
 *
 * @result 0 indicates success
 *
 * @ingroup rtseis_ipps_dft
 *
 */
int rtseis_ippsDFTR2C_applyInverseDFT64f(const int lenft,
                                         const void *xptr,
                                         const int ny,
                                         double y[],
                                         struct ippsDFTR2C_struct *ippsDFTR2C)
{
    double complex *x = (double complex *) xptr;
    double complex *x64;
    float complex *x32;
    float *y32;
    IppsFFTSpec_R_64f *pFFTSpec;
    IppsDFTSpec_R_64f *pDFTSpec;
    const Ipp64f *pSrc;
    Ipp64f *pDst;
    Ipp8u *pBuf;
    IppStatus status;
    int i, ierr;
    ierr = 0;
    if (lenft < 0)
    {
        RTSEIS_ERRMSG("No points to transform: %d", lenft);
        return -1;
    }
    if (!ippsDFTR2C->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsDFTR2C not initialized");
        return -1;
    }
    if (ny > ippsDFTR2C->length)
    {
        RTSEIS_ERRMSG("Output length=%d exceeds limit=%d",
                    ny, ippsDFTR2C->length);
        return -1;
    }
    // Handle the other precision
    if (ippsDFTR2C->precision == RTSEIS_FLOAT)
    {
        x32 = (float complex *) ippsMalloc_32f(2*lenft);
        #pragma omp simd aligned(x32:64)
        for (i=0; i<lenft; i++)
        {
            x32[i] = (float complex) x[i];
        }
        //ippsConvert_64fc32fc(x, x32, lenft);
        y32 = (float *) ippsMalloc_32f(ny);
/*
        ierr = rtseis_ippsDFTR2C_applyForwardDFT32f(n, x32, ny,
                                                    y32, ippsDFTR2C);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Error applying forward DFT");
        }
        else
*/
        {
            ippsConvert_32f64f(y32, y, ny);
        }
        ippsFree(x32);
        ippsFree(y32);
        return ierr;
    }
    if (ippsDFTR2C->precision != RTSEIS_DOUBLE)
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    // Need to pad in Fourier domain
    if (lenft < ippsDFTR2C->lenft)
    {
        x64 = (double complex *) ippsDFTR2C->ywork;
        ippsCopy_64fc((const Ipp64fc *) x, (Ipp64fc *) x64, lenft);
        status = ippsZero_64fc((Ipp64fc *) &x64[lenft],
                               ippsDFTR2C->lenft - lenft);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error padding signal");
            return -1;
        }
        pSrc = (Ipp64f *) x64;
    }
    // Otherwise just set the pointer
    else
    {
        pSrc = (Ipp64f *) x;
    }
    // Need to truncate
    pDst = (Ipp64f *) y;
    if (ny < ippsDFTR2C->length){pDst = (Ipp64f *) ippsDFTR2C->xwork;}
    // Apply the transform
    pBuf = (Ipp8u *) ippsDFTR2C->pBuf;
    if (ippsDFTR2C->luseFFT)
    {
        pFFTSpec = (IppsFFTSpec_R_64f *) ippsDFTR2C->ftHandle;
        status = ippsFFTInv_CCSToR_64f(pSrc, pDst, pFFTSpec, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error applying FFT");
            return -1;
        }
        pFFTSpec= NULL;
    }
    else
    {
        pDFTSpec = (IppsDFTSpec_R_64f *) ippsDFTR2C->ftHandle;
        status = ippsDFTInv_CCSToR_64f(pSrc, pDst, pDFTSpec, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error applying DFT");
            return -1;
        }
        pDFTSpec = NULL;
    }
    // Truncate
    if (ny < ippsDFTR2C->length)
    {
        ippsCopy_64f(pDst, y, ny);
        pDst = NULL;
    }
    pBuf = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Fourier transforms a real signal.
 *
 * @param[in] n      The length of the signal to transform.  This must be 
 *                   0 < n <= ippsDFTR2C->length.  If n < ippsDFTR2C->length
 *                   then x will be 0 padded.
 * @param[in] x      Signal to transform.  This is an array of dimension [n].
 * @param[in] maxy   Max space allocated to y.  If -1 then this is a space 
 *                   inquiry.
 *
 * @param[out] y     The Fourier transformed version of x.  This is an array
 *                   of dimension [maxy].
 * @param[out] ny    Length of signal y.
 *
 * @param[in,out] ippsDFTR2C  On input this is the initialized DFT structure.
 * @param[in,out] ippsDFTR2C  On exit the workspace arrays have been updated.
 *
 * @copyright ISTI distributed under the Apache 2 license.
 *
 * @ingroup rtseis_ipps_dft
 *
 */
int rtseis_ippsDFTR2C_applyForwardDFT32f(const int n,
                                         const float x[],
                                         const int maxy,
                                         int *ny,
                                         void *yptr,
                                         struct ippsDFTR2C_struct *ippsDFTR2C)
{
    float complex *y = (float complex *) yptr;
    double *x64;
    float *x32;
    double complex *y64;
    IppsFFTSpec_R_32f *pFFTSpec;
    IppsDFTSpec_R_32f *pDFTSpec;
    const Ipp32f *pSrc;
    Ipp32f *pDst;
    Ipp8u *pBuf;
    IppStatus status;
    int i, ierr;
    *ny = 0;
    if (n < 0)
    {
        RTSEIS_ERRMSG("No points to transform: %d", n);
        return -1;
    }
    if (!ippsDFTR2C->linit)
    {
        RTSEIS_ERRMSG("%s", "ippsDFTR2C not initialized");
        return -1;
    }
    *ny = ippsDFTR2C->lenft; //ippsDFTR2C->length/2 + 1;
    // Space inquiry
    if (maxy < 0){return 0;}
    // Handle the other precision
    if (ippsDFTR2C->precision == RTSEIS_FLOAT)
    {
        x64 = ippsMalloc_64f(n);
        ippsConvert_32f64f(x, x64, n);
        y64 = (double complex *) ippsMalloc_64f(*ny*2);
        ierr = rtseis_ippsDFTR2C_applyForwardDFT64f(n, x64, maxy, ny,
                                                    y64, ippsDFTR2C);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Error calling DFT forward");
        }
        else
        {
            #pragma omp simd aligned(y64:64)
            for (i=0; i<*ny; i++)
            {
                y[i] = (float complex) y64[i];
            }
        }
        ippsFree(x64);
        ippsFree(y64);
        return ierr;
    }
    if (ippsDFTR2C->precision != RTSEIS_DOUBLE)
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    // Need to pad
    if (n < ippsDFTR2C->length)
    {
        x32 = (float *) ippsDFTR2C->xwork;
        ippsCopy_32f(x, x32, n);
        status = ippsZero_32f(&x32[n], ippsDFTR2C->length - n);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error padding signal");
            return -1;
        }
        pSrc = x32;
        pDst = (Ipp32f *) y;
    }
    // Otherwise just set the pointer
    else
    {
        pSrc = x;
        pDst = (Ipp32f *) y;
    }
    pBuf = (Ipp8u *) ippsDFTR2C->pBuf;
    // Apply the transform 
    if (ippsDFTR2C->luseFFT)
    {
        pFFTSpec = (IppsFFTSpec_R_32f *) ippsDFTR2C->ftHandle;
        status = ippsFFTFwd_RToCCS_32f(pSrc, pDst, pFFTSpec, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error applying FFT");
            return -1;
        }
        pFFTSpec= NULL;
    }
    else
    {
        pDFTSpec = (IppsDFTSpec_R_32f *) ippsDFTR2C->ftHandle;
        status = ippsDFTFwd_RToCCS_32f(pSrc, pDst, pDFTSpec, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error applying DFT");
            return -1;
        }
        pDFTSpec = NULL;
    }
    pBuf = NULL;
    return 0;
}
