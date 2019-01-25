#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#define RTSEIS_LOGGING 1
#include "rtseis/utils/transforms.hpp"
#include "rtseis/log.h"
#include <ipps.h>

using namespace RTSeis::Utils::Transforms;

DFTR2C::DFTR2C(void)
{
    clear();
    return;
}

DFTR2C::DFTR2C(const DFTR2C &dftr2c)
{
    *this = dftr2c;
    return;
}

DFTR2C& DFTR2C::operator=(const DFTR2C &dftr2c)
{
    if (&dftr2c == this){return *this;}
    clear();
    if (!dftr2c.isInitialized()){return *this;}
    int ierr = initialize(dftr2c.length_, dftr2c.ldoFFT_,
                          dftr2c.precision_);
    if (ierr != 0)
    {
        clear();
        return *this;
    }
    if (precision_ == RTSeis::Precision::DOUBLE)
    {

    }
    else
    {

    }
    return *this;
}

void DFTR2C::clear(void)
{
    if (ftHandle_ != nullptr){ippsFree(ftHandle_);}
    if (ftBuffer_ != nullptr){ippsFree(ftBuffer_);}
    if (work_ != nullptr){ippsFree(work_);}
    ftHandle_ = nullptr;
    ftBuffer_ = nullptr;
    work_ = nullptr;
    precision_ = RTSeis::Precision::DOUBLE;
    length_ = 0;
    lenft_ = 0;
    nwork_ = 0;
    bufferSize_ = 0; 
    specSize_ = 0;
    order_ = 0;
    ldoFFT_ = false;
    isInitialized_ = false;
    return;
}


int DFTR2C::initialize(const int length,
                       const bool ldoFFT,
                       const RTSeis::Precision precision)
{
    clear();
    // Check the inputs
    if (length < 2)
    {
        if (length < 2){RTSEIS_ERRMSG("Length=%d must be at least 2", length);}
        return -1;
    }
    // Decide between DFT and FFT
    order_ =-1;
    length_ = length;
    double dlen = static_cast<double> (length_);
    int orderWork = static_cast<int> (std::round(std::log2(dlen)));
    int n2 = static_cast<int> (std::pow(2, orderWork));
    if (ldoFFT || n2 == length)
    {
        // Forcing the FFT
        if (n2 != length)
        {
            // Pad 
            if (n2 < length)
            {
                orderWork = orderWork + 1;
                n2 = static_cast<int> (std::pow(2, orderWork));
            }
            if (n2 < length)
            {
                RTSEIS_ERRMSG("%s", "Algorithmic failure");
                clear();
                return -1;
            }
        }
        ldoFFT_ = true;
        order_ = orderWork;
        length_ = n2;
    }
    // Initialize the appropriate transform
    IppStatus status;
    int sizeInit;
    if (precision == RTSeis::Precision::DOUBLE)
    {
        if (ldoFFT_)
        {
            status = ippsFFTGetSize_R_64f(order_,
                                          IPP_FFT_DIV_INV_BY_N,
                                          ippAlgHintNone,
                                          &specSize_,
                                          &sizeInit,
                                          &bufferSize_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                clear();
                return -1;
            }
            IppsFFTSpec_R_64f *pFFTSpec = NULL;
            Ipp8u *pSpecBuffer = ippsMalloc_8u(sizeInit);
            Ipp8u *pSpec = ippsMalloc_8u(specSize_);
            status = ippsFFTInit_R_64f(&pFFTSpec,
                                       order_,
                                       IPP_FFT_DIV_INV_BY_N,
                                       ippAlgHintNone,
                                       pSpec,
                                       pSpecBuffer);
            if (pSpecBuffer){ippsFree(pSpecBuffer);}
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to initialize FFT");
                ippsFree(pSpec);
                clear();
                return -1;
            }
            ftHandle_ = pFFTSpec;
            ftBuffer_ = ippsMalloc_8u(bufferSize_);
        }
        else
        {
            status = ippsDFTGetSize_R_64f(length_,
                                          IPP_FFT_DIV_INV_BY_N,
                                          ippAlgHintNone,
                                          &specSize_,
                                          &sizeInit,
                                          &bufferSize_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                clear();
                return -1;
            }
            Ipp8u *pSpecBuffer = ippsMalloc_8u(sizeInit);
            IppsDFTSpec_R_64f *pSpec
             = reinterpret_cast<IppsDFTSpec_R_64f *> (ippsMalloc_8u(specSize_));
            status = ippsDFTInit_R_64f(length_,
                                       IPP_FFT_DIV_INV_BY_N,
                                       ippAlgHintNone,
                                       pSpec,
                                       pSpecBuffer);
            if (pSpecBuffer){ippsFree(pSpecBuffer);}
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to initialize DFT");
                clear();
                return -1;
            }
            ftHandle_ = pSpecBuffer;
            ftBuffer_ = ippsMalloc_8u(bufferSize_);
        }
        nwork_ = std::max(length_, 2*lenft_);
        Ipp64f *work = ippsMalloc_64f(nwork_);
        ippsZero_64f(work, nwork_);
        work_ = work;
    }
    else
    {
        if (ldoFFT_)
        {
            status = ippsFFTGetSize_R_32f(order_,
                                          IPP_FFT_DIV_INV_BY_N,
                                          ippAlgHintNone,
                                          &specSize_,
                                          &sizeInit,
                                          &bufferSize_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                clear();
                return -1;
            }
            IppsFFTSpec_R_32f *pFFTSpec = NULL;
            Ipp8u *pSpecBuffer = ippsMalloc_8u(sizeInit);
            Ipp8u *pSpec = ippsMalloc_8u(specSize_);
            status = ippsFFTInit_R_32f(&pFFTSpec,
                                       order_,
                                       IPP_FFT_DIV_INV_BY_N,
                                       ippAlgHintNone,
                                       pSpec,
                                       pSpecBuffer);
            if (pSpecBuffer){ippsFree(pSpecBuffer);}
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to initialize FFT");
                ippsFree(pSpec);
                clear();
                return -1;
            }
            ftHandle_ = pFFTSpec;
            ftBuffer_ = ippsMalloc_8u(bufferSize_);
        }
        else
        {
            status = ippsDFTGetSize_R_32f(length_,
                                          IPP_FFT_DIV_INV_BY_N,
                                          ippAlgHintNone,
                                          &specSize_,
                                          &sizeInit,
                                          &bufferSize_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get buffer sizes");
                clear();
                return -1;
            }
            Ipp8u *pSpecBuffer = ippsMalloc_8u(sizeInit);
            IppsDFTSpec_R_32f *pSpec
             = reinterpret_cast<IppsDFTSpec_R_32f *> (ippsMalloc_8u(specSize_));
            status = ippsDFTInit_R_32f(length_,
                                       IPP_FFT_DIV_INV_BY_N,
                                       ippAlgHintNone,
                                       pSpec,
                                       pSpecBuffer);
            if (pSpecBuffer){ippsFree(pSpecBuffer);}
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to initialize DFT");
                clear();
                return -1;
            }
            ftHandle_ = pSpecBuffer;
            ftBuffer_ = ippsMalloc_8u(bufferSize_);
        }
        nwork_ = std::max(length_, 2*lenft_);
        Ipp32f *work = ippsMalloc_32f(nwork_);
        ippsZero_32f(work, nwork_);
        work_ = work;
    }
    lenft_ = length/2 + 1;
    precision_ = precision;
    isInitialized_ = true;
    return 0;
}

int DFTR2C::inverseTransform(const int lenft,
                             const std::complex<double> x[],
                             const int maxy, double y[])
{
    if (!isInitialized_)
    {   
        RTSEIS_ERRMSG("%s", "Class is not intiialized");
        return -1; 
    }   
    if (x == nullptr || y == nullptr)
    {   
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1;
    }   
    if (lenft > lenft_ || maxy < length_)
    {
        if (maxy < length_)
        {
            RTSEIS_ERRMSG("maxy = %d must be at least %d", maxy, length_);
        }
        if (lenft > lenft_)
        {
            RTSEIS_ERRMSG("lenft = %d cannot exceed %d", lenft, lenft_);
        }
        return -1;
    }
    // In this case the entire signal would be zero-padded
    if (lenft <= 0)
    {
        Ipp64f *pDst = static_cast<Ipp64f *> (y);
        ippsZero_64f(pDst, length_); 
        return 0;
    }
    // Get handles
    Ipp64f *pDst = static_cast<Ipp64f *> (y);
    Ipp8u *pBuf = static_cast<Ipp8u *> (ftBuffer_);
    // No need to zero-pad
    if (lenft == lenft_)
    {
        // Get handle to data
        const Ipp64f *pSrc = reinterpret_cast<const Ipp64f *> (x);
        // FFT
        if (ldoFFT_)
        {
            IppsFFTSpec_R_64f *pSpec
                = static_cast<IppsFFTSpec_R_64f *> (ftHandle_);
            IppStatus status = ippsFFTInv_CCSToR_64f(pSrc, pDst, pSpec, pBuf);
            if (status != ippStsNoErr)
            {   
                RTSEIS_ERRMSG("%s", "Error applying FFT");
                return -1; 
            }       
        }
        // DFT
        else
        {
            IppsDFTSpec_R_64f *pSpec
                = static_cast<IppsDFTSpec_R_64f *> (ftHandle_);
            IppStatus status = ippsDFTInv_CCSToR_64f(pSrc, pDst, pSpec, pBuf);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying DFT");
                return -1;
            }
        }
    }
    else
    {
        // Copy data and zero pad
        const Ipp64f *xTemp = reinterpret_cast<const Ipp64f *> (x);
        Ipp64f *pSrc = static_cast<Ipp64f *> (work_);
        ippsCopy_64f(xTemp, pSrc, 2*lenft);
        ippsZero_64f(&pSrc[2*lenft], 2*(lenft_-lenft));
        // FFT
        if (ldoFFT_)
        {
            IppsFFTSpec_R_64f *pSpec
                = static_cast<IppsFFTSpec_R_64f *> (ftHandle_);
            IppStatus status = ippsFFTInv_CCSToR_64f(pSrc, pDst, pSpec, pBuf);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying FFT");
                return -1;
            }
        }
        // DFT
        else
        {
            IppsDFTSpec_R_64f *pSpec
                = static_cast<IppsDFTSpec_R_64f *> (ftHandle_);
            IppStatus status = ippsDFTInv_CCSToR_64f(pSrc, pDst, pSpec, pBuf);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying DFT");
                return -1;
            }
        }
    }
    return 0;
} 

int DFTR2C::forwardTransform(const int n, const double x[],
                             const int maxy, std::complex<double> y[])
{
    if (!isInitialized_)
    {
        RTSEIS_ERRMSG("%s", "Class is not intiialized");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1;
    }
    if (maxy < lenft_ || n > length_)
    {
        if (maxy < lenft_)
        {
            RTSEIS_ERRMSG("maxy = %d must be at least %d", maxy, lenft_);
        }
        if (n > length_)
        {
            RTSEIS_ERRMSG("n = %d cannot exceed %d", n, length_);
        }
        return -1;
    }
    // In this case the entire signal would be zero-padded
    if (n <= 0)
    {
        Ipp64f *pDst = reinterpret_cast<Ipp64f *> (y);
        ippsZero_64f(pDst, 2*lenft_); 
        return 0;
    }
    // Get handles
    Ipp64f *pDst = reinterpret_cast<Ipp64f *> (y);
    Ipp8u *pBuf = static_cast<Ipp8u *> (ftBuffer_);
    // No need to zero-pad
    if (n == length_)
    {
        // Get handle to data
        const Ipp64f *pSrc = static_cast<const Ipp64f *> (x);
        // FFT
        if (ldoFFT_)
        {
            IppsFFTSpec_R_64f *pSpec
                = static_cast<IppsFFTSpec_R_64f *> (ftHandle_);
            IppStatus status = ippsFFTFwd_RToCCS_64f(pSrc, pDst, pSpec, pBuf);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying FFT");
                return -1; 
            }       
        }
        // DFT
        else
        {
            IppsDFTSpec_R_64f *pSpec
                = static_cast<IppsDFTSpec_R_64f *> (ftHandle_);
            IppStatus status = ippsDFTFwd_RToCCS_64f(pSrc, pDst, pSpec, pBuf);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying DFT");
                return -1;
            }
        }
    }
    // Zero-pad
    else
    {
        // Copy data and zero pad
        Ipp64f *pSrc = static_cast<Ipp64f *> (work_);
        ippsCopy_64f(x, pSrc, n);
        ippsZero_64f(&pSrc[n], length_-n);
        // FFT
        if (ldoFFT_)
        {
            IppsFFTSpec_R_64f *pSpec
                = static_cast<IppsFFTSpec_R_64f *> (ftHandle_);
            IppStatus status = ippsFFTFwd_RToCCS_64f(pSrc, pDst, pSpec, pBuf);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying FFT");
                return -1;
            }
        }
        // DFT
        else
        {
            IppsDFTSpec_R_64f *pSpec
                = static_cast<IppsDFTSpec_R_64f *> (ftHandle_);
            IppStatus status = ippsDFTFwd_RToCCS_64f(pSrc, pDst, pSpec, pBuf);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error applying DFT");
                return -1;
            }
        }
    }
    return 0;
}

int DFTR2C::getTransformLength(void) const
{
    if (!isInitialized_)
    {
        RTSEIS_ERRMSG("%s", "Class is not intiialized");
        return -1;
    }
    return lenft_;
}
 
int DFTR2C::getMaximumInputSignalLength(void) const
{
    if (!isInitialized_)
    {
        RTSEIS_ERRMSG("%s", "Class is not intiialized");
        return -1;
    }
    return length_;
}

bool DFTR2C::isInitialized(void) const
{
    return isInitialized_;
}
