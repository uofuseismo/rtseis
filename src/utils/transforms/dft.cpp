#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/transforms/transforms.hpp"
#include "rtseis/log.h"
#include <ipps.h>

using namespace RTSeis::Utilities::Transforms;

class DFTR2C::DFTImpl
{
    public:
        /// Default constructor
        DFTImpl(void)
        {
            return; 
        }
        /// Copy constructor
        DFTImpl(const DFTImpl &dft)
        {
            *this = dft;
            return;
        }
        /// Default destructor
        ~DFTImpl(void)
        {
            clear();
            return;
        }
        /// (Deep) copy operator
        DFTImpl& operator=(const DFTImpl &dftr2c)
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
            if (specSize_ > 0)
            {
                Ipp8u *pIn = nullptr;
                Ipp8u *pOut = nullptr;
                if (ldoFFT_)
                {
                    if (precision_ == RTSeis::Precision::DOUBLE)
                    {
                        pIn  = reinterpret_cast<Ipp8u *> (dftr2c.pFFTSpec64_);
                        pOut = reinterpret_cast<Ipp8u *> (pFFTSpec64_);
                        ippsCopy_8u(pIn, pOut, specSize_);
                    }
                    else
                    {
                        pIn  = reinterpret_cast<Ipp8u *> (dftr2c.pFFTSpec32_);
                        pOut = reinterpret_cast<Ipp8u *> (pFFTSpec32_);
                        ippsCopy_8u(pIn, pOut, specSize_);
                    }
                }
                else
                {
                    if (precision_ == RTSeis::Precision::DOUBLE)
                    {
                        pIn  = reinterpret_cast<Ipp8u *> (dftr2c.pDFTSpec64_);
                        pOut = reinterpret_cast<Ipp8u *> (pDFTSpec64_);
                        ippsCopy_8u(pIn, pOut, specSize_);
                    }
                    else
                    {
                        pIn  = reinterpret_cast<Ipp8u *> (dftr2c.pDFTSpec32_);
                        pOut = reinterpret_cast<Ipp8u *> (pDFTSpec32_);
                        ippsCopy_8u(pIn, pOut, specSize_);
                    }
                }
            }
            if (bufferSize_ > 0)
            {
                ippsCopy_8u(dftr2c.pBuf_, pBuf_, bufferSize_);
            }
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                if (nwork_ > 0)
                {
                    ippsCopy_64f(dftr2c.work64f_, work64f_, nwork_);
                }
            }
            else
            {
                if (nwork_ > 0)
                {
                    ippsCopy_32f(dftr2c.work32f_, work32f_, nwork_);
                }
            }   
            return *this;
        }
        /// Releases memory on the module
        void clear(void)
        {
            if (pFFTSpec64_ != nullptr){ippsFree(pFFTSpec64_);}
            if (pDFTSpec64_ != nullptr){ippsFree(pDFTSpec64_);}
            if (pFFTSpec32_ != nullptr){ippsFree(pFFTSpec32_);}
            if (pDFTSpec32_ != nullptr){ippsFree(pDFTSpec32_);}
            if (pBuf_ != nullptr){ippsFree(pBuf_);}
            if (work64f_ != nullptr){ippsFree(work64f_);}
            if (work32f_ != nullptr){ippsFree(work32f_);}
            pFFTSpec64_ = nullptr;
            pDFTSpec64_ = nullptr;
            pFFTSpec32_ = nullptr;
            pDFTSpec32_ = nullptr;
            pBuf_ = nullptr;
            work64f_ = nullptr;
            work32f_ = nullptr;
            length_ = 0;
            lenft_ = 0;
            nwork_ = 0;
            bufferSize_ = 0;
            specSize_ = 0;
            order_ = 0;
            precision_ = RTSeis::Precision::DOUBLE;
            ldoFFT_ = false;
            linit_ = false;
            return;
        }
        /// Initializes the DFT
        int initialize(const int length,
                       const bool ldoFFT,
                       const RTSeis::Precision precision)
        {
            clear();
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
            lenft_ = length_/2 + 1;
            nwork_ = std::max(length_, 2*lenft_);
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
                    pFFTSpec64_ = NULL;
                    Ipp8u *pSpecBuffer = ippsMalloc_8u(sizeInit);
                    Ipp8u *pSpec = ippsMalloc_8u(specSize_);
                    status = ippsFFTInit_R_64f(&pFFTSpec64_,
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
                    pBuf_ = ippsMalloc_8u(bufferSize_);
                }
                // DFT
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
                    pDFTSpec64_ = reinterpret_cast<IppsDFTSpec_R_64f *>
                                                  (ippsMalloc_8u(specSize_));
                    status = ippsDFTInit_R_64f(length_,
                                               IPP_FFT_DIV_INV_BY_N,
                                               ippAlgHintNone,
                                               pDFTSpec64_,
                                               pSpecBuffer);
                    if (pSpecBuffer){ippsFree(pSpecBuffer);}
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Failed to initialize DFT");
                        clear();
                        return -1;
                    }
                    pBuf_ = ippsMalloc_8u(bufferSize_);
                } // End DFT/FFT check
                work64f_ = ippsMalloc_64f(nwork_);
                ippsZero_64f(work64f_, nwork_);
            }
            // Initialize float precision
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
                    pFFTSpec32_ = NULL;
                    Ipp8u *pSpecBuffer = ippsMalloc_8u(sizeInit);
                    Ipp8u *pSpec = ippsMalloc_8u(specSize_);
                    status = ippsFFTInit_R_32f(&pFFTSpec32_,
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
                    pBuf_ = ippsMalloc_8u(bufferSize_);
                }
                // DFT
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
                    pDFTSpec32_ = reinterpret_cast<IppsDFTSpec_R_32f *>
                                                  (ippsMalloc_8u(specSize_));
                    status = ippsDFTInit_R_32f(length_,
                                               IPP_FFT_DIV_INV_BY_N,
                                               ippAlgHintNone,
                                               pDFTSpec32_,
                                               pSpecBuffer);
                    if (pSpecBuffer){ippsFree(pSpecBuffer);}
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Failed to initialize DFT");
                        clear();
                        return -1;
                    }
                    pBuf_ = ippsMalloc_8u(bufferSize_);
                } // End DFT/FFT check
                work32f_ = ippsMalloc_32f(nwork_);
                ippsZero_32f(work32f_, nwork_);
            } // End check on precision
            precision_ = precision;
            linit_ = true;
            return 0;
        }
        /// Length of the output transform
        int getTransformLength(void) const
        {
            return lenft_;
        }
        /// Length of inverse transform
        int getInverseTransformLength(void) const
        {
            return length_;
        }
        /// Maximum length of input signal
        int getMaximumInputSignalLength(void) const
        {
            return length_;
        }
        /// Determines if class is initialized
        bool isInitialized(void) const
        {
            return linit_;
        }
        /// Forward transform
        int forwardTransform(const int n, const double x[],
                             std::complex<double> y[])
        {
            // In this case the entire signal would be zero-padded
            if (n <= 0)
            {
                Ipp64f *pDst = reinterpret_cast<Ipp64f *> (y);
                ippsZero_64f(pDst, 2*lenft_);
                return 0;
            }
            // Get handle to output
            IppStatus status;
            Ipp64f *pDst = reinterpret_cast<Ipp64f *> (y);
            // No need to zero-pad
            if (n == length_)
            {
                // FFT
                if (ldoFFT_)
                {
                    status = ippsFFTFwd_RToCCS_64f(x, pDst,
                                                   pFFTSpec64_, pBuf_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Error applying FFT");
                        return -1;
                    }
                }
                // DFT
                else
                {
                    status = ippsDFTFwd_RToCCS_64f(x, pDst,
                                                   pDFTSpec64_, pBuf_);
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
                ippsCopy_64f(x, work64f_, n);
                ippsZero_64f(&work64f_[n], length_-n);
                // FFT
                if (ldoFFT_)
                {
                    status = ippsFFTFwd_RToCCS_64f(work64f_, pDst,
                                                   pFFTSpec64_, pBuf_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Error applying FFT");
                        return -1;
                    }
                }
                // DFT
                else
                {
                    status = ippsDFTFwd_RToCCS_64f(work64f_, pDst,
                                                   pDFTSpec64_, pBuf_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Error applying DFT");
                        return -1;
                    }
                }
           }
           return 0;
        }
        /// Inverse transform
        int inverseTransform(const int lenft,
                             const std::complex<double> x[],
                             double y[])
        {
            // In this case the entire signal would be zero-padded
            if (lenft <= 0)
            {
                Ipp64f *pDst = static_cast<Ipp64f *> (y);
                ippsZero_64f(pDst, length_);
                return 0;
            }
            // Get handle to output
            IppStatus status;
            //Ipp64f *pDst = static_cast<Ipp64f *> (y);
            // No need to zero-pad
            if (lenft == lenft_)
            {
                // Get handle to data
                const Ipp64f *pSrc = reinterpret_cast<const Ipp64f *> (x);
                // FFT
                if (ldoFFT_)
                {
                    status = ippsFFTInv_CCSToR_64f(pSrc, y, pFFTSpec64_, pBuf_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Error applying inverse FFT");
                        return -1; 
                    }
                }
                // DFT
                else
                {
                    status = ippsDFTInv_CCSToR_64f(pSrc, y, pDFTSpec64_, pBuf_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Error applying inverse DFT");
                        return -1;
                    }
                }
            }
            else
            {
                // Copy data and zero pad
                const Ipp64f *xTemp = reinterpret_cast<const Ipp64f *> (x);
                ippsCopy_64f(xTemp, work64f_, 2*lenft);
                ippsZero_64f(&work64f_[2*lenft], 2*(lenft_-lenft));
                // FFT
                if (ldoFFT_)
                {
                    status = ippsFFTInv_CCSToR_64f(work64f_, y,
                                                   pFFTSpec64_, pBuf_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Error applying inverse FFT");
                        return -1;
                    }
                }
                // DFT
                else
                {
                    status = ippsDFTInv_CCSToR_64f(work64f_, y,
                                                   pDFTSpec64_, pBuf_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Error applying inverse DFT");
                        return -1;
                    }
                }
            }
            return 0;
        }
    private:
        /// State structure for double FFT
        IppsFFTSpec_R_64f *pFFTSpec64_ = nullptr;
        /// State structure for double DFT
        IppsDFTSpec_R_64f *pDFTSpec64_ = nullptr;
        /// Workspace for input signals.  This has dimension [nwork_].
        Ipp64f *work64f_ = nullptr;
        /// State structure for double FFT
        IppsFFTSpec_R_32f *pFFTSpec32_ = nullptr;
        /// State structure for float FFT
        IppsDFTSpec_R_32f *pDFTSpec32_ = nullptr; 
        /// Workspace for input signals.  This has dimension [nwork_].
        Ipp32f *work32f_ = nullptr;
        /// Workspace for DFT or FFT
        Ipp8u *pBuf_ = nullptr;
        /// The maximum length of the input signal.
        int length_ = 0;
        /// The length of the Fourier transformed data.
        int lenft_ = 0;  
        /// Workspace for temporary signals.  This equals
        /// 2*std::max(length_, 2*lenft_).
        int nwork_ = 0;
        /// The length of the DFT/FFT buffer.
        int bufferSize_ = 0;
        /// The length of the DFT/FFT pointer.
        int specSize_ = 0;
        /// Specified length of FFT is 2**order.
        int order_ = 0;
        /// Precision of module.
        RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
        /// Flag indicating that the FFT is to be computed.
        bool ldoFFT_ = false;
        /// Flag indicating the module is initialized. 
        bool linit_ = false;
};

DFTR2C::DFTR2C(void) :
    pDFT_(new DFTImpl())
{
    return;
}

DFTR2C::~DFTR2C(void)
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
    if (pDFT_){pDFT_->clear();}
    clear();
    pDFT_ = std::unique_ptr<DFTImpl> (new DFTImpl(*dftr2c.pDFT_));
    return *this;
}

void DFTR2C::clear(void)
{
    pDFT_->clear();
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
    int ierr = pDFT_->initialize(length, ldoFFT, precision);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed ot initialize DFT");
        return -1;
    }
    return 0;
}

int DFTR2C::inverseTransform(const int lenft,
                             const std::complex<double> x[],
                             const int maxy, double y[])
{
    if (!pDFT_->isInitialized())
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
    int ftLen = getTransformLength();
    int length = getInverseTransformLength();
    if (lenft > ftLen || maxy < length)
    {
        if (maxy < length)
        {
            RTSEIS_ERRMSG("maxy = %d must be at least %d", maxy, length);
        }
        if (lenft > ftLen)
        {
            RTSEIS_ERRMSG("lenft = %d cannot exceed %d", lenft, ftLen);
        }
        return -1;
    }
    int ierr = pDFT_->inverseTransform(lenft, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute inverse transform");
        return -1;
    }
    return 0;
} 

int DFTR2C::forwardTransform(const int n, const double x[],
                             const int maxy, std::complex<double> y[])
{
    if (!pDFT_->isInitialized())
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
    int lenft = getTransformLength();
    int length = getInverseTransformLength();
    if (maxy < lenft || n > length)
    {
        if (maxy < lenft)
        {
            RTSEIS_ERRMSG("maxy = %d must be at least %d", maxy, lenft);
        }
        if (n > length)
        {
            RTSEIS_ERRMSG("n = %d cannot exceed %d", n, length);
        }
        return -1;
    }
    int ierr = pDFT_->forwardTransform(n, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply forward transform");
        return -1;
    }
    return 0;
}

int DFTR2C::getTransformLength(void) const
{
    if (!pDFT_->isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class is not intiialized");
        return -1;
    }
    return pDFT_->getTransformLength();
}

int DFTR2C::getInverseTransformLength(void) const
{
    if (!pDFT_->isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class is not intiialized");
        return -1;
    }
    return pDFT_->getInverseTransformLength();
}
 
int DFTR2C::getMaximumInputSignalLength(void) const
{
    if (!pDFT_->isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class is not intiialized");
        return -1;
    }
    return pDFT_->getMaximumInputSignalLength();
}

bool DFTR2C::isInitialized(void) const
{
    return pDFT_->isInitialized();
}
