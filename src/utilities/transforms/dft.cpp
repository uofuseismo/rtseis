#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/transforms/enums.hpp"
#include "rtseis/utilities/transforms/dft.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::Transforms;

class DFT::DFTImpl
{
public:
    /// Default constructor
    DFTImpl() = default;
    /// Copy constructor
    DFTImpl(const DFTImpl &dft)
    {
        *this = dft;
    }
    /// Default destructor
    ~DFTImpl()
    {
        clear();
    }
    /// (Deep) copy operator
    DFTImpl& operator=(const DFTImpl &dft)
    {
        if (&dft == this){return *this;}
        clear();
        if (!dft.isInitialized()){return *this;}
        int ierr = initialize(dft.length_, dft.ldoFFT_,
                              dft.precision_);
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
                    pIn  = reinterpret_cast<Ipp8u *> (dft.pFFTSpec64_);
                    pOut = reinterpret_cast<Ipp8u *> (pFFTSpec64_);
                    ippsCopy_8u(pIn, pOut, specSize_);
                }
                else
                {
                    pIn  = reinterpret_cast<Ipp8u *> (dft.pFFTSpec32_);
                    pOut = reinterpret_cast<Ipp8u *> (pFFTSpec32_);
                    ippsCopy_8u(pIn, pOut, specSize_);
                }
            }
            else
            {
                if (precision_ == RTSeis::Precision::DOUBLE)
                {
                    pIn  = reinterpret_cast<Ipp8u *> (dft.pDFTSpec64_);
                    pOut = reinterpret_cast<Ipp8u *> (pDFTSpec64_);
                    ippsCopy_8u(pIn, pOut, specSize_);
                }
                else
                {
                    pIn  = reinterpret_cast<Ipp8u *> (dft.pDFTSpec32_);
                    pOut = reinterpret_cast<Ipp8u *> (pDFTSpec32_);
                    ippsCopy_8u(pIn, pOut, specSize_);
                }
            }
        }
        if (bufferSize_ > 0)
        {
            ippsCopy_8u(dft.pBuf_, pBuf_, bufferSize_);
        }
        if (precision_ == RTSeis::Precision::DOUBLE)
        {
            if (nwork_ > 0)
            {
                ippsCopy_64fc(dft.work64fc_, work64fc_, nwork_);
            }
        }
        else
        {
            if (nwork_ > 0)
            {
                ippsCopy_32fc(dft.work32fc_, work32fc_, nwork_);
            }
        }   
        return *this;
    }
    /// Releases memory on the module
    void clear()
    {
        if (pFFTSpec64_ != nullptr){ippsFree(pFFTSpec64_);}
        if (pDFTSpec64_ != nullptr){ippsFree(pDFTSpec64_);}
        if (pFFTSpec32_ != nullptr){ippsFree(pFFTSpec32_);}
        if (pDFTSpec32_ != nullptr){ippsFree(pDFTSpec32_);}
        if (pBuf_ != nullptr){ippsFree(pBuf_);}
        if (work64fc_ != nullptr){ippsFree(work64fc_);}
        if (work32fc_ != nullptr){ippsFree(work32fc_);}
        pFFTSpec64_ = nullptr;
        pDFTSpec64_ = nullptr;
        pFFTSpec32_ = nullptr;
        pDFTSpec32_ = nullptr;
        pBuf_ = nullptr;
        work64fc_ = nullptr;
        work32fc_ = nullptr;
        length_ = 0;
        lenft_ = 0;
        nwork_ = 0;
        bufferSize_ = 0;
        specSize_ = 0;
        order_ = 0;
        precision_ = RTSeis::Precision::DOUBLE;
        ldoFFT_ = false;
        linit_ = false;
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
        lenft_ = length_;
        nwork_ = 2*lenft_;
        // Initialize the appropriate transform
        IppStatus status;
        int sizeInit;
        if (precision == RTSeis::Precision::DOUBLE)
        {
            if (ldoFFT_)
            {
                status = ippsFFTGetSize_C_64fc(order_,
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
                pFFTSpec64_ = nullptr;
                Ipp8u *pSpecBuffer = ippsMalloc_8u(sizeInit);
                Ipp8u *pSpec = ippsMalloc_8u(specSize_);
                status = ippsFFTInit_C_64fc(&pFFTSpec64_,
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
                status = ippsDFTGetSize_C_64fc(length_,
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
                pDFTSpec64_ = reinterpret_cast<IppsDFTSpec_C_64fc *>
                                              (ippsMalloc_8u(specSize_));
                status = ippsDFTInit_C_64fc(length_,
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
            work64fc_ = ippsMalloc_64fc(nwork_);
            ippsZero_64fc(work64fc_, nwork_);
        }
        // Initialize float precision
        else
        {
            if (ldoFFT_)
            {
                status = ippsFFTGetSize_C_32fc(order_,
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
                pFFTSpec32_ = nullptr;
                Ipp8u *pSpecBuffer = ippsMalloc_8u(sizeInit);
                Ipp8u *pSpec = ippsMalloc_8u(specSize_);
                status = ippsFFTInit_C_32fc(&pFFTSpec32_,
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
                status = ippsDFTGetSize_C_32fc(length_,
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
                pDFTSpec32_ = reinterpret_cast<IppsDFTSpec_C_32fc *>
                                              (ippsMalloc_8u(specSize_));
                status = ippsDFTInit_C_32fc(length_,
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
            work32fc_ = ippsMalloc_32fc(nwork_);
            ippsZero_32fc(work32fc_, nwork_);
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
    int forwardTransform(const int n, const std::complex<double> x[],
                         std::complex<double> y[])
    {
        // In this case the entire signal would be zero-padded
        if (n <= 0)
        {
            auto *pDst = reinterpret_cast<Ipp64fc *> (y);
            ippsZero_64fc(pDst, lenft_);
            return 0;
        }
        // Get handle to output
        IppStatus status;
        auto *pSrc = reinterpret_cast<const Ipp64fc *> (x);
        auto *pDst = reinterpret_cast<Ipp64fc *> (y);
        // No need to zero-pad
        if (n == length_)
        {
            // FFT
            if (ldoFFT_)
            {
                status = ippsFFTFwd_CToC_64fc(pSrc, pDst,
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
                status = ippsDFTFwd_CToC_64fc(pSrc, pDst,
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
            ippsCopy_64fc(pSrc, work64fc_, n);
            ippsZero_64fc(&work64fc_[n], length_-n);
            // FFT
            if (ldoFFT_)
            {
                status = ippsFFTFwd_CToC_64fc(work64fc_, pDst,
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
                status = ippsDFTFwd_CToC_64fc(work64fc_, pDst,
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
    /// Forward transform
    int forwardTransform(const int n, const std::complex<float> x[],
                         std::complex<float> y[])
    {
        // In this case the entire signal would be zero-padded
        if (n <= 0)
        {
            auto *pDst = reinterpret_cast<Ipp32fc *> (y);
            ippsZero_32fc(pDst, lenft_);
            return 0;
        }
        // Get handle to output
        IppStatus status;
        auto *pSrc = reinterpret_cast<const Ipp32fc *> (x);
        auto *pDst = reinterpret_cast<Ipp32fc *> (y);
        // No need to zero-pad
        if (n == length_)
        {
            // FFT
            if (ldoFFT_)
            {
                status = ippsFFTFwd_CToC_32fc(pSrc, pDst,
                                              pFFTSpec32_, pBuf_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Error applying FFT");
                    return -1;
                }
            }
            // DFT
            else
            {
                status = ippsDFTFwd_CToC_32fc(pSrc, pDst,
                                              pDFTSpec32_, pBuf_);
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
            ippsCopy_32fc(pSrc, work32fc_, n);
            ippsZero_32fc(&work32fc_[n], length_-n);
            // FFT
            if (ldoFFT_)
            {
                status = ippsFFTFwd_CToC_32fc(work32fc_, pDst,
                                              pFFTSpec32_, pBuf_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Error applying FFT");
                    return -1;
                }
            }
            // DFT
            else
            {
                status = ippsDFTFwd_CToC_32fc(work32fc_, pDst,
                                              pDFTSpec32_, pBuf_);
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
                         std::complex<double> y[])
    {
        // In this case the entire signal would be zero-padded
        if (lenft <= 0)
        {
            auto *pDst = reinterpret_cast<Ipp64fc *> (y);
            ippsZero_64fc(pDst, length_);
            return 0;
        }
        IppStatus status;
        auto *pSrc = reinterpret_cast<const Ipp64fc *> (x);
        auto *pDst = reinterpret_cast<Ipp64fc *> (y);
        // No need to zero-pad
        if (lenft == lenft_)
        {
            // FFT
            if (ldoFFT_)
            {
                status = ippsFFTInv_CToC_64fc(pSrc, pDst, pFFTSpec64_, pBuf_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Error applying inverse FFT");
                    return -1; 
                }
            }
            // DFT
            else
            {
                status = ippsDFTInv_CToC_64fc(pSrc, pDst, pDFTSpec64_, pBuf_);
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
            ippsCopy_64fc(pSrc, work64fc_, lenft);
            ippsZero_64fc(&work64fc_[lenft], lenft_-lenft);
            // FFT
            if (ldoFFT_)
            {
                status = ippsFFTInv_CToC_64fc(work64fc_, pDst,
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
                status = ippsDFTInv_CToC_64fc(work64fc_, pDst,
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
    /// Inverse transform
    int inverseTransform(const int lenft,
                         const std::complex<float> x[],
                         std::complex<float> y[])
    {
        // In this case the entire signal would be zero-padded
        if (lenft <= 0)
        {
            auto *pDst = reinterpret_cast<Ipp32fc *> (y);
            ippsZero_32fc(pDst, length_);
            return 0;
        }
        IppStatus status;
        auto *pSrc = reinterpret_cast<const Ipp32fc *> (x);
        auto *pDst = reinterpret_cast<Ipp32fc *> (y);
        // No need to zero-pad
        if (lenft == lenft_)
        {
            // FFT
            if (ldoFFT_)
            {
                status = ippsFFTInv_CToC_32fc(pSrc, pDst, pFFTSpec32_, pBuf_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Error applying inverse FFT");
                    return -1; 
                }
            }
            // DFT
            else
            {
                status = ippsDFTInv_CToC_32fc(pSrc, pDst, pDFTSpec32_, pBuf_);
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
            ippsCopy_32fc(pSrc, work32fc_, lenft);
            ippsZero_32fc(&work32fc_[lenft], lenft_-lenft);
            // FFT
            if (ldoFFT_)
            {
                status = ippsFFTInv_CToC_32fc(work32fc_, pDst,
                                              pFFTSpec32_, pBuf_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Error applying inverse FFT");
                    return -1;
                }
            }
            // DFT
            else
            {
                status = ippsDFTInv_CToC_32fc(work32fc_, pDst,
                                              pDFTSpec32_, pBuf_);
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
    IppsFFTSpec_C_64fc *pFFTSpec64_ = nullptr;
    /// State structure for double DFT
    IppsDFTSpec_C_64fc *pDFTSpec64_ = nullptr;
    /// Workspace for input signals.  This has dimension [nwork_].
    Ipp64fc *work64fc_ = nullptr;
    /// State structure for double FFT
    IppsFFTSpec_C_32fc *pFFTSpec32_ = nullptr;
    /// State structure for float FFT
    IppsDFTSpec_C_32fc *pDFTSpec32_ = nullptr; 
    /// Workspace for input signals.  This has dimension [nwork_].
    Ipp32fc *work32fc_ = nullptr;
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

DFT::DFT() :
    pImpl(std::make_unique<DFTImpl>())
{
}

DFT::~DFT() = default;

DFT::DFT(const DFT &dft)
{
    *this = dft;
}

/*
DFT::DFT(DFT &&dft) noexcept
{
    *this = std::Move(dft);
}
*/

DFT& DFT::operator=(const DFT &dft)
{
    if (&dft == this){return *this;}
    if (pImpl){pImpl->clear();}
    clear();
    pImpl = std::make_unique<DFTImpl> (*dft.pImpl);
    return *this;
}

/*
DFT& DFT::operator=(DFT &&dft) noexcept
{
    if (&dft == this){return *this;}
    pImpl = std::move(hilbert.pImpl);
    return *this;
}
*/

void DFT::clear() noexcept
{
    pImpl->clear();
}


void DFT::initialize(
    const int length,
    const FourierTransformImplementation implementation,
    const RTSeis::Precision precision)
{
    clear();
    // Check the inputs
    if (length < 2)
    {
        RTSEIS_THROW_IA("Length=%d must be at least 2", length);
    }
    bool ldoFFT = false;
    if (implementation == FourierTransformImplementation::FFT)
    {
        ldoFFT = true;
    }
    int ierr = pImpl->initialize(length, ldoFFT, precision);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed ot initialize DFT");
    }
    return;
}

void DFT::inverseTransform(const int lenft,
                           const std::complex<double> x[],
                           const int maxy, std::complex<double> y[])
{
    if (!pImpl->isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not intiialized");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    int ftLen = getTransformLength();
    int length = getInverseTransformLength();
    if (lenft > ftLen || maxy < length)
    {
        if (maxy < length)
        {
            RTSEIS_THROW_IA("maxy = %d must be at least %d", maxy, length);
        }
        RTSEIS_THROW_IA("lenft = %d cannot exceed %d", lenft, ftLen);
    }
    int ierr = pImpl->inverseTransform(lenft, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute inverse transform");
    }
    return;
} 

void DFT::forwardTransform(const int n, const std::complex<double> x[],
                           const int maxy, std::complex<double> y[])
{
    if (!pImpl->isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not intiialized");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    int lenft = getTransformLength();
    int length = getInverseTransformLength();
    if (maxy < lenft || n > length)
    {
        if (maxy < lenft)
        {
            RTSEIS_THROW_IA("maxy = %d must be at least %d", maxy, lenft);
        }
        RTSEIS_THROW_IA("n = %d cannot exceed %d", n, length);
    }
    int ierr = pImpl->forwardTransform(n, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply forward transform");
        return;
    }
    return;
}

int DFT::getTransformLength() const
{
    if (!pImpl->isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not initialized");
    }
    return pImpl->getTransformLength();
}

int DFT::getInverseTransformLength() const
{
    if (!pImpl->isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not intiialized");
    }
    return pImpl->getInverseTransformLength();
}
 
int DFT::getMaximumInputSignalLength() const
{
    if (!pImpl->isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not intiialized");
    }
    return pImpl->getMaximumInputSignalLength();
}

bool DFT::isInitialized() const noexcept
{
    return pImpl->isInitialized();
}
