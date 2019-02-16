#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ipps.h>
#include <ippversion.h>
#include <ippcore.h>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/filters.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::Filters;

class IIRFilter::IIRFilterImpl
{
    public:
        /// Default constructor
        IIRFilterImpl(void)
        {
            return;
        }
        /// Copy constructor
        IIRFilterImpl(const IIRFilterImpl &iir)
        {
            *this = iir;
            return;
        }
        /// Default destructor
        ~IIRFilterImpl(void) 
        {
            clear();
            return;
        }
        /// Copy operator
        IIRFilterImpl& operator=(const IIRFilterImpl &iir)
        {
            if (&iir == this){return *this;}
            clear();
            if (!iir.linit_){return *this;}
            int ierr = initialize(iir.nbRef_, iir.bRef_,
                                  iir.naRef_, iir.aRef_,
                                  iir.mode_, iir.precision_,
                                  iir.implementation_);
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Failed to initialize filter");
                clear();
                return *this;
            }
            if (bufferSize_ > 0)
            {
                ippsCopy_8u(iir.pBuf_, pBuf_, bufferSize_);
            }
            ippsCopy_64f(iir.zi_, zi_, nbDly_);
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(iir.pBufIPP64f_, pBufIPP64f_, bufIPPLen_);
                ippsCopy_64f(iir.pDlySrc64f_, pDlySrc64f_, nbDly_);
                ippsCopy_64f(iir.pDlyDst64f_, pDlyDst64f_, nbDly_);
            }
            else
            {
                ippsCopy_32f(iir.pBufIPP32f_, pBufIPP32f_, bufIPPLen_);
                ippsCopy_32f(iir.pDlySrc32f_, pDlySrc32f_, nbDly_);
                ippsCopy_32f(iir.pDlyDst32f_, pDlyDst32f_, nbDly_);
            }
            return *this;
        }
        /// Releases memory on the module
        void clear(void)
        {
            if (pTaps64f_ != nullptr){ippsFree(pTaps64f_);}
            if (pBufIPP64f_ != nullptr){ippsFree(pBufIPP64f_);}
            if (pDlySrc64f_ != nullptr){ippsFree(pDlySrc64f_);}
            if (pDlyDst64f_ != nullptr){ippsFree(pDlyDst64f_);}
            if (pTaps32f_ != nullptr){ippsFree(pTaps32f_);}
            if (pBufIPP32f_ != nullptr){ippsFree(pBufIPP32f_);}
            if (pDlySrc32f_ != nullptr){ippsFree(pDlySrc32f_);}
            if (pDlyDst32f_ != nullptr){ippsFree(pDlyDst32f_);}
            if (pBuf_ != nullptr){ippsFree(pBuf_);}
            if (bRef_ != nullptr){ippsFree(bRef_);}
            if (aRef_ != nullptr){ippsFree(aRef_);}
            if (bNorm64f_ != nullptr){ippsFree(bNorm64f_);}
            if (aNorm64f_ != nullptr){ippsFree(aNorm64f_);}
            if (bNorm32f_ != nullptr){ippsFree(bNorm32f_);}
            if (aNorm32f_ != nullptr){ippsFree(aNorm32f_);}
            if (zi_ != nullptr){ippsFree(zi_);}
            pIIRState64f_ = nullptr;
            pTaps64f_ = nullptr;
            pBufIPP64f_ = nullptr;
            pDlySrc64f_ = nullptr;
            pDlyDst64f_ = nullptr;
            pIIRState32f_ = nullptr;;
            pTaps32f_ = nullptr;
            pBufIPP32f_ = nullptr;
            pDlySrc32f_ = nullptr;
            pDlyDst32f_ = nullptr;
            pBuf_ = nullptr;
            bRef_ = nullptr;
            aRef_ = nullptr;
            bNorm64f_ = nullptr;
            aNorm64f_ = nullptr;
            zi_ = nullptr;
            nbDly_ = 0;
            bufIPPLen_ = 0;
            order_ = 0;
            nbRef_ = 0;
            naRef_ = 0;
            bufferSize_ = 0;
            implementation_ = IIRFilter::Implementation::DF2_FAST;
            mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
            precision_ = RTSeis::Precision::DOUBLE;
            linit_ = false;
            return;
        }
        //====================================================================//
        /// Initializes the filter
        int initialize(const int nb, const double b[],
                       const int na, const double a[],
                       const RTSeis::ProcessingMode mode,
                       const RTSeis::Precision precision,
                       const IIRFilter::Implementation implementation)
        {
            clear();
            // May have to force the implementation to be slow
            IppStatus status;
            IIRFilter::Implementation impUse = implementation;
            // IPP changed implementation in 2018.  I'm not supporting it.
            // High order IIR filters aren't a good idea.
            #if IPP_VERSION_MAJOR < 2019
            impUse = IIRFilter::Implementation::DF2_SLOW;
            #endif
            if (impUse == IIRFilter::Implementation::DF2_FAST &&
                std::max(nb, na) - 1 > 8)
            {
                Ipp64u featureMask, enabledMask;
                status = ippGetCpuFeatures(&featureMask, NULL);
                if (status == ippStsNoErr)
                {
                    enabledMask = ippGetEnabledCpuFeatures();
                    if ((featureMask & ippCPUID_AVX2) &&
                        (enabledMask & ippCPUID_AVX2))
                    {
                        impUse = IIRFilter::Implementation::DF2_SLOW;
                    }
                    if ((featureMask & ippCPUID_AVX512F) &&
                        (enabledMask & ippCPUID_AVX512F))
                    {
                        impUse = IIRFilter::Implementation::DF2_SLOW;
                    }
                }
                else
                {
                    RTSEIS_ERRMSG("%s", "Failed to get processor info");
                    return -1;
                }
            }
            if (impUse != implementation)
            {
                RTSEIS_WARNMSG("%s", "Overriding implementation to DF2_SLOW");
            }
            // Set sizes
            double a0 = a[0];
            nbRef_ = nb;
            naRef_ = na;
            order_ = std::max(nbRef_, naRef_) - 1;
            bufIPPLen_ = 2*std::max(order_+1, 1024);
            nbDly_ = std::max(8, order_ + 1);
            // Copy and normalize the filter coefficients
            bRef_ = ippsMalloc_64f(nbRef_);
            ippsCopy_64f(b, bRef_, nbRef_);
            aRef_ = ippsMalloc_64f(naRef_);
            ippsCopy_64f(a, aRef_, naRef_);
            zi_ = ippsMalloc_64f(nbDly_);
            ippsZero_64f(zi_, nbDly_);
            // Alllocate space
            if (precision == RTSeis::Precision::DOUBLE)
            {
                if (impUse == Implementation::DF2_FAST)
                {
                    status = ippsIIRGetStateSize_64f(order_, &bufferSize_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Failed to get state size");
                        clear();
                        return -1;
                    }
                    // Set the workspace
                    pBufIPP64f_ = ippsMalloc_64f(bufIPPLen_);
                    ippsZero_64f(pBufIPP64f_, bufIPPLen_);
                    pDlySrc64f_ = ippsMalloc_64f(nbDly_);
                    ippsZero_64f(pDlySrc64f_, nbDly_);
                    pDlyDst64f_ = ippsMalloc_64f(nbDly_);
                    ippsZero_64f(pDlyDst64f_, nbDly_);
                    pBuf_ = ippsMalloc_8u(bufferSize_);
                    // Set the (normalized) filter taps
                    pTaps64f_ = ippsMalloc_64f(2*(order_ + 1));
                    ippsZero_64f(pTaps64f_, 2*(order_ + 1));
                    ippsDivC_64f(bRef_, a0, &pTaps64f_[0],        nbRef_);
                    ippsDivC_64f(aRef_, a0, &pTaps64f_[order_+1], naRef_);
                    // Initialize the filter
                    status = ippsIIRInit_64f(&pIIRState64f_, pTaps64f_, order_,
                                             pDlySrc64f_, pBuf_);

                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
                        clear();
                        return -1;
                    }
                    // Set the delay line
                    status = ippsIIRSetDlyLine_64f(pIIRState64f_, pBufIPP64f_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Failed to set delay line");
                        clear();
                        return -1;
                    }
                }
                else
                {
                    bNorm64f_ = ippsMalloc_64f(order_+1);
                    ippsZero_64f(bNorm64f_, order_+1);
                    ippsDivC_64f(bRef_,a0,  bNorm64f_, nbRef_);
                    aNorm64f_ = ippsMalloc_64f(order_+1);
                    ippsZero_64f(aNorm64f_, order_+1);
                    ippsDivC_64f(aRef_, a0, aNorm64f_, naRef_);
                    // Initial/final conditions
                    pDlySrc64f_ = ippsMalloc_64f(nbDly_);
                    ippsZero_64f(pDlySrc64f_, nbDly_);
                    pDlyDst64f_ = ippsMalloc_64f(nbDly_);
                    ippsZero_64f(pDlyDst64f_, nbDly_);
                }
            }
            else
            {
                if (impUse == Implementation::DF2_FAST)
                {
                    status = ippsIIRGetStateSize_32f(order_, &bufferSize_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Failed to get state size");
                        clear();
                        return -1;
                    }
                    // Set the workspace
                    pBufIPP32f_ = ippsMalloc_32f(bufIPPLen_);
                    ippsZero_32f(pBufIPP32f_, bufIPPLen_);
                    pDlySrc32f_ = ippsMalloc_32f(nbDly_);
                    ippsZero_32f(pDlySrc32f_, nbDly_);
                    pDlyDst32f_ = ippsMalloc_32f(nbDly_); 
                    ippsZero_32f(pDlyDst32f_, nbDly_);
                    pBuf_ = ippsMalloc_8u(bufferSize_);
                    // Set the (normalized) filter taps
                    pTaps32f_ = ippsMalloc_32f(2*(order_ + 1));
                    ippsZero_32f(pTaps32f_, 2*(order_ + 1));
                    ippsConvert_64f32f(bRef_, &pTaps32f_[0],        nbRef_);
                    ippsConvert_64f32f(aRef_, &pTaps32f_[order_+1], naRef_);
                    float a04 = static_cast<float> (a0);
                    ippsDivC_32f_I(a04, &pTaps32f_[0],        nbRef_);
                    ippsDivC_32f_I(a04, &pTaps32f_[order_+1], naRef_);
                    // Initialize the filter
                    status = ippsIIRInit_32f(&pIIRState32f_, pTaps32f_, order_,
                                             pDlySrc32f_, pBuf_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
                        clear();
                        return -1;
                    }
                    // Set the delay line
                    status = ippsIIRSetDlyLine_32f(pIIRState32f_, pBufIPP32f_);
                    if (status != ippStsNoErr)
                    {
                        RTSEIS_ERRMSG("%s", "Failed to set delay line");
                        clear();
                        return -1;
                    }
                }
                else
                {
                    float a04 = static_cast<float> (a0);
                    bNorm32f_ = ippsMalloc_32f(order_+1);
                    ippsZero_32f(bNorm32f_, order_+1);
                    ippsConvert_64f32f(bRef_, bNorm32f_, nbRef_);
                    ippsDivC_32f_I(a04, bNorm32f_, nbRef_);
                    aNorm32f_ = ippsMalloc_32f(order_+1);
                    ippsZero_32f(aNorm32f_, order_+1);
                    ippsConvert_64f32f(aRef_, aNorm32f_, naRef_);
                    ippsDivC_32f_I(a04, aNorm32f_, naRef_);
                    // Initial/final conditions
                    pDlySrc32f_ = ippsMalloc_32f(nbDly_);
                    ippsZero_32f(pDlySrc32f_, nbDly_);
                    pDlyDst32f_ = ippsMalloc_32f(nbDly_);
                    ippsZero_32f(pDlyDst32f_, nbDly_);
                }
            }
            implementation_ = impUse;
            mode_ = mode;
            precision_ = precision;
            linit_ = true;
            return 0;
        }
        /// Determines if the filter is initialized
        bool isInitialized(void) const
        {
            return linit_;
        }
        /// Gets the length of the initial conditions
        int getInitialConditionLength(void) const
        {
            return order_;
        }
        /// Sets the initial conditions
        int setInitialConditions(const int nz, const double zi[])
        {
            resetInitialConditions();
            int nzRef = getInitialConditionLength();
            if (nzRef != nz)
            {
                RTSEIS_ERRMSG("%s", "Shouldn't be here");
            }
            if (nzRef == 0){return 0;}
            ippsCopy_64f(zi, zi_, nzRef);
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(zi_, pBufIPP64f_, nzRef); 
                ippsIIRSetDlyLine_64f(pIIRState64f_, pBufIPP64f_);
            }
            else
            {
                ippsConvert_64f32f(zi_, pBufIPP32f_, nzRef);
                ippsIIRSetDlyLine_32f(pIIRState32f_, pBufIPP32f_);
            }
            return 0;
        } 
        /// Resets the initial conditions
        int resetInitialConditions(void)
        {
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                ippsZero_64f(pBufIPP64f_, bufIPPLen_);
                ippsZero_64f(pDlySrc64f_, nbDly_);
                ippsZero_64f(pDlyDst64f_, nbDly_);
                if (order_ > 0)
                {
                    ippsCopy_64f(zi_, pDlySrc64f_, order_);
                }
                ippsIIRSetDlyLine_64f(pIIRState64f_, pBufIPP64f_);
            }
            else
            {
                ippsZero_32f(pBufIPP32f_, bufIPPLen_);
                ippsZero_32f(pDlySrc32f_, nbDly_);
                ippsZero_32f(pDlyDst32f_, nbDly_);
                if (order_ > 0)
                {
                    ippsConvert_64f32f(zi_, pDlySrc32f_, order_);
                }
                ippsIIRSetDlyLine_32f(pIIRState32f_, pBufIPP32f_);
            }
            return 0;
        } 
        /// Applies the filter
        int apply(const int n, const double x[], double y[])
        {
            if (n <= 0){return 0;}
            if (precision_ == RTSeis::Precision::FLOAT)
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
            if (implementation_ == Implementation::DF2_FAST)
            {
                status = ippsIIR_64f(x, y, n, pIIRState64f_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Failed to set delay line");
                    return -1;
                }
                if (mode_ == RTSeis::ProcessingMode::REAL_TIME)
                {
                    ippsIIRGetDlyLine_64f(pIIRState64f_, pBufIPP64f_);
                }
            }
            else
            {
                iirDF2Transpose(n, x, y);
            }
            return 0; 
        }
        /// Applies the filter
        int apply(const int n, const float x[], float y[])
        {
            if (n <= 0){return 0;}
            if (precision_ == RTSeis::Precision::DOUBLE)
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
            if (implementation_ == Implementation::DF2_FAST)
            {
                status = ippsIIR_32f(x, y, n, pIIRState32f_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Failed to set delay line");
                    return -1;
                }
                if (mode_ == RTSeis::ProcessingMode::REAL_TIME)
                {
                    ippsIIRGetDlyLine_32f(pIIRState32f_, pBufIPP32f_);
                }
            }
            else
            {
                iirDF2Transpose(n, x, y);
            }
            return 0;
        }
        /// A more numerically robust yet slower filter implementation
        int iirDF2Transpose(const int n, const double x[], double y[])
        {
            Ipp64f *b = bNorm64f_;
            Ipp64f *a = aNorm64f_;
            Ipp64f *vi = pDlySrc64f_;
            Ipp64f *v  = pDlyDst64f_;
            // Loop on samples
            for (int i=0; i<n; i++)
            {
                #pragma omp simd
                for (int j=order_; j>=1; j--){v[j] = vi[j-1];}
                double Xi = 0;
                double Yi = 0;
                #pragma omp simd reduction(+:Xi,Yi)
                for (int j=1; j<=order_; j++)
                {
                    Xi = Xi - a[j]*v[j];
                    Yi = Yi + b[j]*v[j];
                }
                Xi   = Xi + x[i];
                v[0] = Xi;
                y[i] = Yi + b[0]*v[0];
                #pragma omp simd
                for (int j=0; j<=order_; j++){vi[j] = v[j];}
            }
            if (mode_ == RTSeis::ProcessingMode::POST_PROCESSING)
            {
                #pragma omp simd
                for (int j=0; j<=order_; j++){vi[j] = 0;}
            }
            return 0;
        }
        /// A more numerically robust yet slower filter implementation
        int iirDF2Transpose(const int n, const float x[], float y[])
        {
            Ipp32f *b = bNorm32f_;
            Ipp32f *a = aNorm32f_;
            Ipp32f *vi = pDlySrc32f_;
            Ipp32f *v  = pDlyDst32f_;
            // Loop on samples
            for (int i=0; i<n; i++)
            {
                #pragma omp simd
                for (int j=order_; j>=1; j--){v[j] = vi[j-1];}
                float Xi = 0;
                float Yi = 0;
                #pragma omp simd reduction(+:Xi,Yi)
                for (int j=1; j<=order_; j++)
                {
                    Xi = Xi - a[j]*v[j];
                    Yi = Yi + b[j]*v[j];
                }   
                Xi   = Xi + x[i];
                v[0] = Xi; 
                y[i] = Yi + b[0]*v[0];
                #pragma omp simd
                for (int j=0; j<=order_; j++){vi[j] = v[j];}
            }   
            if (mode_ == RTSeis::ProcessingMode::POST_PROCESSING)
            {   
                #pragma omp simd
                for (int j=0; j<=order_; j++){vi[j] = 0;} 
            }
            return 0;
        }
    private:
        /// IIR filtering state
        IppsIIRState_64f *pIIRState64f_ = nullptr;
        /// The Filter taps.  This has dimension [2*(order_+1)].
        Ipp64f *pTaps64f_ = nullptr;
        /// Every 1024 points sampled we need to rollover the IIR filter
        /// in IPP.  This has dimension [bufIPPLen_].
        Ipp64f *pBufIPP64f_ = nullptr;
        /// Holds the input IIR filter delay line.  This has dimension
        /// [nbDly_].
        Ipp64f *pDlySrc64f_ = nullptr;
        /// Holds the output IIR filter delay line.  This has dimension
        /// [nbDly_].
        Ipp64f *pDlyDst64f_ = nullptr;
        /// IIR filtering state
        IppsIIRState_32f *pIIRState32f_ = nullptr;;
        /// The Filter taps.  This has dimension [2*(order_+1)].
        Ipp32f *pTaps32f_ = nullptr;
        /// Holds the input IIR filter conditions.  This has dimension
        /// [bufIPPLen_].
        Ipp32f *pBufIPP32f_ = nullptr;
        /// Holds the input IIR filter delay line.  This has dimension
        /// [nbDly_].
        Ipp32f *pDlySrc32f_ = nullptr;
        /// Holds the output IIR filter delay line.  This has dimension
        /// [nbDly_].
        Ipp32f *pDlyDst32f_ = nullptr;
        /// The workspace buffer for the filter.
        Ipp8u *pBuf_ = nullptr;
        /// The reference filter numerator coefficients.
        /// This has dimension [nbRef_].
        Ipp64f *bRef_ = nullptr;
        /// The reference filter denominator coefficients.
        /// This has dimension [naRef_].
        Ipp64f *aRef_ = nullptr;
        /// These are the normalized numerator coefficients.
        /// This has dimension [order_+1].
        Ipp64f *bNorm64f_ = nullptr;
        /// These are the normalized denominator coefficients.
        /// This has dimension [order_+1].
        Ipp64f *aNorm64f_ = nullptr;
        /// These are the normalized numerator coefficients.
        /// This has dimension [order_+1].
        Ipp32f *bNorm32f_ = nullptr;
        /// These are the normalized denominator coefficients.
        /// This has dimension [order_+1].
        Ipp32f *aNorm32f_ = nullptr;
        /// Holds a copy of the initial conditions
        Ipp64f *zi_ = nullptr;
        /// The length of the delay line.  This is length order_ + 1.
        int nbDly_ = 0;
        /// The length of the IPP buffer workspace.
        int bufIPPLen_ = 0;
        /// The filter order = max(nbRef_, naRef_) - 1. 
        int order_ = 0;
        /// Reference number of numerator coefficients
        int nbRef_ = 0; 
        /// Reference number of denominator coefficients
        int naRef_ = 0;
        /// The length of the pBuf.
        int bufferSize_ = 0;
        /// Filter implementation
        Implementation implementation_ = IIRFilter::Implementation::DF2_FAST;
        /// Processing mode
        RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        /// Precision of filter application
        RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
        /// Flag indicating the module is initialized
        bool linit_ = false;
};

//============================================================================//

IIRFilter::IIRFilter(void) :
    pIIR_(new IIRFilterImpl())
{
    return;
}

IIRFilter::IIRFilter(const IIRFilter &iir)
{
    *this = iir;
    return;
}

IIRFilter& IIRFilter::operator=(const IIRFilter &iir)
{
    if (&iir == this){return *this;}
    if (pIIR_){pIIR_->clear();}
    pIIR_ = std::unique_ptr<IIRFilterImpl> (new IIRFilterImpl(*iir.pIIR_));
    return *this;
}

IIRFilter::~IIRFilter(void)
{
    clear();
    return;
}

void IIRFilter::clear(void)
{
    pIIR_->clear();
    return;
}

int IIRFilter::initialize(const int nb, const double b[],
                          const int na, const double a[],
                          const RTSeis::ProcessingMode mode,
                          const RTSeis::Precision precision,
                          const Implementation implementation)
{
    clear();
    if (nb < 1 || b == nullptr || na < 1 || a == nullptr)
    {
        if (nb < 1){RTSEIS_ERRMSG("%s", "No numerator coefficients");}
        if (na < 1){RTSEIS_ERRMSG("%s", "No denominator coefficients");}
        if (b == nullptr){RTSEIS_ERRMSG("%s", "b is NULL");}
        if (a == nullptr){RTSEIS_ERRMSG("%s", "a is NULL");}
        return -1;
    }
    if (a[0] == 0)
    {
        RTSEIS_ERRMSG("%s", "a[0] cannot be zero");
        return -1;
    }
    int ierr = pIIR_->initialize(nb, b, na, a, mode, precision, implementation);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        clear();
        return -1;
    }
    return 0;
}

int IIRFilter::getInitialConditionLength(void) const
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    return pIIR_->getInitialConditionLength();
}

int IIRFilter::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    int nzRef = pIIR_->getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef){RTSEIS_ERRMSG("nz=%d should equal %d", nz, nzRef);}
        if (zi == nullptr){RTSEIS_ERRMSG("%s", "zi is NULL");}
        return -1;
    }
    pIIR_->setInitialConditions(nz, zi);
    return 0;
}

int IIRFilter::resetInitialConditions(void)
{
    if (!isInitialized())
    {   
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1; 
    }   
    int ierr = pIIR_->resetInitialConditions();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to reset initial conditions");
        return -1;
    }
    return 0;
}

int IIRFilter::apply(const int n, const float x[], float y[])
{
    if (n <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
    }
    int ierr = pIIR_->apply(n, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
    }
    return 0;
}

int IIRFilter::apply(const int n, const double x[], double y[])
{
    if (n <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
    } 
    int ierr = pIIR_->apply(n, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
    }
    return 0;
}

bool IIRFilter::isInitialized(void) const
{
    int len = pIIR_->isInitialized();
    return len;
}
