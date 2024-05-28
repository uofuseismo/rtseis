#include <stdexcept>
#include <iostream>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#include <ippversion.h>
#include <ippcore.h>
#include <ipptypes.h>
#endif
#include "rtseis/enums.hpp"
#include "rtseis/filterImplementations/iirFilter.hpp"
#include "rtseis/filterRepresentations/ba.hpp"

using namespace RTSeis::FilterImplementations;

template<RTSeis::ProcessingMode E, class T>
class IIRFilter<E, T>::IIRFilterImpl
{
public:
    /// Default constructor
    IIRFilterImpl() = default;
    /// Copy constructor
    IIRFilterImpl(const IIRFilterImpl &iir)
    {
        *this = iir;
    }
    /// Default destructor
    ~IIRFilterImpl()
    {
        clear();
    }
    /// Copy operator
    IIRFilterImpl& operator=(const IIRFilterImpl &iir)
    {
        if (&iir == this){return *this;}
        clear();
        if (!iir.linit_){return *this;}
        int ierr = initialize(iir.nbRef_, iir.bRef_,
                              iir.naRef_, iir.aRef_,
                              iir.implementation_);
        if (ierr != 0)
        {
            std::cerr << "Failed to initialize filter in impl c'tor"
                      << std::endl;
            clear();
            return *this;
        }
        if (bufferSize_ > 0)
        {
            ippsCopy_8u(iir.pBuf_, pBuf_, bufferSize_);
        }
        ippsCopy_64f(iir.zi_, zi_, nbDly_);
        if (mPrecision == RTSeis::Precision::DOUBLE)
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
    void clear() noexcept
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
        pIIRState32f_ = nullptr;
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
        implementation_ = IIRDFImplementation::DF2_FAST;
        linit_ = false;
    }
    //========================================================================//
    /// Initializes the filter
    int initialize(const int nb, const double b[],
                   const int na, const double a[],
                   const IIRDFImplementation implementation)
    {
        clear();
        // May have to force the implementation to be slow
        IppStatus status;
        IIRDFImplementation impUse = implementation;
        // IPP changed implementation in 2018.  I'm not supporting it.
        // High order IIR filters aren't a good idea.
        #if IPP_VERSION_MAJOR < 2019
        impUse = IIRDFImplementation::DF2_SLOW;
        #endif
        if (impUse == IIRDFImplementation::DF2_FAST &&
            std::max(nb, na) - 1 > 8)
        {
            Ipp64u featureMask, enabledMask;
            status = ippGetCpuFeatures(&featureMask, nullptr);
            if (status == ippStsNoErr)
            {
                enabledMask = ippGetEnabledCpuFeatures();
                if ((featureMask & ippCPUID_AVX2) &&
                    (enabledMask & ippCPUID_AVX2))
                {
                    impUse = IIRDFImplementation::DF2_SLOW;
                }
                if ((featureMask & ippCPUID_AVX512F) &&
                    (enabledMask & ippCPUID_AVX512F))
                {
                    impUse = IIRDFImplementation::DF2_SLOW;
                }
            }
            else
            {
                std::cerr << "Failed to get processor info" << std::endl;
                return -1;
            }
        }
        // Actually a problem with order 0 filters and IPP
        if (std::max(nb, na) - 1 == 0){impUse = IIRDFImplementation::DF2_SLOW;}
        if (impUse != implementation)
        {
            std::cout << "Overriding IIR implementation to DF2_SLOW"
                      << std::endl;
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
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            if (impUse == IIRDFImplementation::DF2_FAST)
            {
                status = ippsIIRGetStateSize_64f(order_, &bufferSize_);
                if (status != ippStsNoErr)
                {
                    std::cerr << "Failed to get state size" << std::endl;
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
                    std::cerr << "Failed to initialize filter in double"
                              << std::endl;
                    clear();
                    return -1;
                }
                // Set the delay line
                status = ippsIIRSetDlyLine_64f(pIIRState64f_, pBufIPP64f_);
                if (status != ippStsNoErr)
                {
                    std::cerr << "Failed to set delay line in double"
                              << std::endl;
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
            if (impUse == IIRDFImplementation::DF2_FAST)
            {
                status = ippsIIRGetStateSize_32f(order_, &bufferSize_);
                if (status != ippStsNoErr)
                {
                    std::cerr << "Failed to get state size in float"
                              << std::endl;
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
                auto a04 = static_cast<float> (a0);
                ippsDivC_32f_I(a04, &pTaps32f_[0],        nbRef_);
                ippsDivC_32f_I(a04, &pTaps32f_[order_+1], naRef_);
                // Initialize the filter
                status = ippsIIRInit_32f(&pIIRState32f_, pTaps32f_, order_,
                                         pDlySrc32f_, pBuf_);
                if (status != ippStsNoErr)
                {
                    std::cerr << "Failed to initialize float filter"
                              << std::endl;
                    clear();
                    return -1;
                }
                // Set the delay line
                status = ippsIIRSetDlyLine_32f(pIIRState32f_, pBufIPP32f_);
                if (status != ippStsNoErr)
                {
                    std::cerr << "Failed to set delay line in float"
                              << std::endl;
                    clear();
                    return -1;
                }
            }
            else
            {
                auto a04 = static_cast<float> (a0);
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
        linit_ = true;
        return 0;
    }
    /// Determines if the filter is initialized
    [[nodiscard]] bool isInitialized() const noexcept
    {
        return linit_;
    }
    /// Gets the length of the initial conditions
    [[nodiscard]] int getInitialConditionLength() const
    {
        return order_;
    }
    /// Sets the initial conditions
    void setInitialConditions(const int nz, const double zi[]) noexcept
    {
        resetInitialConditions();
        int nzRef = getInitialConditionLength();
#ifndef NDEBUG
        assert(nzRef == nz);
#endif
        if (nzRef == 0){return;}
        ippsCopy_64f(zi, zi_, nzRef);
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            ippsCopy_64f(zi_, pBufIPP64f_, nzRef); 
            ippsIIRSetDlyLine_64f(pIIRState64f_, pBufIPP64f_);
        }
        else
        {
            ippsConvert_64f32f(zi_, pBufIPP32f_, nzRef);
            ippsIIRSetDlyLine_32f(pIIRState32f_, pBufIPP32f_);
        }
    }
    /// Resets the initial conditions
    void resetInitialConditions() noexcept
    {
        if (mPrecision == RTSeis::Precision::DOUBLE)
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
    }
    /// Applies the filter
    [[nodiscard]] int apply(const int n, const double x[], double y[])
    {
        if (n <= 0){return 0;}
        if (mPrecision == RTSeis::Precision::FLOAT)
        {
            Ipp32f *x32 = ippsMalloc_32f(n);
            Ipp32f *y32 = ippsMalloc_32f(n);
            ippsConvert_64f32f(x, x32, n);
            int ierr = apply(n, x32, y32);
            ippsFree(x32);
            if (ierr != 0)
            {
                std::cerr << "Failed to apply float filter in double"
                          << std::endl;
                ippsFree(y32);
                return -1;
            }
            ippsConvert_32f64f(y32, y, n);
            ippsFree(y32);
            return 0;
        }
        IppStatus status;
        if (implementation_ == IIRDFImplementation::DF2_FAST)
        {
            status = ippsIIR_64f(x, y, n, pIIRState64f_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to set delay line in double"
                          << std::endl;
                return -1;
            }
            if (mMode == RTSeis::ProcessingMode::REAL_TIME)
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
    [[nodiscard]] int apply(const int n, const float x[], float y[])
    {
        if (n <= 0){return 0;}
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            Ipp64f *x64 = ippsMalloc_64f(n);
            Ipp64f *y64 = ippsMalloc_64f(n);
            ippsConvert_32f64f(x, x64, n); 
            int ierr = apply(n, x64, y64);
            ippsFree(x64);
            if (ierr != 0)
            {
                std::cerr << "Failed to apply float double in float"
                          << std::endl;
                ippsFree(y64);
                return -1; 
            }
            ippsConvert_64f32f(y64, y, n); 
            ippsFree(y64);
            return 0;
        }
        IppStatus status;
        if (implementation_ == IIRDFImplementation::DF2_FAST)
        {
            status = ippsIIR_32f(x, y, n, pIIRState32f_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to set delay line in float"
                          << std::endl;
                return -1;
            }
            if (mMode == RTSeis::ProcessingMode::REAL_TIME)
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
        if (mMode == RTSeis::ProcessingMode::POST)
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
        if (mMode == RTSeis::ProcessingMode::POST)
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
    /// Holds the input IIR filter delay line.  This has dimension [nbDly_].
    Ipp64f *pDlySrc64f_ = nullptr;
    /// Holds the output IIR filter delay line.  This has dimension [nbDly_].
    Ipp64f *pDlyDst64f_ = nullptr;
    /// IIR filtering state
    IppsIIRState_32f *pIIRState32f_ = nullptr;;
    /// The Filter taps.  This has dimension [2*(order_+1)].
    Ipp32f *pTaps32f_ = nullptr;
    /// Holds the input IIR filter conditions.  This has dimension [bufIPPLen_].
    Ipp32f *pBufIPP32f_ = nullptr;
    /// Holds the input IIR filter delay line.  This has dimension [nbDly_].
    Ipp32f *pDlySrc32f_ = nullptr;
    /// Holds the output IIR filter delay line.  This has dimension [nbDly_].
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
    IIRDFImplementation implementation_ = IIRDFImplementation::DF2_FAST;
    /// Real-time vs. post-processing.
    const RTSeis::ProcessingMode mMode = E;
    /// Precision of filter application
    const RTSeis::Precision mPrecision
        = (sizeof(T) == sizeof(double)) ? RTSeis::Precision::DOUBLE :
          RTSeis::Precision::FLOAT;
    /// Flag indicating the module is initialized
    bool linit_ = false;
};

//============================================================================//

/// C'tor
template<RTSeis::ProcessingMode E, class T>
IIRFilter<E, T>::IIRFilter() :
    pImpl(std::make_unique<IIRFilterImpl> ())
{
}

/// Copy c'tor
template<RTSeis::ProcessingMode E, class T>
[[maybe_unused]]
IIRFilter<E, T>::IIRFilter(const IIRFilter &iir)
{
    *this = iir;
}

/// Move c'tor
template<RTSeis::ProcessingMode E, class T>
[[maybe_unused]]
IIRFilter<E, T>::IIRFilter(IIRFilter &&iir) noexcept
{
    *this = std::move(iir);
}

/// Copy assignment
template<RTSeis::ProcessingMode E, class T>
IIRFilter<E, T>& IIRFilter<E, T>::operator=(const IIRFilter &iir)
{
    if (&iir == this){return *this;}
    if (pImpl){pImpl->clear();}
    pImpl = std::make_unique<IIRFilterImpl> (*iir.pImpl);
    return *this;
}

/// Move assignment
template<RTSeis::ProcessingMode E, class T>
IIRFilter<E, T>& IIRFilter<E, T>::operator=(IIRFilter &&iir) noexcept
{
    if (&iir == this){return  *this;}
    pImpl = std::move(iir.pImpl);
    return *this;
}

/// Destructor
template<RTSeis::ProcessingMode E, class T>
IIRFilter<E, T>::~IIRFilter()
{
    clear();
}

/// Resets the class
template<RTSeis::ProcessingMode E, class T>
void IIRFilter<E, T>::clear() noexcept
{
    pImpl->clear();
}

/// Initialization
template<RTSeis::ProcessingMode E, class T>
void IIRFilter<E, T>::initialize(const RTSeis::FilterRepresentations::BA &ba,
                                 const IIRDFImplementation implementation)
{
    auto b = ba.getNumeratorCoefficients();
    if (b.empty()){throw std::invalid_argument("No numerator coefficients");}
    auto a = ba.getDenominatorCoefficients();
    if (a.empty()){throw std::invalid_argument("No denominator coefficients");}
    initialize(b.size(), b.data(),
               a.size(), a.data(),
               implementation);
}

/// Initialization
template<RTSeis::ProcessingMode E, class T>
void IIRFilter<E, T>::initialize(const int nb, const double b[],
                                 const int na, const double a[],
                                 const IIRDFImplementation implementation)
{
    clear();
    if (nb < 1 || b == nullptr || na < 1 || a == nullptr)
    {
        if (nb < 1){throw std::invalid_argument("No numerator coefficients");}
        if (na < 1)
        {
            throw std::invalid_argument("No denominator coefficients");
        }
        if (b == nullptr){throw std::invalid_argument("b is NULL");}
        throw std::invalid_argument("a is NULL");
    }
    if (a[0] == 0){throw std::invalid_argument("a[0] cannot be zero");}
#ifndef NDEBUG
    int ierr = pImpl->initialize(nb, b, na, a, implementation);
    assert(ierr == 0);
#else
    pImpl->initialize(nb, b, na, a, implementation);
#endif
}

/// Initial conditions
template<RTSeis::ProcessingMode E, class T>
int IIRFilter<E, T>::getInitialConditionLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->getInitialConditionLength();
}

template<RTSeis::ProcessingMode E, class T>
void IIRFilter<E, T>::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int nzRef = pImpl->getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef)
        {
            auto errmsg = " nz = " + std::to_string(nz)
                        + " must equal " + std::to_string(nzRef);
            throw std::invalid_argument(errmsg);
        }
        throw std::invalid_argument("zi is NULL");
    }
    pImpl->setInitialConditions(nz, zi);
}

template<RTSeis::ProcessingMode E, class T>
void IIRFilter<E, T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pImpl->resetInitialConditions();
}

template<RTSeis::ProcessingMode E, class T>
void IIRFilter<E, T>::apply(const int n, const T x[], T *yIn[])
{
    if (n <= 0){return;} // Nothing to do
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    T *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    } 
#ifndef NDEBUG
    int ierr = pImpl->apply(n, x, y);
    assert(ierr == 0);
#else
    auto error = pImpl->apply(n, x, y);
    if (error != 0)
    {
        throw std::runtime_error("Failed to apply filter");
    }
#endif
}

template<RTSeis::ProcessingMode E, class T>
bool IIRFilter<E, T>::isInitialized() const noexcept
{
    return pImpl->isInitialized();
}

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class RTSeis::FilterImplementations::IIRFilter<RTSeis::ProcessingMode::POST, double>;
template class RTSeis::FilterImplementations::IIRFilter<RTSeis::ProcessingMode::REAL_TIME, double>;
template class RTSeis::FilterImplementations::IIRFilter<RTSeis::ProcessingMode::POST, float>;
template class RTSeis::FilterImplementations::IIRFilter<RTSeis::ProcessingMode::REAL_TIME, float>;
