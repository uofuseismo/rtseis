#include <cstdio>
#include <cmath>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/enums.hpp"
#include "private/throw.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterImplementations;

template<class T>
class FIRFilter<T>::FIRImpl
{
public:
    /// Default constructor
    FIRImpl()
    {
    }
    /// Copy constructor
    FIRImpl(const FIRImpl &fir)
    {
        *this = fir;
    }
    /// Default destructor.
    ~FIRImpl()
    {
        clear();
    }
    /// (Deep) copy operator
    FIRImpl& operator=(const FIRImpl &fir)
    {
        if (&fir == this){return *this;}
        if (!fir.linit_){return *this;}
        int ierr = initialize(fir.tapsLen_, fir.tapsRef_,
                              fir.mode_, fir.precision_,
                              fir.implementation_);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Initialization failed");
            clear();
            return *this;
        }
        // Now copy the filter states
        if (bufferSize_ > 0)
        {
            ippsCopy_8u(fir.pBuf_, pBuf_, bufferSize_);
        }   
        // Copy the initial conditions
        if (order_ > 0){ippsCopy_64f(fir.zi_, zi_, order_);}
        if (nwork_ > 0)
        {
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(fir.dlysrc64_, dlysrc64_, nwork_); 
                ippsCopy_64f(fir.dlydst64_, dlydst64_, nwork_);
            }
            else
            {
                ippsCopy_32f(fir.dlysrc32_, dlysrc32_, nwork_);
                ippsCopy_32f(fir.dlydst32_, dlydst32_, nwork_);
            }
        }
        return *this;
    }
    /// Clears memory off the module.
    void clear() noexcept
    {
        if (pSpec64_ != nullptr){ippsFree(pSpec64_);}
        if (pTaps64_ != nullptr){ippsFree(pTaps64_);}
        if (dlysrc64_ != nullptr){ippsFree(dlysrc64_);}
        if (dlydst64_ != nullptr){ippsFree(dlydst64_);}
        if (pSpec32_ != nullptr){ippsFree(pSpec32_);}
        if (pTaps32_ != nullptr){ippsFree(pTaps32_);}
        if (dlysrc32_ != nullptr){ippsFree(dlysrc32_);}
        if (dlydst32_ != nullptr){ippsFree(dlydst32_);}
        if (pBuf_ != nullptr){ippsFree(pBuf_);}
        if (tapsRef_ != nullptr){ippsFree(tapsRef_);}
        if (zi_ != nullptr){ippsFree(zi_);}
        pSpec64_ = nullptr;
        pTaps64_ = nullptr;
        dlysrc64_ = nullptr;
        dlydst64_ = nullptr;
        pSpec32_ = nullptr;
        pTaps32_ = nullptr;
        dlysrc32_ = nullptr;
        dlydst32_ = nullptr;
        pBuf_ = nullptr;
        tapsRef_ = nullptr;
        zi_ = nullptr;
        tapsLen_ = 0;
        nwork_ = 0;
        bufferSize_ = 0;
        specSize_ = 0;
        order_ = 0;
        implementation_ = FIRImplementation::DIRECT;
        precision_ = RTSeis::Precision::DOUBLE;
        mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        linit_ = false;
    }
    //========================================================================//
    /// Initializes the filter 
    int initialize(const int nb, const double b[],
                   const RTSeis::ProcessingMode mode,
                   const RTSeis::Precision precision,
                   const FIRImplementation implementation)
    {
        clear();
        // Figure out sizes and save some basic info
        tapsLen_ = nb;
        order_ = nb - 1;
        nwork_ = std::max(128, order_+1);
        tapsRef_ = ippsMalloc_64f(nb);
        ippsCopy_64f(b, tapsRef_, nb);
        zi_ = ippsMalloc_64f(std::max(1, order_));
        ippsZero_64f(zi_, std::max(1, order_));
        // Determine the algorithm type
        IppAlgType algType = ippAlgDirect;
        if (implementation == FIRImplementation::FFT){algType = ippAlgFFT;}
        if (implementation == FIRImplementation::AUTO){algType = ippAlgAuto;}
        // Initialize FIR filter
        if (precision == RTSeis::Precision::DOUBLE)
        {
            dlysrc64_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(dlysrc64_, nwork_);
            dlydst64_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(dlydst64_, nwork_);
            pTaps64_ = ippsMalloc_64f(tapsLen_);
            ippsCopy_64f(b, pTaps64_, tapsLen_);
            IppStatus status = ippsFIRSRGetSize(tapsLen_, ipp64f,
                                                &specSize_, &bufferSize_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error getting state size");
                clear();
                return -1; 
            }
            pSpec64_ = reinterpret_cast<IppsFIRSpec_64f *>
                       (ippsMalloc_8u(specSize_));
            pBuf_ = ippsMalloc_8u(bufferSize_);
            status = ippsFIRSRInit_64f(pTaps64_, tapsLen_,
                                       algType, pSpec64_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error initializing state structure");
                clear();
                return -1; 
            }
        }
        else
        {
            dlysrc32_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(dlysrc32_, nwork_);
            dlydst32_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(dlydst32_, nwork_);
            pTaps32_ = ippsMalloc_32f(tapsLen_);
            ippsConvert_64f32f(b, pTaps32_, tapsLen_);
            IppStatus status = ippsFIRSRGetSize(tapsLen_, ipp32f,
                                                &specSize_, &bufferSize_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error getting state size");
                clear();
                return -1; 
            }
            pSpec32_ = reinterpret_cast<IppsFIRSpec_32f *>
                       (ippsMalloc_8u(specSize_));
            pBuf_ = ippsMalloc_8u(bufferSize_);
            status = ippsFIRSRInit_32f(pTaps32_, tapsLen_,
                                       algType, pSpec32_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error initializing state structure");
                clear();
                return -1;
            }
        }
        implementation_ = implementation;
        precision_ = precision;
        mode_ = mode;
        linit_ = true;
        return 0;
    }
    /// Determines if the module is initialized.
    [[nodiscard]] bool isInitialized() const
    {
        return linit_;
    }
    /// Determines the length of the initial conditions.
    [[nodiscard]] int getInitialConditionLength() const
    {
        return order_;
    }
    /// Gets a copy of the initial conditions
    int getInitialConditions(const int nz, double zi[]) const
    {
        int nzRef = getInitialConditionLength();
        if (nz != nzRef){RTSEIS_WARNMSG("%s", "Shouldn't happen");}
        if (nzRef > 0){ippsCopy_64f(zi_, zi, nzRef);}
        return 0;
    }
    /// Sets the initial conditions
    int setInitialConditions(const int nz, const double zi[])
    {
        resetInitialConditions();
        int nzRef = getInitialConditionLength();
        if (nz != nzRef){RTSEIS_WARNMSG("%s", "Shouldn't happen");}
        if (nzRef > 0)
        {
            ippsCopy_64f(zi, zi_, nzRef);
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(zi_, dlysrc64_, nzRef);
            }
            else
            {
                ippsConvert_64f32f(zi_, dlysrc32_, nzRef);
            }
        }
        return 0;
    }
    /// Resets the initial conditions
    void resetInitialConditions() noexcept
    {
        if (order_ > 0)
        {   
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(zi_, dlysrc64_, order_);
            }
            else
            {
                ippsConvert_64f32f(zi_, dlysrc32_, order_);
            }
        }
    }
    /// Applies the filter
    int apply(const int n, const double x[], double y[])
    {
        if (n <= 0){return 0;} // Nothing to do
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
        IppStatus status = ippsFIRSR_64f(x, y, n, pSpec64_,
                                         dlysrc64_, dlydst64_, pBuf_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to apply FIR filter");
            return -1; 
        }   
        if (mode_ == RTSeis::ProcessingMode::REAL_TIME && order_ > 0)
        {
            ippsCopy_64f(dlydst64_, dlysrc64_, order_);
        }
        return 0;
    }
    /// Applies the filter
    int apply(const int n, const float x[], float y[])
    {
        if (n <= 0){return 0;} // Nothing to do
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
        IppStatus status = ippsFIRSR_32f(x, y, n, pSpec32_, 
                                         dlysrc32_, dlydst32_, pBuf_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to apply FIR filter");
            return -1;
        }
        if (mode_ == RTSeis::ProcessingMode::REAL_TIME && order_ > 0)
        {
            ippsCopy_32f(dlydst32_, dlysrc32_, order_);
        }
        return 0;
    }
private:
    /// The filter state.
    IppsFIRSpec_64f *pSpec64_ = nullptr;
    /// The filter taps.  This has dimension [tapsLen_].
    Ipp64f *pTaps64_ = nullptr;
    /// The input delay line.  This has dimension [nwork_].
    Ipp64f *dlysrc64_ = nullptr;
    /// The output delay line.  This has dimension [nwork_].
    Ipp64f *dlydst64_ = nullptr;
    /// The filter state.
    IppsFIRSpec_32f *pSpec32_ = nullptr;
    /// The filter taps.  This has dimension [tapsLen_].
    Ipp32f *pTaps32_ = nullptr;
    /// The input delay line.  This has dimension [nwork_].
    Ipp32f *dlysrc32_ = nullptr;
    /// The output delay line.  This has dimension [nwork_].
    Ipp32f *dlydst32_ = nullptr;
    /// Workspace.  This has dimension [bufferSize_].
    Ipp8u *pBuf_ = nullptr;
    /// A copy of the input taps. This has dimension [tapsLen_].
    double *tapsRef_ = nullptr;
    /// A copy of the initial conditions.  This has dimension [order_].
    double *zi_ = nullptr;
    /// The number of taps.
    int tapsLen_ = 0;
    /// The length of the delay line which is max(128, order+1).
    int nwork_ = 0;
    /// Size of pBuf.
    int bufferSize_ = 0;
    /// Size of state.
    int specSize_ = 0;
    /// Filter order.
    int order_ = 0;
    /// Implementation.
    FIRImplementation implementation_ = FIRImplementation::DIRECT;
    /// By default the module does post-procesing.
    RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
    /// The default module implementation.
    RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
    /// Flag indicating the module is initialized.
    bool linit_ = false;
};

//============================================================================//

/// C'tor
template<class T>
FIRFilter<T>::FIRFilter() :
    pFIR_(std::make_unique<FIRImpl> ())
{
}

/// Copy c'tor
template<class T>
FIRFilter<T>::FIRFilter(const FIRFilter &fir)
{
    *this = fir;
}

/// Move c'tor
template<class T>
FIRFilter<T>::FIRFilter(FIRFilter &&fir) noexcept
{
    *this = std::move(fir);
}

/// Destructor
template<class T>
FIRFilter<T>::~FIRFilter() = default;

/// Copy assignment
template<class T>
FIRFilter<T>& FIRFilter<T>::operator=(const FIRFilter &fir)
{
    if (&fir == this){return *this;}
    if (pFIR_){pFIR_->clear();}
    pFIR_ = std::make_unique<FIRImpl> (*fir.pFIR_);
    return *this;
}

/// Move assignment
template<class T>
FIRFilter<T>& FIRFilter<T>::operator=(FIRFilter &&fir) noexcept
{
    if (&fir == this){return *this;}
    if (pFIR_){pFIR_->clear();}
    pFIR_ = std::move(fir.pFIR_);
    return *this;
}

/// Initialization
template<>
void FIRFilter<double>::initialize(const int nb, const double b[],
                                   const RTSeis::ProcessingMode mode,
                                   FIRImplementation implementation)
{
    clear();
    // Checks
    if (nb < 1 || b == nullptr)
    {
        if (nb < 1){RTSEIS_THROW_IA("%s", "No b coefficients");}
        RTSEIS_THROW_IA("%s", "b is NULL");
    }
    constexpr RTSeis::Precision precision = RTSeis::Precision::DOUBLE;
#ifdef DEBUG
    int ierr = pFIR_->initialize(nb, b, mode, precision, implementation);
    assert(ierr == 0);
#else
    pFIR_->initialize(nb, b, mode, precision, implementation);
#endif
}

template<>
void FIRFilter<float>::initialize(const int nb, const double b[],
                                  const RTSeis::ProcessingMode mode,
                                  FIRImplementation implementation)
{
    clear();
    // Checks
    if (nb < 1 || b == nullptr)
    {
        if (nb < 1){RTSEIS_THROW_IA("%s", "No b coefficients");}
        RTSEIS_THROW_IA("%s", "b is NULL");
    }
    constexpr RTSeis::Precision precision = RTSeis::Precision::FLOAT;
#ifdef DEBUG
    int ierr = pFIR_->initialize(nb, b, mode, precision, implementation);
    assert(ierr == 0);
#else
    pFIR_->initialize(nb, b, mode, precision, implementation);
#endif
}

template<class T>
void FIRFilter<T>::clear() noexcept
{
    pFIR_->clear();
}
 
/// Initial conditions
template<class T>
void FIRFilter<T>::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    int nzRef = pFIR_->getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef){RTSEIS_THROW_IA("nz=%d should equal %d", nz, nzRef);}
        RTSEIS_THROW_IA("%s", "zi is NULL");
    }
    pFIR_->setInitialConditions(nz, zi);
}

template<class T>
void FIRFilter<T>::resetInitialConditions()
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    pFIR_->resetInitialConditions();
}

/// Filter application
template<class T>
void FIRFilter<T>::apply(const int n, const T x[], T *yIn[])
{
    if (n <= 0){return;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    T *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "Error x is NULL");}
        RTSEIS_THROW_IA("%s", "Error y is NULL");
    }
#ifdef DEBUG
    int ierr = pFIR_->apply(n, x, y);
    assert(ierr == 0);
#else
    pFIR_->apply(n, x, y);
#endif
}

/// Utility routine for initial conditon length
template<class T>
int FIRFilter<T>::getInitialConditionLength() const
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class is not yet initialized");
        return -1;
    }
    int len = pFIR_->getInitialConditionLength();
    return len;
}

/// Initialized?
template<class T>
bool FIRFilter<T>::isInitialized() const noexcept
{
    bool linit = pFIR_->isInitialized();
    return linit;
}

/// Get initial conditions
template<class T>
void FIRFilter<T>::getInitialConditions(const int nz, double *ziOut[]) const
{
    auto zi = *ziOut;
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class is not yet initialized");
    }
    int nzLen = getInitialConditionLength();
    if (nzLen > nz)
    {
        RTSEIS_THROW_IA("nz = %d must be at least %d", nz, nzLen);
    }
    pFIR_->getInitialConditions(nz, zi);
}

/// Template instantiation
template class RTSeis::Utilities::FilterImplementations::FIRFilter<double>;
template class RTSeis::Utilities::FilterImplementations::FIRFilter<float>;
