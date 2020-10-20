#include <iostream>
#include <cmath>
#include <string>
#ifndef NDEBUG
#include <cassert>
#endif
#include <ipps.h>
#include "rtseis/enums.hpp"
#include "rtseis/utilities/filterImplementations/iiriirFilter.hpp"

using namespace RTSeis::Utilities::FilterImplementations;

template<class T>
class IIRIIRFilter<T>::IIRIIRImpl
{
public:
    /// Default constructor
    /*
    IIRIIRImpl()
    {
    }
    /// Copy constructor
    IIRIIRImpl(const IIRIIRImpl &iiriir)
    {
        *this = iiriir;
    }
    */
    /// Destructor
    ~IIRIIRImpl()
    {
        clear();
    }
    /// (Deep) Copy operator
    IIRIIRImpl& operator=(const IIRIIRImpl &iiriir)
    {
        if (&iiriir == this){return *this;}
        if (!iiriir.linit_ ){return *this;}
        // Reinitialize the filter
        int ierr = initialize(iiriir.nbRef_, iiriir.bRef_,
                              iiriir.naRef_, iiriir.aRef_);
        if (ierr != 0)
        {
            std::cerr << "Failed to initialize filter in impl copy assignment"
                      << std::endl;
            clear();
            return *this;
        }
        // Now copy the filter states
        if (bufferSize_ > 0)
        {
            ippsCopy_8u(iiriir.pBuf_, pBuf_, bufferSize_);
        }
        // Copy the initial conditions
        if (nwork_ > 0)
        {
            ippsCopy_64f(iiriir.zi_, zi_, nwork_);
            if (iiriir.mPrecision == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(iiriir.dlysrc64_, dlysrc64_, nwork_);
            }
            else
            {
                ippsCopy_32f(iiriir.dlysrc32_, dlysrc32_, nwork_);
            }
        }
        lhaveZI_ = iiriir.lhaveZI_;
        return *this; 
    }
    /// Releases memory on the module.
    void clear()
    {
        if (pTaps64_ != nullptr){ippsFree(pTaps64_);}
        if (dlysrc64_ != nullptr){ippsFree(dlysrc64_);}
        if (pTaps32_ != nullptr){ippsFree(pTaps32_);}
        if (dlysrc32_ != nullptr){ippsFree(dlysrc32_);}
        if (pBuf_ != nullptr){ippsFree(pBuf_);}
        if (bRef_ != nullptr){ippsFree(bRef_);}
        if (aRef_ != nullptr){ippsFree(aRef_);}
        if (zi_ != nullptr){ippsFree(zi_);}
        pState64_ = nullptr;
        pTaps64_ = nullptr;
        dlysrc64_ =nullptr;
        pState32_ = nullptr;
        pTaps32_ = nullptr;
        dlysrc32_ = nullptr;
        pBuf_ = nullptr;
        bRef_ = nullptr;
        aRef_ = nullptr;
        zi_ = nullptr;
        nwork_ = 0;
        bufferSize_ = 0;
        order_ = 0;
        nbRef_ = 0;
        naRef_ = 0;
        lhaveZI_ = false;
        linit_ = false;
    }
    /// Sets the initial conditions
    [[nodiscard]]
    int initialize(const int nb, const double b[],
                   const int na, const double a[])
    {
        clear();
        // Copy basics
        nbRef_ = nb;
        naRef_ = na;
        order_ = std::max(nb, na) - 1;
        nwork_ = std::max(32, order_);
        bRef_ = ippsMalloc_64f(nb);
        ippsCopy_64f(b, bRef_, nb);
        aRef_ = ippsMalloc_64f(na);
        ippsCopy_64f(a, aRef_, na);
        zi_ = ippsMalloc_64f(nwork_);
        ippsZero_64f(zi_, nwork_);
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            // Workspace query 
            IppStatus status = ippsIIRIIRGetStateSize_64f(order_,
                                                          &bufferSize_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to get buffer size in double"
                          << std::endl;
                clear();
                return -1;
            }
            pTaps64_ = ippsMalloc_64f(2*(order_+1));
            pBuf_ = ippsMalloc_8u(bufferSize_);
            dlysrc64_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(dlysrc64_, nwork_);
            // Copy the filter taps 
            ippsZero_64f(pTaps64_, 2*(order_+1));
            ippsCopy_64f(b, &pTaps64_[0],        nb);
            ippsCopy_64f(a, &pTaps64_[order_+1], na);
            // Initialize the IIR state 
            status = ippsIIRIIRInit_64f(&pState64_, pTaps64_, order_,
                                        nullptr, pBuf_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to initialize IIRIIR state in double"
                          << std::endl;
                clear();
                return -1;
            }
        }
        else
        {
            // Workspace query 
            IppStatus status = ippsIIRIIRGetStateSize_32f(order_,
                                                          &bufferSize_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to get buffer size in float" << std::endl;
                clear();
                return -1; 
            }
            pTaps32_ = ippsMalloc_32f(2*(order_+1));
            pBuf_ = ippsMalloc_8u(bufferSize_);
            dlysrc32_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(dlysrc32_, nwork_);
            // Copy the filter taps 
            ippsZero_32f(pTaps32_, 2*(order_+1));
            ippsConvert_64f32f(b, &pTaps32_[0],        nb);
            ippsConvert_64f32f(a, &pTaps32_[order_+1], na);
            // Initialize the IIR state 
            status = ippsIIRIIRInit_32f(&pState32_, pTaps32_, order_,
                                        nullptr, pBuf_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to initialize IIRIIR state in float"
                          << std::endl;
                clear();
                return -1; 
            }
        }
        lhaveZI_ = false;
        linit_ = true;
        return 0;
    }
    /// Determines if the module is initialized
    [[nodiscard]] bool isInitialized() const
    {
        return linit_;
    }
    /// Gets the length of the initial conditions
    [[nodiscard]] int getInitialConditionLength() const
    {
        return order_;
    }
    /// Gets the filter order
    [[nodiscard]] [[maybe_unused]] int getFilterOrder() const
    {
        return order_;
    }
    /// Sets the intitial conditions length
    void setInitialConditions(const int nz, const double zi[])
    {
        resetInitialConditions();
        int nzRef = getInitialConditionLength();
#ifndef NDEBUG
        assert(nz == nzRef);
        //if (nz != nzRef){RTSEIS_WARNMSG("%s", "Shouldn't happen");}
#endif
        ippsCopy_64f(zi, zi_, nzRef);
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            ippsCopy_64f(zi, dlysrc64_, nzRef);
        }
        else
        {
            ippsConvert_64f32f(zi, dlysrc32_, nzRef);
        }
        lhaveZI_ = true;
    }
    /// Resets the initial conditions to those set in setInitialConditions.
    /// Note, the filter final coefficients are never extracted so the
    /// original filter initial conditions are already set.
    void resetInitialConditions() noexcept
    {
    }
    /// Applies the filter
    int apply(const int n, const double x[], double y[])
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
                std::cerr << "Failed to float apply filter in double"
                          << std::endl;
                ippsFree(y32);
                return -1;
            }
            ippsConvert_32f64f(y32, y, n);
            ippsFree(y32);
            return 0;
        }
        // Set a delay line if the user desires it.  Note, the
        // initialization sets a NULL delay line.
        IppStatus status;
        if (lhaveZI_)
        {
            status = ippsIIRIIRSetDlyLine_64f(pState64_, dlysrc64_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to set delay line in double" << std::endl;
                return -1;
            }
        }
        // Apply the filter
        status = ippsIIRIIR_64f(x, y, n, pState64_);
        if (status != ippStsNoErr)
        {
            std::cerr << "Failed to apply filter in double" << std::endl;
            return -1;
        }
        // Undo the action of setting a delay line
        if (lhaveZI_){ippsIIRIIRSetDlyLine_64f(pState64_, nullptr);}
        return 0; 
    }
    /// Applies the filter
    int apply(const int n, const float x[], float y[])
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
                std::cerr << "Failed to double apply filter in float"
                          << std::endl;
                ippsFree(y64);
                return -1;
            }
            ippsConvert_64f32f(y64, y, n);
            ippsFree(y64);
            return 0;
        }
        // Set a delay line if the user desires it.  Note, the
        // initialization sets a NULL delay line.
        IppStatus status;
        if (lhaveZI_)
        {
            status = ippsIIRIIRSetDlyLine_32f(pState32_, dlysrc32_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to set delay line in float" << std::endl;
                return -1;
            }
        }
        // Apply the filter
        status = ippsIIRIIR_32f(x, y, n, pState32_);
        if (status != ippStsNoErr)
        {
            std::cerr << "Failed to apply filter in double" << std::endl;
            return -1;
        }
        // Undo the action of setting a delay line
        if (lhaveZI_){ippsIIRIIRSetDlyLine_32f(pState32_, nullptr);}
        return 0;
    }
//private:
    /// The IIR filter state. 
    IppsIIRState_64f *pState64_ = nullptr;
    /// The IIR filter taps.  This has dimension [2*(order_+1)].
    Ipp64f *pTaps64_ = nullptr;
    /// The initial conditions.  This has dimension [nwork_]. 
    Ipp64f *dlysrc64_ = nullptr;
    /// The IIR filter state.
    IppsIIRState_32f *pState32_ = nullptr;
    /// The FIR filter taps.  This has dimension [2*(order_+1)].
    Ipp32f *pTaps32_ = nullptr;
    /// The initial conditions.  This has dimension [nwork_].
    Ipp32f *dlysrc32_ = nullptr;
    /// Workspace for the filter.  This has dimension [bufferSize_].
    Ipp8u *pBuf_ = nullptr;
    /// A saved copy of the numerator coefficients.
    /// This has dimension [nbRef_].
    double *bRef_ = nullptr;
    /// A saved copy of the denominator coefficients.
    /// This has dimension [naRef_].
    double *aRef_ = nullptr;
    /// A saved copy of the initial conditions. This has dimension [nwork_].
    double *zi_ = nullptr;
    /// The initial condtiions.  This has dimension [nwork_].
    int nwork_ = 0;
    /// The buffer size.
    int bufferSize_ = 0;
    /// The filter order.  This equals max(nbRef_, naRef_) - 1.
    int order_ = 0;
    /// The number of numerator coefficients.
    int nbRef_ = 0;
    /// The number of denominator coefficients.
    int naRef_ = 0;
    /// Flag indicating that the initial conditions have been set.
    bool lhaveZI_ = false;
    /// The default module implementation.
    /// Single or double precision.
    const RTSeis::Precision mPrecision
        = (sizeof(T) == sizeof(double)) ? RTSeis::Precision::DOUBLE :
          RTSeis::Precision::FLOAT;
    /// The default processing mode is post-processing only.
    //const RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
    /// Flag indicating that the filter is initialized.
    bool linit_ = false;
};

/// C'tor
template<class T>
IIRIIRFilter<T>::IIRIIRFilter() :
    pIIRIIR_(std::make_unique<IIRIIRImpl> ())
{
}

/// Copy c'tor
template<class T>
[[maybe_unused]]
IIRIIRFilter<T>::IIRIIRFilter(const IIRIIRFilter &iiriir)
{
    *this = iiriir;
}

/// Move c'tor
template<class T>
[[maybe_unused]]
IIRIIRFilter<T>::IIRIIRFilter(IIRIIRFilter &&iiriir) noexcept
{
    *this = std::move(iiriir);
}

/// Copy assignment
template<class T>
IIRIIRFilter<T>& IIRIIRFilter<T>::operator=(const IIRIIRFilter &iiriir)
{
    if (&iiriir == this){return *this;}
    if (pIIRIIR_){pIIRIIR_->clear();}
    pIIRIIR_ = std::make_unique<IIRIIRImpl>(*iiriir.pIIRIIR_);
    return *this;
}

/// Move assignment
template<class T>
IIRIIRFilter<T>& IIRIIRFilter<T>::operator=(IIRIIRFilter &&iiriir) noexcept
{
    if (&iiriir == this){return *this;}
    pIIRIIR_ = std::move(iiriir.pIIRIIR_);
    return *this;
}

/// Destructor
template<class T>
IIRIIRFilter<T>::~IIRIIRFilter()
{
    pIIRIIR_->clear();
}

/// Clears the class
template<class T>
void IIRIIRFilter<T>::clear() noexcept
{
    pIIRIIR_->clear();
}

/// Initialize the filter
template<>
void IIRIIRFilter<double>::initialize(const int nb, const double b[],
                                      const int na, const double a[])
{
    clear();
    // Check inputs
    if (nb < 1 || na < 1 || b == nullptr || a == nullptr)
    {
        if (nb < 1){throw std::invalid_argument("No b coefficients");}
        if (na < 1){throw std::invalid_argument("No a coefficients");}
        if (b == nullptr){throw std::invalid_argument("b is NULL");}
        throw std::invalid_argument("a is NULL");
    }
    if (a[0] == 0){throw std::invalid_argument("a[0] cannot equal 0");}
#ifndef NDEBUG
    int ierr = pIIRIIR_->initialize(nb, b, na, a);
    assert(ierr == 0);
#else
    int ierr = pIIRIIR_->initialize(nb, b, na, a);
    if (ierr != 0){throw std::runtime_error("Failed to initialize filter");}
#endif
}

/*
/// Initialize the filter (float)
template<>
void IIRIIRFilter<float>::initialize(const int nb, const double b[],
                                     const int na, const double a[])
{
    clear();
    // Check inputs
    if (nb < 1 || na < 1 || b == nullptr || a == nullptr)
    {
        if (nb < 1){RTSEIS_THROW_IA("%s", "No b coefficients");}
        if (na < 1){RTSEIS_THROW_IA("%s", "No a coefficients");}
        if (b == nullptr){RTSEIS_THROW_IA("%s", "b is NULL");}
        throw std::invalid_argument("a is NULL");
    }
    if (a[0] == 0)
    {
        throw std::invalid_argument("a[0] cannot equal 0");
    }
    constexpr RTSeis::Precision precision = RTSeis::Precision::FLOAT;
#ifndef NDEBUG
    int ierr = pIIRIIR_->initialize(nb, b, na, a, precision);
    assert(ierr == 0);
#else
    pIIRIIR_->initialize(nb, b, na, a, precision);
#endif
}
*/

/// Initial conditions
template<class T>
void IIRIIRFilter<T>::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int nzRef = pIIRIIR_->getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef)
        {
            auto errmsg = "nz = " + std::to_string(nz)
                + " must equal " + std::to_string(nzRef);
            throw std::invalid_argument(errmsg);
        }
        throw std::invalid_argument("zi is NULL");
    }
    pIIRIIR_->setInitialConditions(nz, zi);
}

/// Initialized?
template<class T>
bool IIRIIRFilter<T>::isInitialized() const noexcept
{
    bool linit = pIIRIIR_->isInitialized();
    return linit;
}

/// Apply
template<class T>
void IIRIIRFilter<T>::apply(const int n, const T x[], T *yIn[])
{
    if (n <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    T *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
#ifdef DEBUG
    int ierr = pIIRIIR_->apply(n, x, y);
    assert(ierr == 0);
#else
    pIIRIIR_->apply(n, x, y);
#endif
}

/// Reset initial conditions
template<class T>
void IIRIIRFilter<T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pIIRIIR_->resetInitialConditions();
}

/// Get initial condition length
template<class T>
int IIRIIRFilter<T>::getInitialConditionLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int len = pIIRIIR_->getInitialConditionLength();
    return len;
}

/// Get filter order
template<class T>
[[maybe_unused]]
int IIRIIRFilter<T>::getFilterOrder() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int len = pIIRIIR_->getFilterOrder();
    return len;
}

/// Template instantiation
template class RTSeis::Utilities::FilterImplementations::IIRIIRFilter<double>;
template class RTSeis::Utilities::FilterImplementations::IIRIIRFilter<float>;
