#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#ifndef NDEBUG
#include <cassert>
#endif
#include <ipps.h>
#include "rtseis/enums.hpp"
#include "rtseis/filterImplementations/sosFilter.hpp"

using namespace RTSeis::FilterImplementations;

template<RTSeis::ProcessingMode E, class T>
class SOSFilter<E, T>::SOSFilterImpl
{
public:
    /// Default constructor
    SOSFilterImpl() = default;

    /// Copy constructor
    SOSFilterImpl(const SOSFilterImpl &sos)
    {
        *this = sos;
    }
    /// (Deep) copy operator
    SOSFilterImpl& operator=(const SOSFilterImpl &sos)
    {
        if (&sos == this){return *this;}
        if (!sos.mInitialized){return *this;}
        // Reinitialize the filter
        initialize(sos.nsections_, sos.bsRef_, sos.asRef_);
        // Now copy the filter states
        if (bufferSize_ > 0)
        {
            ippsCopy_8u(sos.pBuf_, pBuf_, bufferSize_);
        }
        // Copy the initial conditions
        if (nsections_ > 0){ippsCopy_64f(sos.zi_, zi_, 2*nsections_);}
        // And the delay lines
        if (nwork_ > 0)
        {
            if (mPrecision == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(sos.dlySrc64f_, dlySrc64f_, nwork_);
                ippsCopy_64f(sos.dlyDst64f_, dlyDst64f_, nwork_);
            }
            else
            {
                ippsCopy_32f(sos.dlySrc32f_, dlySrc32f_, nwork_); 
                ippsCopy_32f(sos.dlyDst32f_, dlyDst32f_, nwork_);
            }
        }
        return *this;
    }
    /// Default constructor
    ~SOSFilterImpl()
    {
        clear();
    }
    /// Clears the memory off the module
    void clear()
    {
        if (pTaps64f_ != nullptr){ippsFree(pTaps64f_);}
        if (dlySrc64f_ != nullptr){ippsFree(dlySrc64f_);}
        if (dlyDst64f_ != nullptr){ippsFree(dlyDst64f_);}
        if (pTaps32f_ != nullptr){ippsFree(pTaps32f_);}
        if (dlySrc32f_ != nullptr){ippsFree(dlySrc32f_);}
        if (dlyDst32f_ != nullptr){ippsFree(dlyDst64f_);}
        if (pBuf_ != nullptr){ippsFree(pBuf_);}
        if (bsRef_ != nullptr){ippsFree(bsRef_);}
        if (asRef_ != nullptr){ippsFree(asRef_);} 
        if (zi_ != nullptr){ippsFree(zi_);}
        pState64f_ = nullptr;
        pTaps64f_ = nullptr;
        dlySrc64f_ = nullptr;
        dlyDst64f_ = nullptr;
        pState32f_ = nullptr;
        pTaps32f_ = nullptr;
        dlySrc32f_ = nullptr;
        dlyDst32f_ = nullptr; 
        pBuf_ = nullptr;
        bsRef_ = nullptr;
        asRef_ = nullptr;
        zi_ = nullptr;
        nsections_ = 0;
        tapsLen_ = 0;
        nwork_ = 0;
        bufferSize_ = 0;
        mInitialized = false;
    }
    //========================================================================//
    int initialize(const int ns,
                   const double bs[],
                   const double as[])
    {
        clear();
        // Figure out sizes and copy the inputs
        nsections_ = ns;
        tapsLen_ = 6*nsections_;
        nwork_ = std::max(128, 2*nsections_);
        bsRef_ = ippsMalloc_64f(3*nsections_);
        ippsCopy_64f(bs, bsRef_, 3*nsections_);
        asRef_ = ippsMalloc_64f(3*nsections_);
        ippsCopy_64f(as, asRef_, 3*nsections_);
        zi_ = ippsMalloc_64f(2*nsections_);
        ippsZero_64f(zi_, 2*nsections_);
        IppStatus status;
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            status = ippsIIRGetStateSize_BiQuad_64f(nsections_,
                                                    &bufferSize_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to get state size" << std::endl;
                clear();
                return -1;
            }
            pBuf_ = ippsMalloc_8u(bufferSize_);
            pTaps64f_ = ippsMalloc_64f(tapsLen_);
            for (int i=0; i<nsections_; i++)
            {
                pTaps64f_[6*i+0] = bs[3*i+0];
                pTaps64f_[6*i+1] = bs[3*i+1];
                pTaps64f_[6*i+2] = bs[3*i+2];
                pTaps64f_[6*i+3] = as[3*i+0];
                pTaps64f_[6*i+4] = as[3*i+1];
                pTaps64f_[6*i+5] = as[3*i+2];
            }
            dlySrc64f_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(dlySrc64f_, nwork_);
            dlyDst64f_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(dlyDst64f_, nwork_);
            status = ippsIIRInit_BiQuad_64f(&pState64f_, pTaps64f_,
                                            nsections_,
                                            dlySrc64f_, pBuf_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to initialize biquad filter" << std::endl;
                clear();
                return -1;
            }
        }
        else
        {
            status = ippsIIRGetStateSize_BiQuad_32f(nsections_,
                                                    &bufferSize_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to get state size" << std::endl;
                clear();
                return -1; 
            }
            pBuf_ = ippsMalloc_8u(bufferSize_);
            pTaps32f_ = ippsMalloc_32f(tapsLen_);
            for (int i=0; i<nsections_; i++)
            {
                pTaps32f_[6*i+0] = static_cast<float> (bs[3*i+0]);
                pTaps32f_[6*i+1] = static_cast<float> (bs[3*i+1]);
                pTaps32f_[6*i+2] = static_cast<float> (bs[3*i+2]);
                pTaps32f_[6*i+3] = static_cast<float> (as[3*i+0]);
                pTaps32f_[6*i+4] = static_cast<float> (as[3*i+1]);
                pTaps32f_[6*i+5] = static_cast<float> (as[3*i+2]);
            }
            dlySrc32f_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(dlySrc32f_, nwork_);
            dlyDst32f_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(dlyDst32f_, nwork_);
            status = ippsIIRInit_BiQuad_32f(&pState32f_, pTaps32f_,
                                            nsections_,
                                            dlySrc32f_, pBuf_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to initialized biquad filter" << std::endl;
                clear();
                return -1;
            }
        }
        mInitialized = true;
        return 0;
    }
    /// Determines the length of the initial conditions
    [[nodiscard]] int getInitialConditionLength() const
    {
        return 2*nsections_;
    }
    /// Gets the number of sections
    [[nodiscard]] int getNumberOfSections() const
    {
        return nsections_;
    }
    /// Sets the initial conditions
    void setInitialConditions(const int nz, const double zi[]) noexcept
    {
        resetInitialConditions();
        int nzRef = getInitialConditionLength();
#ifndef NDEBUG
        assert(nz == nzRef);
        //if (nz != nzRef){std::cerr << "Shouldn't be here" << std::endl;}
#endif
        ippsCopy_64f(zi, zi_, nzRef);
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            ippsCopy_64f(zi_, dlySrc64f_, nzRef);
        }
        else
        {
            ippsConvert_64f32f(zi_, dlySrc32f_, nzRef);
        }
    }
    /// Resets the initial conditions
    void resetInitialConditions() noexcept
    {
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            ippsCopy_64f(zi_, dlySrc64f_, 2*nsections_);
        }
        else
        {
            ippsConvert_64f32f(zi_, dlySrc32f_, 2*nsections_);
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
        // Get a pointer to the filter state and set the initial conditions
        IppStatus status = ippsIIRSetDlyLine_64f(pState64f_, dlySrc64f_);
        if (status != ippStsNoErr)
        {
            std::cerr << "Failed to set delay line in double" << std::endl;
            return -1;
        }
        // Apply the filters
        status = ippsIIR_64f(x, y, n, pState64f_);
        if (status != ippStsNoErr)
        {
            std::cerr << "Failed to apply filter in double" << std::endl;
            return -1;
        }
        if (mMode == RTSeis::ProcessingMode::REAL_TIME)
        {
            status = ippsIIRGetDlyLine_64f(pState64f_, dlyDst64f_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to apply real-time filter in double"
                          << std::endl;
                return -1;
            }
            ippsCopy_64f(dlyDst64f_, dlySrc64f_, 2*nsections_);
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
                std::cerr << "Failed to apply double filter in float"
                          << std::endl;
                ippsFree(y64);
                return -1; 
            }
            ippsConvert_64f32f(y64, y, n); 
            ippsFree(y64);
            return 0;
        }
        // Get a pointer to the filter state and set the initial conditions
        IppStatus status = ippsIIRSetDlyLine_32f(pState32f_, dlySrc32f_);
        if (status != ippStsNoErr)
        {
            std::cerr << "Failed to set delay line in float" << std::endl;
            return -1;
        }
        // Apply the filters
        status = ippsIIR_32f(x, y, n, pState32f_);
        if (status != ippStsNoErr)
        {
            std::cerr << "Failed to apply filter in float" << std::endl;
            return -1;
        }
        if (mMode == RTSeis::ProcessingMode::REAL_TIME)
        {
            status = ippsIIRGetDlyLine_32f(pState32f_, dlyDst32f_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to get delay line" << std::endl;
                return -1;
            }
            ippsCopy_32f(dlyDst32f_, dlySrc32f_, 2*nsections_);
        }
        return 0;
    }
///private:
    IppsIIRState_64f *pState64f_ = nullptr;
    /// Filter taps.  This has dimension [tapsLen_].
    Ipp64f *pTaps64f_ = nullptr;
    /// Initial conditions. This has dimension [nwork_].
    Ipp64f *dlySrc64f_ = nullptr;
    /// Final conditions.  This has dimension [nwork_].
    Ipp64f *dlyDst64f_ = nullptr;
    /// Final conditions.  This has dimension [nwork_].
    /// Handle on filter state. 
    IppsIIRState_32f *pState32f_ = nullptr;
    /// Filter taps.  This has dimension [tapsLen_].
    Ipp32f *pTaps32f_ = nullptr;
    /// Initial conditions. This has dimension [nwork_].
    Ipp32f *dlySrc32f_ = nullptr;
    /// Final conditions.  This has dimension [nwork_].
    Ipp32f *dlyDst32f_ = nullptr;
    /// The workspace buffer.
    Ipp8u *pBuf_ = nullptr;
    /// A copy of the numerator filter coefficients.  This has
    /// dimension [3 x nsections_].
    double *bsRef_ = nullptr;
    /// A copy of the denominator filter coefficients.  This has
    /// dimension [3 x nsections_].
    double *asRef_ = nullptr; 
    /// A copy of the initial conditions.  This has dimension
    /// [2 x nsections_].
    double *zi_ = nullptr;
    /// The number of sections.
    int nsections_ = 0;
    /// The number of filter taps.  This equals 6*nsections_.
    int tapsLen_ = 0;
    /// Workspace for the delay lines.
    int nwork_ = 0;
    /// Size of workspace buffer.
    int bufferSize_ = 0;
    /// Real-time or post-processing.
    const RTSeis::ProcessingMode mMode = E;
    /// Single or double precision.
    const RTSeis::Precision mPrecision
        = (sizeof(T) == sizeof(double)) ? RTSeis::Precision::DOUBLE :
          RTSeis::Precision::FLOAT;
    /// Flag indicating the module is intiialized
    bool mInitialized = false;
};

//============================================================================//

/// C'tor
template<RTSeis::ProcessingMode E, class T>
SOSFilter<E, T>::SOSFilter() :
    pImpl(std::make_unique<SOSFilterImpl>())
{
}

/// Copy c'tor
template<RTSeis::ProcessingMode E, class T>
SOSFilter<E, T>::SOSFilter(const SOSFilter &sos)
{
    *this = sos;
}

/// Move c'tor
template<RTSeis::ProcessingMode E, class T>
SOSFilter<E, T>::SOSFilter(SOSFilter &&sos) noexcept
{
    *this = std::move(sos);
}

/// Destructor
template<RTSeis::ProcessingMode E, class T>
SOSFilter<E, T>::~SOSFilter() = default;

//// Clear the class
template<RTSeis::ProcessingMode E, class T>
void SOSFilter<E, T>::clear() noexcept
{
    pImpl->clear();
}

/// Copy assignment
template<RTSeis::ProcessingMode E, class T>
SOSFilter<E, T>& SOSFilter<E, T>::operator=(const SOSFilter &sos)
{
    if (&sos == this){return *this;}
    if (pImpl){pImpl->clear();}
    pImpl = std::make_unique<SOSFilterImpl> (*sos.pImpl);
    return *this;
}

/// Move assignment
template<RTSeis::ProcessingMode E, class T>
SOSFilter<E, T>& SOSFilter<E, T>::operator=(SOSFilter &&sos) noexcept
{
    if (&sos == this){return *this;}
    pImpl = std::move(sos.pImpl);
    return *this;
}

/// Initialization
template<RTSeis::ProcessingMode E, class T>
void SOSFilter<E, T>::initialize(const int ns,
                                 const double bs[],
                                 const double as[])
{
    clear();
    // Checks
    if (ns < 1 || bs == nullptr || as == nullptr)
    {
        if (ns < 1){throw std::invalid_argument("No sections");}
        if (bs == nullptr){throw std::invalid_argument("bs is NULL");}
        throw std::invalid_argument("as is NULL");
    }
    // Verify the highest order coefficients make sense
    for (auto i=0; i<ns; i++)
    {
        if (bs[3*i] == 0.0)
        {
            throw std::invalid_argument("Leading bs coefficient of section "
                                      + std::to_string(i) + " is zero");
        }
        if (as[3*i] == 0.0)
        {
            throw std::invalid_argument("Leading as coefficient of section "
                                      + std::to_string(i) + " is zero");
        }
    }
    auto ierr = pImpl->initialize(ns, bs, as);
#ifndef NDEBUG
    assert(ierr == 0);
#endif
    if (ierr != 0)
    {
        clear();
        throw std::runtime_error("Failed to initialize sos filter");
    }
}

/*
template<>
void SOSFilter<float>::initialize(const int ns,
                                  const double bs[],
                                  const double as[],
                                  const RTSeis::ProcessingMode mode)
{
    clear();
    // Checks
    if (ns < 1 || bs == nullptr || as == nullptr)
    {
        if (ns < 1){throw std::invalid_argument("No sections");}
        if (bs == nullptr){throw std::invalid_argument("bs is NULL");}
        throw std::invalid_argument("as is NULL");
    }
    // Verify the highest order coefficients make sense
    for (auto i=0; i<ns; i++)
    {
        if (bs[3*i] == 0.0)
        {
            throw std::invalid_argument("Leading bs coefficient of section " + std::to_string(i) + " is zero");
        }
        if (as[3*i] == 0.0)
        {
            throw std::invalid_argument("Leading as coefficient of section " + std::to_string(i) + " is zero");
        }
    }
    auto ierr = pImpl->initialize(ns, bs, as, mode, RTSeis::Precision::FLOAT);
#ifdef DEBUG
    assert(ierr == 0);
#endif
    if (ierr != 0)
    {
        clear();
        throw std::runtime_error("Failed to initialize sos filter");
    }
}
*/

/// Set initial conditions
template<RTSeis::ProcessingMode E, class T>
void SOSFilter<E, T>::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    resetInitialConditions();
    auto nzRef = pImpl->getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef)
        {
            auto errmsg = "nz = " + std::to_string(nz) + " must equal "
                        + std::to_string(nzRef);
            throw std::invalid_argument(errmsg);
        }
        throw std::invalid_argument("zi is NULL");
    }
    pImpl->setInitialConditions(nz, zi);
}

/// Reset initial conditions
template<RTSeis::ProcessingMode E, class T>
void SOSFilter<E, T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pImpl->resetInitialConditions();
}

/// Apply filter
template<RTSeis::ProcessingMode E, class T>
void SOSFilter<E, T>::apply(const int n, const T x[], T *yIn[]) 
{
    if (n <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    auto y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
#ifdef DEBUG
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

/*
template<>
void SOSFilter<float>::apply(const int n, const float x[], float *yIn[])
{
    if (n <= 0){return;}
    if (!isInitialized())
    {
        throw std::runtime_error("Class not initialized");
    }
    float *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::runtime_error("x is NULL");}
        throw std::runtime_error("y is NULL");
    }
#ifdef DEBUG
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
*/

/// Get initial conditions
template<RTSeis::ProcessingMode E, class T>
int SOSFilter<E, T>::getInitialConditionLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->getInitialConditionLength();
}

/// Get number of sections
template<RTSeis::ProcessingMode E, class T>
[[maybe_unused]]
int SOSFilter<E, T>::getNumberOfSections() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->getNumberOfSections();
}

/// Initialized?
template<RTSeis::ProcessingMode E, class T>
bool SOSFilter<E, T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

///--------------------------------------------------------------------------///
///                             Template Instantiation                       ///
///--------------------------------------------------------------------------///
template class RTSeis::FilterImplementations::SOSFilter<RTSeis::ProcessingMode::POST, double>;
template class RTSeis::FilterImplementations::SOSFilter<RTSeis::ProcessingMode::REAL_TIME, double>;
template class RTSeis::FilterImplementations::SOSFilter<RTSeis::ProcessingMode::POST, float>;
template class RTSeis::FilterImplementations::SOSFilter<RTSeis::ProcessingMode::REAL_TIME, float>;
