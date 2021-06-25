#include <iostream>
#include <ipps.h>
#ifndef NDEBUG
#include <cassert>
#endif
#include "rtseis/enums.hpp"
#include "rtseis/filterImplementations/medianFilter.hpp"

using namespace RTSeis::FilterImplementations;

template<RTSeis::ProcessingMode E, class T>
class MedianFilter<E, T>::MedianFilterImpl
{
public:
    /// Default constructor.
    MedianFilterImpl() = default;
    /// Copy constructor.
    MedianFilterImpl(const MedianFilterImpl &median)
    {
        *this = median;
    }
    /// (Deep) Copy operator
    MedianFilterImpl& operator=(const MedianFilterImpl &median)
    {
        if (&median == this){return *this;}
        if (!median.mInitialized){return *this;}
        // Reinitialize the filter
        int ierr = initialize(median.maskSize_);
        if (ierr != 0)
        {
            std::cerr << "Failed to initialize median filter in impl c'tor"
                      << std::endl;
            clear();
            return *this;
        }
        // Now copy the filter states
        if (bufferSize_ > 0)
        {
            ippsCopy_8u(median.pBuf_, pBuf_, bufferSize_);
        }
        if (nwork_ > 0)
        {
            ippsCopy_64f(median.zi_, zi_, nwork_);
            if (mPrecision == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(median.dlysrc64_, dlysrc64_, nwork_);
                ippsCopy_64f(median.dlydst64_, dlydst64_, nwork_);
            }
            else
            {
                ippsCopy_32f(median.dlysrc32_, dlysrc32_, nwork_);
                ippsCopy_32f(median.dlydst32_, dlydst32_, nwork_);
            }
        }
        return *this;
    }
    /// Destructor
    ~MedianFilterImpl()
    {
        clear();
    }
    //========================================================================//
    /// Clears the filter/releases the memory
    void clear() noexcept
    {
        if (dlysrc64_ != nullptr){ippsFree(dlysrc64_);}
        if (dlydst64_ != nullptr){ippsFree(dlydst64_);}
        if (dlysrc32_ != nullptr){ippsFree(dlysrc32_);}
        if (dlydst32_ != nullptr){ippsFree(dlydst32_);}
        if (pBuf_     != nullptr){ippsFree(pBuf_);}
        if (zi_       != nullptr){ippsFree(zi_);}
        dlysrc64_ = nullptr;
        dlydst64_ = nullptr;
        dlysrc32_ = nullptr;
        dlydst32_ = nullptr;
        pBuf_ = nullptr;
        zi_ = nullptr;
        maskSize_ = 0;
        bufferSize_ = 0;
        mInitialized = false;
    }
    /// Initializes the filter
    int initialize(const int n)
    {
        clear();
        maskSize_ = n; // This better be odd by this point
        // Set the space
        nwork_ = std::max(8, maskSize_ - 1);
        zi_ = ippsMalloc_64f(nwork_);
        ippsZero_64f(zi_, nwork_);
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            IppStatus status = ippsFilterMedianGetBufferSize(maskSize_,
                                                             ipp64f,
                                                             &bufferSize_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Error getting double buffer size" << std::endl;
                clear();
                return -1;
            }
            dlysrc64_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(dlysrc64_, nwork_);
            dlydst64_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(dlydst64_, nwork_);
            pBuf_ = ippsMalloc_8u(bufferSize_);
            ippsZero_8u(pBuf_, bufferSize_);
        }
        else
        {
            IppStatus status = ippsFilterMedianGetBufferSize(maskSize_,
                                                             ipp32f,
                                                             &bufferSize_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Error getting float buffer size" << std::endl;
                clear();
                return -1;
            }
            dlysrc32_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(dlysrc32_, nwork_);
            dlydst32_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(dlydst32_, nwork_);
            pBuf_ = ippsMalloc_8u(bufferSize_);
            ippsZero_8u(pBuf_, bufferSize_);
        }
        mInitialized = true;
        return 0;
    }
    /// Determines the length of the initial conditions.
    [[nodiscard]] int getInitialConditionLength() const
    {
        int len = maskSize_ - 1;
        return len;
    }
    /// Gets the group delay.
    [[nodiscard]] [[maybe_unused]] int getGroupDelay() const
    {
        int grpDelay = maskSize_/2;
        return grpDelay;
    }
    /// Set the initial conditions
    int setInitialConditions(const int nz, const double zi[])
    {
        resetInitialConditions();
        int nzRef = getInitialConditionLength();
#ifndef NDEBUG
        assert(nzRef == nz);
#endif
        if (nzRef > 0){ippsCopy_64f(zi, zi_, nzRef);}
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {   
            if (nzRef > 0){ippsCopy_64f(zi, dlysrc64_, nzRef);}
        }
        else
        {
            if (nzRef > 0){ippsConvert_64f32f(zi, dlysrc32_, nzRef);}
        }   
        return 0;
    }
    /// Resets the initial conditions
    int resetInitialConditions()
    {
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            if (nwork_ > 0){ippsCopy_64f(zi_, dlysrc64_, nwork_);}
        }
        else
        {
            if (nwork_ > 0){ippsConvert_64f32f(zi_, dlysrc32_, nwork_);}
        }   
        return 0;
    }
    /// Apply the filter
    int apply(const int n, const double x[], double y[])
    {
        if (n <= 0){return 0;} // Nothing to do
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
        if (mMode == RTSeis::ProcessingMode::REAL_TIME)
        {
            status = ippsFilterMedian_64f(x, y, n, maskSize_,
                                          dlysrc64_, dlydst64_, pBuf_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to apply real-time filter in double"
                          << std::endl;
                return -1;
            }
            if (maskSize_ > 1)
            {
                ippsCopy_64f(dlydst64_, dlysrc64_, maskSize_-1);
            }
        }
        else
        {
            status = ippsFilterMedian_64f(x, y, n, maskSize_,
                                          dlysrc64_, nullptr, pBuf_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to apply post-processing filter in double"
                          << std::endl;
                return -1;
            }
        }
        return 0;
    }
    /// Apply the filter
    int apply(const int n, const float x[], float y[])
    {
        if (n <= 0){return 0;} // Nothing to do
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
        IppStatus status;
        if (mMode == RTSeis::ProcessingMode::REAL_TIME)
        {
            status = ippsFilterMedian_32f(x, y, n, maskSize_,
                                          dlysrc32_, dlydst32_, pBuf_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to apply real-time filter in float"
                          << std::endl;
                return -1;
            }
            if (maskSize_ > 1)
            {
                ippsCopy_32f(dlydst32_, dlysrc32_, maskSize_-1);
            }
        }
        else
        {
            status = ippsFilterMedian_32f(x, y, n, maskSize_,
                                          dlysrc32_, nullptr, pBuf_);
            if (status != ippStsNoErr)
            {
                std::cerr << "Failed to apply post-processing filter in float"
                          << std::endl;
                return -1;
            }
        }
        return 0;
    }
//private:
    /// Delay line source vector.  This has dimension [nwork_].
    Ipp64f *dlysrc64_ = nullptr;
    /// Delay line destination vector.  This has dimension [nwork_].
    Ipp64f *dlydst64_ = nullptr;
    /// Delay line source vector.  This has dimension [nwork_].
    Ipp32f *dlysrc32_ = nullptr;
    /// Delay line destination vector.  This has dimension [nwork_].
    Ipp32f *dlydst32_ = nullptr;
    /// Workspace for median filter.  This has dimension [bufferSize_].
    Ipp8u *pBuf_ = nullptr;
    /// A reference of the saved initial conditions.  This has 
    /// dimension [nwork_] though only the first maskSize_  - 1
    /// points are valid.
    Ipp64f *zi_ = nullptr;
    /// The median filter window length.
    int maskSize_ = 0;
    /// The workspace for the delay lines.
    int nwork_ = 0;
    /// The size of the workspace buffer.
    int bufferSize_ = 0;
    /// Real-time vs. post-processing.
    const RTSeis::ProcessingMode mMode = E;
    /// The default module implementation.
    const RTSeis::Precision mPrecision
         = (sizeof(T) == sizeof(double)) ? RTSeis::Precision::DOUBLE :
                                           RTSeis::Precision::FLOAT;
    /// Flag indicating the module is initialized.
    bool mInitialized = false;
};

/// C'tor
template<RTSeis::ProcessingMode E, class T>
MedianFilter<E, T>::MedianFilter() :
    pImpl(std::make_unique<MedianFilterImpl> ())
{
}

/// Copy c'tor
template<RTSeis::ProcessingMode E, class T>
MedianFilter<E, T>::MedianFilter(const MedianFilter &median)
{
    *this = median;
}

/// Move c'tor
template<RTSeis::ProcessingMode E, class T>
[[maybe_unused]]
MedianFilter<E, T>::MedianFilter(MedianFilter &&median) noexcept
{
    *this = std::move(median);
}

/// Copy assignment
template<RTSeis::ProcessingMode E, class T>
MedianFilter<E, T>& MedianFilter<E, T>::operator=(const MedianFilter &median)
{
    if (&median == this){return *this;}
    if (pImpl){pImpl->clear();}
    pImpl = std::make_unique<MedianFilterImpl> (*median.pImpl);
    return *this;
}

/// Move assignment
template<RTSeis::ProcessingMode E, class T>
MedianFilter<E, T>&
MedianFilter<E, T>::operator=(MedianFilter &&median) noexcept
{
    if (&median == this){return *this;}
    pImpl = std::move(median.pImpl);
    return *this;
}

/// Destructor
template<RTSeis::ProcessingMode E, class T>
MedianFilter<E, T>::~MedianFilter() = default;

/// Clears the class
template<RTSeis::ProcessingMode E, class T>
void MedianFilter<E, T>::clear() noexcept
{
    pImpl->clear();
}

/// Initialization
template<RTSeis::ProcessingMode E, class T>
void MedianFilter<E, T>::initialize(const int n)
{
    clear();
    // Set the mask size
    if (n < 1)
    {
        auto errmsg = "Mask size = " + std::to_string(n)
                    + " must be positive";
        throw std::invalid_argument(errmsg);
    }
    int maskSize = n;
    if (maskSize%2 == 0)
    {
        maskSize = maskSize + 1;
        std::cout << "n = " << n << " should be odd; setting maskSize to "
                  << maskSize << std::endl;
    }
#ifndef NDEBUG
    int ierr = pImpl->initialize(maskSize);
    assert(ierr == 0);
#else
    pImpl->initialize(maskSize);
#endif
}

/// Sets the initial conditions
template<RTSeis::ProcessingMode E, class T>
void MedianFilter<E, T>::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int nzRef = pImpl->getInitialConditionLength();
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

template<RTSeis::ProcessingMode E, class T>
void MedianFilter<E, T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pImpl->resetInitialConditions();
}

template<RTSeis::ProcessingMode E, class T>
bool MedianFilter<E, T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

template<RTSeis::ProcessingMode E, class T>
int MedianFilter<E, T>::getInitialConditionLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    auto len = pImpl->getInitialConditionLength();
    return len;
}

/// Apply the filter
template<RTSeis::ProcessingMode E, class T>
void MedianFilter<E, T>::apply(const int n, const T x[], T *yIn[])
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
    pImpl->apply(n, x, y);
#endif
}

template<RTSeis::ProcessingMode E, class T>
[[maybe_unused]]
int MedianFilter<E, T>::getGroupDelay() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized\n");}
    return pImpl->getGroupDelay();
}

///--------------------------------------------------------------------------///
///                           Template Instantiation                         ///
///--------------------------------------------------------------------------///
template class RTSeis::FilterImplementations::MedianFilter<RTSeis::ProcessingMode::POST, double>;
template class RTSeis::FilterImplementations::MedianFilter<RTSeis::ProcessingMode::REAL_TIME, double>;
template class RTSeis::FilterImplementations::MedianFilter<RTSeis::ProcessingMode::POST, float>;
template class RTSeis::FilterImplementations::MedianFilter<RTSeis::ProcessingMode::REAL_TIME, float>;
