#include <stdexcept>
#include <iostream>
#ifndef NDEBUG
#include <cassert>
#endif
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/enums.hpp"
#include "rtseis/filterImplementations/downsample.hpp"

using namespace RTSeis::FilterImplementations;

template<RTSeis::ProcessingMode E, class T>
class Downsample<E, T>::DownsampleImpl
{
public:
    /// Default constructor
    DownsampleImpl() = default;

    /// Copy constructor
    DownsampleImpl(const DownsampleImpl &downsample)
    {
        *this = downsample;
    }
    /// (Deep) copy operator
    DownsampleImpl& operator=(const DownsampleImpl &downsample)
    {
        if (&downsample == this){return *this;}
        phase0_ = downsample.phase0_;
        downFactor_ = downsample.downFactor_;
        phase_ = downsample.phase_;
        //mode_ = downsample.mode_;
        //precision_ = downsample.precision_;
        mInitialized = downsample.mInitialized;
        return *this;
    }
    /// Default destructor
    ~DownsampleImpl() = default;
    /// Clears memory and resets the module
    void clear()
    {
        phase0_ = 0;
        downFactor_ = 0;
        phase_ = 0;
        //mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        //precision_ = RTSeis::Precision::DOUBLE;
        mInitialized = false;
    }
    //--------------------------------------------------------------------//
    /// Initialize the filter
    void initialize(const int downFactor) noexcept
    //               const RTSeis::ProcessingMode mode)
    //               const RTSeis::Precision precision)
    {
        clear();
        downFactor_ = downFactor;
        phase0_ = 0;
        phase_ = 0;
        //mode_ = mode;
        //precision_ = precision;
        mInitialized = true;
    }
    /// Estimates the space required to hold the downsampled signal
    [[nodiscard]] int estimateSpace(const int n) const noexcept
    {
        int phase = 0;
        if (mMode ==  RTSeis::ProcessingMode::REAL_TIME){phase = phase_;}
        int pDstLen = (n + downFactor_ - 1 - phase)/downFactor_;
        return pDstLen;  
    }
    /// Gets the downsampling factor
    [[nodiscard]] [[maybe_unused]]
    int getDownsampleFactor() const noexcept
    {
        return downFactor_;
    }
    /// Sets the initial conditions
    void setInitialConditions(const int phase)
    {
        resetInitialConditions();
#ifndef NDEBUG
        assert(phase >= 0 && phase < downFactor_);
        /*
        if (phase < 0 || phase > downFactor_ - 1) 
        {
            return -1;
        }
        */
#endif
        phase0_ = phase;
        phase_  = phase;
    }
    /// Resets the initial conditions
    void resetInitialConditions() noexcept
    {
        phase_ = phase0_;
    }
    /// Apply the downsampler
    int apply(const int nx, const double x[],
              int *nyDown, double y[])
    {
        *nyDown = 0;
        if (nx == 0){return 0;} 
        // There's really nothing to do
        if (downFactor_ == 1)
        {
            ippsCopy_64f(x, y, nx);
            *nyDown = nx; 
            return 0;
        }
        // Apply downsampler
        int pDstLen = estimateSpace(nx);
        int phase = 0;
        int pEst = pDstLen;
        if (mMode == RTSeis::ProcessingMode::REAL_TIME){phase = phase_;}
        IppStatus status = ippsSampleDown_64f(x, nx, y, &pEst,
                                              downFactor_, &phase);
        if (status != ippStsNoErr)
        {
            std::cerr << "Failed to downsample signal in double" << std::endl;
            return -1; 
        }
        // TODO - this is hack'ish.  I think I need to try cacheing vectors
        // that are too small to downsample then hit them on the next
        // iteration. 
        // Weird unhandled exception
        if (pEst < 0 || pEst >= nx){pEst = pDstLen;}
        *nyDown = pEst;
        if (mMode == RTSeis::ProcessingMode::REAL_TIME){phase_ = phase;}
        return 0;
    }
    /// Apply the downsampler
    int apply(const int nx, const float x[],
              int *nyDown, float y[])
    {
        *nyDown = 0;
        if (nx == 0){return 0;}
        // There's really nothing to do
        if (downFactor_ == 1)
        {
            ippsCopy_32f(x, y, nx);
            *nyDown = nx;
            return 0;
        }
        // Apply downsampler
        int pDstLen = estimateSpace(nx);
        int phase = 0;
        int pEst = pDstLen;
        if (mMode == RTSeis::ProcessingMode::REAL_TIME){phase = phase_;}
        IppStatus status = ippsSampleDown_32f(x, nx, y, &pEst,
                                              downFactor_, &phase);
        if (status != ippStsNoErr)
        {
            std::cerr << "Failed to downsample signal" << std::endl;
            return -1;
        }
        // Weird unhandled exception
        if (pEst < 0 || pEst >= nx){pEst = pDstLen;}
        *nyDown = pEst;
        if (mMode == RTSeis::ProcessingMode::REAL_TIME){phase_ = phase;}
        return 0;
    }
//private:
    /// Initial conditions for phase.
    int phase0_ = 0;
    /// Downsampling factor.
    int downFactor_ = 0;
    /// The phase.
    int phase_ = 0;
    /// By default the module does post-procesing.
    const RTSeis::ProcessingMode mMode = E; // RTSeis::ProcessingMode::POST_PROCESSING;
    /// The default module implementation.
    //RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
    /// Flag indicating the module is initialized.
    bool mInitialized = false;
};

//============================================================================//
//                               End Implementations                          //
//============================================================================//

template<RTSeis::ProcessingMode E, class T>
Downsample<E, T>::Downsample() :
    pImpl(std::make_unique<DownsampleImpl>())
{
}

template<RTSeis::ProcessingMode E, class T>
Downsample<E, T>::~Downsample() = default;

template<RTSeis::ProcessingMode E, class T>
Downsample<E, T>::Downsample(const Downsample &downsample)
{
    *this = downsample;
}

template<RTSeis::ProcessingMode E, class T>
Downsample<E, T>::Downsample(Downsample &&downsample) noexcept
{
    *this = std::move(downsample); 
}

template<RTSeis::ProcessingMode E, class T>
Downsample<E, T>& Downsample<E, T>::operator=(const Downsample &downsample)
{
    if (&downsample == this){return *this;}
    pImpl = std::make_unique<DownsampleImpl> (*downsample.pImpl);
    //pImpl = std::unique_ptr<DownsampleImpl>
    //               (new DownsampleImpl(*downsample.pImpl));
    return *this;
}

template<RTSeis::ProcessingMode E, class T>
Downsample<E, T>& Downsample<E, T>::operator=(Downsample &&downsample) noexcept
{
    if (&downsample == this){return *this;}
    pImpl = std::move(downsample.pImpl);
    return *this; 
}

/// Initializers
template<RTSeis::ProcessingMode E, class T>
void Downsample<E, T>::initialize(const int downFactor)
{
    clear();
    if (downFactor < 1)
    {
        throw std::invalid_argument("Downsampling factor = "
                                  + std::to_string(downFactor)
                                  + " must be positive");
    }
    pImpl->initialize(downFactor); //, RTSeis::Precision::DOUBLE);
}

/*
template<>
void Downsample<float>::initialize(const int downFactor,
                                   const RTSeis::ProcessingMode mode)
{
    clear();
    if (downFactor < 1)
    {
        throw std::invalid_argument("Downsampling factor=%d must be positive", downFactor);
    }
    pImpl->initialize(downFactor, mode);//, RTSeis::Precision::FLOAT);
}
 */

template<RTSeis::ProcessingMode E, class T>
void Downsample<E, T>::clear() noexcept
{
    if (pImpl){pImpl->clear();}
}

template<RTSeis::ProcessingMode E, class T>
void Downsample<E, T>::setInitialConditions(const int phase)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int downFactor = pImpl->getDownsampleFactor();
    if (phase < 0 || phase > downFactor - 1)
    {
        throw std::invalid_argument("phase = " + std::to_string(phase)
                                  + " must be in range [0,"
                                  + std::to_string(downFactor-1) + "]");
    }
    pImpl->setInitialConditions(phase);
}

template<RTSeis::ProcessingMode E, class T>
void Downsample<E, T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pImpl->resetInitialConditions();
}

template<RTSeis::ProcessingMode E, class T>
void Downsample<E, T>::apply(const int nx, const T x[], const int ny,
                             int *nyDown, T *yIn[])
{
    *nyDown = 0;
    if (nx <= 0){return;} // Nothing to do
    if (!isInitialized())
    {
        throw std::runtime_error("Downsampler not intitialized");
    }
    T *y = *yIn;
    int pDstLen = pImpl->estimateSpace(nx); 
#ifndef NDEBUG
    assert(pDstLen >= 0);
#endif
    if (ny < pDstLen)
    {
        auto errmsg = "ny = " + std::to_string(ny)
                    + " must be at least length = "
                    + std::to_string(pDstLen);
        throw std::invalid_argument(errmsg);
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    int ierr = pImpl->apply(nx, x, nyDown, y);
#ifndef NDEBUG
    assert(ierr == 0);
#endif
    if (ierr != 0)
    {
        throw std::runtime_error("Failed to apply downsampler");
    }
}

template<RTSeis::ProcessingMode E, class T>
bool Downsample<E, T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

template<RTSeis::ProcessingMode E, class T>
int Downsample<E, T>::estimateSpace(const int n) const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (n < 0){throw std::invalid_argument("n cannot be negative");}
    return pImpl->estimateSpace(n);
}
 
template<RTSeis::ProcessingMode E, class T>
[[maybe_unused]]
int Downsample<E, T>::getDownsampleFactor() const noexcept
{
    //if (!isInitialized()){throw std::runtime_error("Class is not initialized");}
    return pImpl->getDownsampleFactor();
}

/// Template instantiation
template class RTSeis::FilterImplementations::Downsample<RTSeis::ProcessingMode::POST, double>;
template class RTSeis::FilterImplementations::Downsample<RTSeis::ProcessingMode::REAL_TIME, double>;
template class RTSeis::FilterImplementations::Downsample<RTSeis::ProcessingMode::POST, float>;
template class RTSeis::FilterImplementations::Downsample<RTSeis::ProcessingMode::REAL_TIME, float>;
