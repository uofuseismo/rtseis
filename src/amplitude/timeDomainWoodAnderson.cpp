#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include "rtseis/amplitude/timeDomainWoodAnderson.hpp"
#include "rtseis/amplitude/timeDomainWoodAndersonParameters.hpp"
#include "rtseis/filterImplementations/iirFilter.hpp"
#include "rtseis/filterImplementations/detrend.hpp"
#include "rtseis/filterImplementations/taper.hpp"
#include "rtseis/utilities/math/convolve.hpp"

using namespace RTSeis::Amplitude;
namespace RConvolve = RTSeis::Utilities::Math::Convolve;

namespace
{

}

/// Implemnetation
template<RTSeis::ProcessingMode E, class T>
class TimeDomainWoodAnderson<E, T>::TimeDomainWoodAndersonImpl
{
public:
    RTSeis::FilterImplementations::IIRFilter<E, T> mWAAccelerationFilter;
    RTSeis::FilterImplementations::IIRFilter<E, T> mWAVelocityFilter;
    RTSeis::FilterImplementations::Detrend<T> mDetrender;
    RTSeis::FilterImplementations::Taper<T> mTaperer;
    TimeDomainWoodAndersonParameters mParameters;
    bool mHighPassFilter = false;
    bool mVelocityFilter = true;
    bool mDemean = false;
    bool mTaper = false;
    bool mInitialized = false;
};

/// C'tor
template<RTSeis::ProcessingMode E, class T>
TimeDomainWoodAnderson<E, T>::TimeDomainWoodAnderson() :
    pImpl(std::make_unique<TimeDomainWoodAndersonImpl> ())
{
}

/// 
template<RTSeis::ProcessingMode E, class T>
void TimeDomainWoodAnderson<E, T>::initialize(
    const TimeDomainWoodAndersonParameters &parameters)
{
    if (!parameters.haveSamplingRate())
    {
        throw std::invalid_argument("Sampling rate not set");
    }
    if (!parameters.haveSimpleResponse())
    {
        throw std::invalid_argument("Simple response not set");
    }
    if (!parameters.haveInputUnits())
    {
        throw std::invalid_argument("Input units not set");
    }
    double df = parameters.getSamplingRate();
    double dt = 1./df;
    double h0 = parameters.getOptimizedDampingConstant();
    double f0 = parameters.getOptimizedNaturalAngularFrequency();
    double g0 = parameters.getOptimizedGain();
    double wdt = (2*M_PI)*(f0*dt);
    double gain = parameters.getSimpleResponse();
    double waCorrectedGain = 1;
    double c1 = 1 + h0*wdt;
    double c2 = 1 + 2*(h0*wdt) + wdt*wdt;
    
    if (parameters.getWoodAndersonGain() == WoodAndersonGain::WA_2080)
    {
        waCorrectedGain = 2080./2800.;
    }
    // Convolve in the high-pass filter so we get this in shot
    std::vector<double> bHighPass(1, 1);
    std::vector<double> aHighPass(1, 1);
    if (parameters.getHighPassFilter() == HighPassFilter::Yes)
    {
        double q = parameters.getHighPassFilterQ();
        bHighPass.resize(2);
        aHighPass.resize(2);
        bHighPass[0] = (1 + q)/2;
        bHighPass[1] =-(1 + q)/2;
        aHighPass[0] = 1;
        aHighPass[1] =-q;
    }
    // Create velocity filters
    if (parameters.getInputUnits() == InputUnits::Velocity)
    {
        pImpl->mVelocityFilter = true;
        double scalar = (g0/gain)*dt; 
        double b0 = waCorrectedGain*(scalar/c2);
        double b1 =-waCorrectedGain*(scalar/c2);
        double a1 =-((2*c1)/c2);
        double a2 = 1/c2;
        std::vector<double> b(3);
        b[0] = b0;
        b[1] = b1;
        b[2] = 0;
        std::vector<double> a(3);
        a[0] = 1;
        a[1] = a1;
        a[2] = a2;
        // Convolve filters 
        auto bFilter = RConvolve::convolve(bHighPass, b,
                                           RConvolve::FULL,
                                           RConvolve::Implementation::DIRECT);
        auto aFilter = RConvolve::convolve(aHighPass, a,
                                           RConvolve::FULL,
                                           RConvolve::Implementation::DIRECT);
        pImpl->mWAVelocityFilter.initialize(bFilter.size(), bFilter.data(),
                                            aFilter.size(), aFilter.data());
    }
    else
    {
        pImpl->mVelocityFilter = false;
        double scalar = (g0/gain)*(dt*dt);
        double b0 = waCorrectedGain*(scalar/c2);
        double b1 = 0;
        double a1 =-((2*c1)/c2);
        double a2 = 1/c2;
        std::vector<double> b(3);
        b[0] = b0;
        b[1] = b1;
        b[2] = 0;
        std::vector<double> a(3);
        a[0] = 1;
        a[1] = a1;
        a[2] = a2;
        // Convolve filters 
        auto bFilter = RConvolve::convolve(bHighPass, b,
                                           RConvolve::FULL,
                                           RConvolve::Implementation::DIRECT);
        auto aFilter = RConvolve::convolve(aHighPass, a,
                                           RConvolve::FULL,
                                           RConvolve::Implementation::DIRECT);
        pImpl->mWAAccelerationFilter.initialize(bFilter.size(), bFilter.data(),
                                                aFilter.size(), aFilter.data());
    }
    // Detrend?
    if constexpr (E == RTSeis::ProcessingMode::POST)
    {
        if (parameters.getDetrendType() == DetrendType::RemoveTrend)
        {
            pImpl->mDetrender.initialize(
                RTSeis::FilterImplementations::DetrendType::LINEAR);
            pImpl->mDemean = true;
        }
        else if (parameters.getDetrendType() == DetrendType::RemoveMean)
        {
            pImpl->mDetrender.initialize(
                RTSeis::FilterImplementations::DetrendType::CONSTANT);
            pImpl->mDemean = true;
        }
    } 
    // Taper
    if constexpr (E == RTSeis::ProcessingMode::POST)
    {
        auto percentage = parameters.getTaperPercentage();
        if (parameters.getWindowType() == WindowType::Sine &&
            percentage >= 0 && percentage <= 100)
        {
            pImpl->mTaperer.initialize(percentage,
                RTSeis::FilterImplementations::TaperWindowType::Sine);
            pImpl->mTaper = true;
        }
    }
    pImpl->mParameters = parameters;
    pImpl->mInitialized = true;
}

/// Initialized?
template<RTSeis::ProcessingMode E, class T>
bool TimeDomainWoodAnderson<E, T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Velocity filter?
template<RTSeis::ProcessingMode E, class T>
bool TimeDomainWoodAnderson<E, T>::isVelocityFilter() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mVelocityFilter;
}

/// Apply post-processing filter
template<RTSeis::ProcessingMode E, class T>
void TimeDomainWoodAnderson<E, T>::apply(const int n, const T x[], T *yPtr[])
{
    if (n < 1){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    if (*yPtr == nullptr){throw std::invalid_argument("y is NULL");}
    if constexpr (E == RTSeis::ProcessingMode::POST)
    {
        const T *xPtr = x;
        std::vector<T> xPrep; 
        // Do pre-processing
        if (pImpl->mDemean || pImpl->mTaper)
        {
            xPrep.resize(n);
            // Detrend
            if (pImpl->mDemean)
            {
                T *xPrepPtr = xPrep.data();
                pImpl->mDetrender.apply(n, x, &xPrepPtr);
            }
            else
            {
                std::copy(x, x + n, xPrep.data());
            }
            // Taper
            if (pImpl->mTaper)
            {
                std::vector<T> yPrep(n); 
                T *yPrepPtr = yPrep.data();
                pImpl->mTaperer.apply(n, xPrep.data(), &yPrepPtr);
                std::copy(yPrep.begin(), yPrep.end(), xPrep.data());
            }
            xPtr = xPrep.data();
        }
        // Do filtering
        if (isVelocityFilter())
        {
            pImpl->mWAVelocityFilter.apply(n, xPtr, yPtr);
        }
        else
        {
            pImpl->mWAAccelerationFilter.apply(n, xPtr, yPtr);
        }
    }
    else // Real time
    {
        if (isVelocityFilter())
        {
            pImpl->mWAVelocityFilter.apply(n, x, yPtr);
        }
        else
        {
            pImpl->mWAAccelerationFilter.apply(n, x, yPtr);
        }
    }
}

/// Reset initial filter conditions
template<RTSeis::ProcessingMode E, class T>
void TimeDomainWoodAnderson<E, T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (isVelocityFilter())
    {
        pImpl->mWAVelocityFilter.resetInitialConditions();
    }
    else
    {
        pImpl->mWAAccelerationFilter.resetInitialConditions();
    }
}
   

/// Destructor
template<RTSeis::ProcessingMode E, class T>
TimeDomainWoodAnderson<E, T>::~TimeDomainWoodAnderson() = default;

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class RTSeis::Amplitude::TimeDomainWoodAnderson<RTSeis::ProcessingMode::POST, double>;
template class RTSeis::Amplitude::TimeDomainWoodAnderson<RTSeis::ProcessingMode::REAL_TIME, double>;
template class RTSeis::Amplitude::TimeDomainWoodAnderson<RTSeis::ProcessingMode::POST, float>;
template class RTSeis::Amplitude::TimeDomainWoodAnderson<RTSeis::ProcessingMode::REAL_TIME, float>;
