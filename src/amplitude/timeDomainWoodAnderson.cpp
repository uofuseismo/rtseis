#include <cmath>
#include <array>
#include "rtseis/amplitude/timeDomainWoodAnderson.hpp"
#include "rtseis/amplitude/timeDomainWoodAndersonParameters.hpp"
#include "rtseis/filterImplementations/iirFilter.hpp"

using namespace RTSeis::Amplitude;

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
    TimeDomainWoodAndersonParameters mParameters;
    bool mVelocityFilter = true;
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
    double dt = 1/df;
    double h0 = parameters.getOptimizedDampingConstant();
    double f0 = parameters.getOptimizedNaturalAngularFrequency();
    double g0 = parameters.getOptimizedGain();
    double wdt = (2*M_PI)*(f0*dt);
    double gain = parameters.getSimpleResponse();
    double waCorrectedGain = 1;
    if (parameters.getWoodAndersonGain() == WoodAndersonGain::WA_2080)
    {
        waCorrectedGain = 2080./2800.;
    }
    if (parameters.getInputUnits() == InputUnits::Velocity)
    {
        pImpl->mVelocityFilter = true;
        double scalar = (g0/gain)*dt; 
        double c1 = 1 + h0*wdt;
        double c2 = 1 + 2*(h0*wdt) + wdt*wdt;
        double b0 = waCorrectedGain*(1/c2);
        double b1 =-waCorrectedGain*(scalar/c2);
        double a0 =-waCorrectedGain*((2*c1)/c2);
        double a1 = waCorrectedGain*(1/c2);
        const std::array<double, 2> b{b0, b1};
        const std::array<double, 2> a{a0, a1};
        pImpl->mWAVelocityFilter.initialize(b.size(), b.data(),
                                            a.size(), a.data());
    }
    else
      {
        pImpl->mVelocityFilter = false;
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


/// Destructor
template<RTSeis::ProcessingMode E, class T>
TimeDomainWoodAnderson<E, T>::~TimeDomainWoodAnderson() = default;

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class RTSeis::Amplitude::TimeDomainWoodAnderson<RTSeis::ProcessingMode::POST, double>;
template class RTSeis::Amplitude::TimeDomainWoodAnderson<RTSeis::ProcessingMode::REAL_TIME, double>;
//template class RTSeis::Amplitude::TimeDomainWoodAnderson<RTSeis::ProcessingMode::POST, float>;
//template class RTSeis::Amplitude::TimeDomainWoodAnderson<RTSeis::ProcessingMode::REAL_TIME, float>;
