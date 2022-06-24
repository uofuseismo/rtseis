#include <stdexcept>
#include "rtseis/amplitude/tauP.hpp"
#include "rtseis/amplitude/tauPParameters.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "rtseis/filterImplementations/sosFilter.hpp"

using namespace RTSeis::Amplitude;

template<RTSeis::ProcessingMode E, class T>
class TauP<E, T>::TauPImpl
{
public:
    RTSeis::FilterImplementations::SOSFilter<E, T> mAccelerationFilter;
    RTSeis::FilterImplementations::SOSFilter<E, T> mVelocityFilter;
    bool mAcceleration{false};
};

/// C'tor
template<RTSeis::ProcessingMode E, class T>
TauP<E, T>::TauP() :
    pImpl(std::make_unique<TauPImpl> ())
{
}

/// Destructor
template<RTSeis::ProcessingMode E, class T>
TauP<E, T>::~TauP() = default;

/// Initialize
template<RTSeis::ProcessingMode E, class T>
void TauP<E, T>::initialize(const TauPParameters &parameters)
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
    if (parameters.getInputUnits() == InputUnits::Acceleration)
    {
        auto accelerationFilter = parameters.getAccelerationFilter();
        pImpl->mAccelerationFilter.initialize(accelerationFilter);
        pImpl->mAcceleration = true;
    }
    else
    {
        pImpl->mAcceleration = false;
    }
    auto velocityFilter = parameters.getVelocityFilter();    
    pImpl->mVelocityFilter.initialize(velocityFilter);
}
