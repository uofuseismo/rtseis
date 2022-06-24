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
    RTSeis::FilterImplementations::SOSFilter<E, T> mHighpassFilter;
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
    auto highpassFilter = parameters.getFilter();
    pImpl->mHighpassFilter.initialize(highpassFilter);
}
