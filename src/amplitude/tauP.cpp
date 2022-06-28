#include <stdexcept>
#include "rtseis/amplitude/tauP.hpp"
#include "rtseis/amplitude/tauPParameters.hpp"
#include "rtseis/filterImplementations/detrend.hpp"
#include "rtseis/filterImplementations/taper.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "rtseis/filterImplementations/sosFilter.hpp"

using namespace RTSeis::Amplitude;

template<RTSeis::ProcessingMode E, class T>
class TauP<E, T>::TauPImpl
{
public:
    RTSeis::FilterImplementations::SOSFilter<E, T> mAccelerationFilter;
    RTSeis::FilterImplementations::SOSFilter<E, T> mVelocityFilter;
    RTSeis::FilterImplementations::Detrend<T> mDetrender;
    RTSeis::FilterImplementations::Taper<T> mTaperer;
    TauPParameters mParameters;
    double mGain{1};
    double mGainInverse{1};
    bool mDemean{false};
    bool mTaper{false};
    bool mVelocity{true};
    bool mInitialized{false};
};

/// C'tor
template<RTSeis::ProcessingMode E, class T>
TauP<E, T>::TauP() :
    pImpl(std::make_unique<TauPImpl> ())
{
}

/// Reset class
template<RTSeis::ProcessingMode E, class T>
void TauP<E, T>::clear() noexcept
{
    pImpl = std::make_unique<TauPImpl> ();
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
    pImpl->mGain = parameters.getSimpleResponse();
    pImpl->mGainInverse = 1./pImpl->mGain; 
    auto velocityFilter = parameters.getVelocityFilter();    
    pImpl->mVelocityFilter.initialize(velocityFilter);
    // Setup acceleration filters
    if (parameters.getInputUnits() == InputUnits::Acceleration)
    {
        auto accelerationFilter = parameters.getAccelerationFilter();
        pImpl->mAccelerationFilter.initialize(accelerationFilter);
        pImpl->mVelocity = false;
    }
    else
    {
        pImpl->mVelocity = true;
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

template<RTSeis::ProcessingMode E, class T>
bool TauP<E, T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

template<RTSeis::ProcessingMode E, class T>
bool TauP<E, T>::isVelocityFilter() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mVelocity;
}

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class RTSeis::Amplitude::TauP<RTSeis::ProcessingMode::POST, double>;
//template class RTSeis::Amplitude::TauP<RTSeis::ProcessingMode::REAL_TIME, double>;
//template class RTSeis::Amplitude::TauP<RTSeis::ProcessingMode::POST, float>;
//template class RTSeis::Amplitude::TauP<RTSeis::ProcessingMode::REAL_TIME, float>;

