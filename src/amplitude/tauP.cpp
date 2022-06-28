#include <stdexcept>
#include <vector>
#include "rtseis/amplitude/tauP.hpp"
#include "rtseis/amplitude/tauPParameters.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterImplementations/detrend.hpp"
#include "rtseis/filterImplementations/iirFilter.hpp"
#include "rtseis/filterImplementations/taper.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"
#include "rtseis/filterImplementations/sosFilter.hpp"
#include "rtseis/utilities/math/convolve.hpp"

using namespace RTSeis::Amplitude;
namespace RConvolve = RTSeis::Utilities::Math::Convolve;

template<RTSeis::ProcessingMode E, class T>
class TauP<E, T>::TauPImpl
{
public:
    // Optionally highpass then integrate to velocity with trapezoid rule 
    // convolved with a RC filter
    RTSeis::FilterImplementations::SOSFilter<E, T> mIntegratorFilter;
    // Lowpass filter the velocity signals
    RTSeis::FilterImplementations::SOSFilter<E, T> mVelocityFilter;
    // Detrend (post-processing)
    RTSeis::FilterImplementations::Detrend<T> mDetrender;
    // Taper (post-processing)
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
    auto samplingRate = parameters.getSamplingRate();
    auto samplingPeriod = 1./samplingRate;
    constexpr RTSeis::FilterDesign::SOSPairing pairing
        = RTSeis::FilterDesign::SOSPairing::NEAREST;
    if (velocityFilter.getNumberOfNumeratorCoefficients() > 0 &&
        velocityFilter.getNumberOfDenominatorCoefficients() > 0)
    {
        auto sos = RTSeis::FilterDesign::IIR::tf2sos(velocityFilter,
                                                     pairing);
        pImpl->mVelocityFilter.initialize(sos);
    }
    // Setup acceleration filters
    if (parameters.getInputUnits() == InputUnits::Acceleration)
    {
        // IIR integrator (trapezoid rule)
        auto q = parameters.getFilterConstantQ();
        auto b = 2/(1 + q);
        auto bg = samplingPeriod/(b*pImpl->mGain);
        std::vector<double> bIntegrator{bg, bg};
        std::vector<double> aIntegrator{1, -q}; 
 
/*
        RTSeis::FilterRepresentations::BA baIntegrator;
        auto accelerationFilter = parameters.getAccelerationFilter();
        if (accelerationFilter.getNumberOfNumeratorCoefficients() > 0 &&
            accelerationFilter.getNumberOfDenominatorCoefficients() > 0)
        {
            auto bHighPass = accelerationFilter.getNumeratorCoefficients();
            auto aHighPass = accelerationFilter.getDenominatorCoefficients();
            // Convolve filters
            auto bWork = RConvolve::convolve(bHighPass, bIntegrator,
                                             RConvolve::FULL,
                                             RConvolve::Implementation::DIRECT);
            auto aWork = RConvolve::convolve(aHighPass, aIntegrator,
                                             RConvolve::FULL,
                                             RConvolve::Implementation::DIRECT);
            baIntegrator.setNumeratorCoefficients(bWork);
            baIntegrator.setDenominatorCoefficients(aWork);
        }
        else
        {
            baIntegrator.setNumeratorCoefficients(bIntegrator);
            baIntegrator.setDenominatorCoefficients(aIntegrator);
        }
        auto sosIntegrator = RTSeis::FilterDesign::IIR::tf2sos(baIntegrator,
                                                               pairing);
        pImpl->mIntegratorFilter.initialize(sosIntegrator);
        pImpl->mVelocity = false;
*/
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

