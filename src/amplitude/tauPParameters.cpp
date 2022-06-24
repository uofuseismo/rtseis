#include <stdexcept>
#include <iostream>
#include <array>
#include <cmath>
#include "rtseis/amplitude/tauPParameters.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterRepresentations/sos.hpp"

using namespace RTSeis::Amplitude;

class TauPParameters::TauPParametersImpl
{
public:
    TauPParametersImpl()
    {
        update();
    }
    // Design a 0.075 Hz highpass filter.  
    // This is to be applied prior to integration.
    void designHighpassFilter()
    {
        constexpr int n = 4; // Brown et al., 2011 use n = 2 
                             // but it appears higher order filters
                             // result in less variance in tau p max
        constexpr auto highpass = RTSeis::FilterDesign::Bandtype::HIGHPASS;
        constexpr auto butterworth
            = RTSeis::FilterDesign::IIRPrototype::BUTTERWORTH;
        constexpr auto digital = RTSeis::FilterDesign::IIRFilterDomain::DIGITAL;
        constexpr auto pairing = RTSeis::FilterDesign::SOSPairing::NEAREST; 
        constexpr double rp = 5;
        constexpr double rs = 20; 
        const double nyquist = 0.5*mSamplingRate;
        const std::array<double, 2> Wn{0.075/nyquist, 0}; 
        mAccelerationFilter
            = FilterDesign::IIR::designSOSIIRFilter(n, Wn.data(), rp, rs,
                                                    highpass, butterworth,
                                                    digital, pairing);
    }
    // Design a 3 Hz highpass filter.  
    // This is to be applied prior to integration.
    void designLowpassFilter()
    {
        constexpr int n = 2; // Brown et al., 2011 use n = 2 
                             // but it appears higher order filters
                             // result in less variance in tau p max
        constexpr auto lowpass = RTSeis::FilterDesign::Bandtype::LOWPASS;
        constexpr auto butterworth
            = RTSeis::FilterDesign::IIRPrototype::BUTTERWORTH;
        constexpr auto digital = RTSeis::FilterDesign::IIRFilterDomain::DIGITAL;
        constexpr auto pairing = RTSeis::FilterDesign::SOSPairing::NEAREST;
        constexpr double rp = 5;
        constexpr double rs = 20; 
        const double nyquist = 0.5*mSamplingRate;
        const std::array<double, 2> Wn{3/nyquist, 0};
        mVelocityFilter
            = FilterDesign::IIR::designSOSIIRFilter(n, Wn.data(), rp, rs,
                                                    lowpass, butterworth,
                                                    digital, pairing);
    }
    // Update smoothing parameter - Brown et al., 2011
    void updateSmoothing()
    {
        mSmoothing = 1. - 1./mSamplingRate;
    }
    // For a new sampling rate updates the smoothing and lowpass filter
    void update()
    {
        designHighpassFilter();
        designLowpassFilter();
        updateSmoothing();
    }
//private:
    RTSeis::FilterRepresentations::SOS mAccelerationFilter;
    RTSeis::FilterRepresentations::SOS mVelocityFilter;
    double mSamplingRate{100};
    double mFilterQ{0.994};
    double mSmoothing{0.99}; //1. - 1./mSamplingRate;
    double mSimpleResponse{0}; // For scaling to units of cm to avoid overlflow
    InputUnits mInputUnits{InputUnits::Velocity};
    DetrendType mDetrendType{DetrendType::RemoveMean};
    bool mHaveInputUnits{false};
    bool mHaveSimpleResponse{false};
};

/// C'tor
TauPParameters::TauPParameters() :
    pImpl(std::make_unique<TauPParametersImpl> ())
{
}

/// Copy c'tor
TauPParameters::TauPParameters(const TauPParameters &parameters)
{
    *this = parameters;
}

/// Move c'tor
TauPParameters::TauPParameters(TauPParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Copy assignment
TauPParameters& TauPParameters::operator=(const TauPParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<TauPParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move assignment
TauPParameters& TauPParameters::operator=(TauPParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Reset class
void TauPParameters::clear() noexcept
{
    pImpl = std::make_unique<TauPParametersImpl> ();
}

/// Destructor
TauPParameters::~TauPParameters() = default;

/// Input units
void TauPParameters::setInputUnits(
    const InputUnits inputUnits) noexcept
{
    pImpl->mInputUnits = inputUnits;
    pImpl->mHaveInputUnits = true;
} 

InputUnits TauPParameters::getInputUnits() const
{
    if (!haveInputUnits()){throw std::runtime_error("Input units not set");}
    return pImpl->mInputUnits; 
}

bool TauPParameters::haveInputUnits() const noexcept
{
    return pImpl->mHaveInputUnits;
}

/// Sampling rate
void TauPParameters::setSamplingRate(const double samplingRate)
{
    if (samplingRate <= 0)
    {
        throw std::runtime_error("Sampling rate must be positive");
    }
    pImpl->mSamplingRate = samplingRate;
    pImpl->update();
}

double TauPParameters::getSamplingRate() const
{
    if (!haveSamplingRate()){throw std::runtime_error("Sampling rate not set");}
    return pImpl->mSamplingRate;
}

bool TauPParameters::haveSamplingRate() const noexcept
{
    return (pImpl->mSamplingRate > 0);
}

/// Simple response
void TauPParameters::setSimpleResponse(const double response)
{
    if (response == 0){throw std::runtime_error("Response cannot be 0");}
    pImpl->mSimpleResponse = response;
}

double TauPParameters::getSimpleResponse() const
{
    if (!haveSimpleResponse())
    {
        throw std::runtime_error("Simple response not set");
    }
    return pImpl->mSimpleResponse;
}

bool TauPParameters::haveSimpleResponse() const noexcept
{
    return (std::abs(pImpl->mSimpleResponse) > 0);
}

/// Acceleration filter
void TauPParameters::setAccelerationFilter(
    const RTSeis::FilterRepresentations::SOS &sos)
{
    pImpl->mAccelerationFilter = sos;
}

RTSeis::FilterRepresentations::SOS
    TauPParameters::getAccelerationFilter() const noexcept
{
    return pImpl->mAccelerationFilter;
}

/// Velocity filter
void TauPParameters::setVelocityFilter(
    const RTSeis::FilterRepresentations::SOS &sos)
{
    pImpl->mVelocityFilter = sos;
}

RTSeis::FilterRepresentations::SOS 
    TauPParameters::getVelocityFilter() const noexcept
{
    return pImpl->mVelocityFilter;
}

/// Smoothing 
void TauPParameters::setSmoothingParameter(const double alpha)
{
    if (std::abs(alpha) >= 1)
    {
        throw std::invalid_argument("|alpha| must be less than 1");
    }
    pImpl->mSmoothing = alpha;
}

double TauPParameters::getSmoothingParameter() const noexcept
{
    return pImpl->mSmoothing;
}

/// Q constant for RC Filter
void TauPParameters::setFilterConstantQ(const double q)
{
    if (std::abs(q) >= 1)
    {
        throw std::invalid_argument("|q| must be less than 1");
    }
    pImpl->mFilterQ = q;
}

double TauPParameters::getFilterConstantQ() const noexcept
{
    return pImpl->mFilterQ;
}

/// Detrend
void TauPParameters::setDetrendType(
    const DetrendType detrend)
{
    pImpl->mDetrendType = detrend;
}

DetrendType TauPParameters::getDetrendType() const noexcept
{
    return pImpl->mDetrendType;
}
