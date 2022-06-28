#include <stdexcept>
#include <iostream>
#include <complex>
#include <array>
#include <cmath>
#include "rtseis/amplitude/tauPParameters.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"
#include "rtseis/utilities/math/convolve.hpp"

using namespace RTSeis::Amplitude;
namespace RConvolve = RTSeis::Utilities::Math::Convolve;

class TauPParameters::TauPParametersImpl
{
public:
    TauPParametersImpl()
    {
        update();
    }
    void designIntegratorAndConvolve()
    {
        // Design a trapezoid-rule integrator (in frequency domain)
        auto dt = 1./mSamplingRate;
        auto b = 2/(1 + mFilterQ);
        //std::vector<double> bIntegrator{dt/(2*b), dt/(2*b)};
        //std::vector<double> aIntegrator{1, mFilterQ};
        std::vector<std::complex<double>> zIntegrator{-1};
        std::vector<std::complex<double>> pIntegrator{-mFilterQ};
        //
        auto k = mAccelerationFilter.getGain();
        auto zHighPass = mAccelerationFilter.getZeros();
        auto pHighPass = mAccelerationFilter.getPoles();
        // Convolve filters
/*
        auto bFilter = RConvolve::convolve(zHighPass, zIntegrator,
                                           RConvolve::FULL,
                                           RConvolve::Implementation::DIRECT);
        auto aFilter = RConvolve::convolve(pHighPass, pIntegrator,
                                           RConvolve::FULL,
                                           RConvolve::Implementation::DIRECT);
*/
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
        constexpr double rp = 5;
        constexpr double rs = 20; 
        const double nyquist = 0.5*mSamplingRate;
        const std::array<double, 2> Wn{0.075/nyquist, 0}; 
        auto mAccelerationFilter
            = FilterDesign::IIR::designZPKIIRFilter(n, Wn.data(), rp, rs,
                                                    highpass, butterworth,
                                                    digital);
        
/*
        auto bHighpass = baHighpass.getNumeratorCoefficients();
        auto aHighpass = baHighpass.getDenominatorCoefficients();
        // Design a trapezoid-rule integrator
        auto dt = 1./mSamplingRate;
        auto b = 2/(1 + mFilterQ);
        std::vector<double> bIntegrator{dt/(2*b), dt/(2*b)};
        std::vector<double> aIntegrator{1, mFilterQ}; 
        // Convolve filters
        auto bFilter = RConvolve::convolve(bHighPass, bIntegrator,
                                           RConvolve::FULL,
                                           RConvolve::Implementation::DIRECT);
        auto aFilter = RConvolve::convolve(aHighPass, aIntegrator,
                                           RConvolve::FULL,
                                           RConvolve::Implementation::DIRECT);
        // Convert to SOS filter
        //ba2sos(baFilter, pairing);
        // Design SOS highpass filter
        mAccelerationFilter
            = FilterDesign::IIR::designSOSIIRFilter(n, Wn.data(), rp, rs,
                                                    highpass, butterworth,
                                                    digital, pairing);
*/
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
        constexpr double rp = 5;
        constexpr double rs = 20; 
        const double nyquist = 0.5*mSamplingRate;
        const std::array<double, 2> Wn{3/nyquist, 0};
        mVelocityFilter
            = FilterDesign::IIR::designBAIIRFilter(n, Wn.data(), rp, rs,
                                                   lowpass, butterworth,
                                                   digital);
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
        designIntegratorAndConvolve();
        designLowpassFilter();
        updateSmoothing();
    }
//private:
    RTSeis::FilterRepresentations::ZPK mAccelerationFilter;
    RTSeis::FilterRepresentations::BA mVelocityFilter;
    double mSamplingRate{100};
    double mFilterQ{0.994};
    double mSmoothing{0.99}; //1. - 1./mSamplingRate;
    double mSimpleResponse{0}; // For scaling to units of cm to avoid overlflow
    double mTaperPct{5};
    InputUnits mInputUnits{InputUnits::Velocity};
    WindowType mWindow{WindowType::Sine};
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

/*
/// Acceleration filter
void TauPParameters::setAccelerationFilter(
    const RTSeis::FilterRepresentations::BA &ba)
{
    pImpl->mAccelerationFilter = ba;
}

RTSeis::FilterRepresentations::BA
    TauPParameters::getAccelerationFilter() const noexcept
{
    return pImpl->mAccelerationFilter;
}
*/

/// Velocity filter
void TauPParameters::setVelocityFilter(
    const RTSeis::FilterRepresentations::BA &ba)
{
    pImpl->mVelocityFilter = ba;
}

RTSeis::FilterRepresentations::BA
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

/// Taper pct
void TauPParameters::setTaper(
    const double pct, const WindowType window)
{
    if (pct < 0 || pct > 100)
    {
        throw std::invalid_argument("Percentage = " + std::to_string(pct)
                                  + " must be in range [0,100]");
    }
    pImpl->mTaperPct = pct;
    pImpl->mWindow = window;
}

double TauPParameters::getTaperPercentage() const noexcept
{
    return pImpl->mTaperPct;
}

WindowType TauPParameters::getWindowType() const noexcept
{
    return pImpl->mWindow;
}


