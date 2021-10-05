#include <iostream>
#include <string>
#include <cassert>
#include <array>
#include "rtseis/amplitude/timeDomainWoodAndersonParameters.hpp"

using namespace RTSeis::Amplitude;

namespace
{
bool isValidSamplingRateWork(const double samplingRate,
                             const double tolerance = 1.e-4)
{
    const std::array<double, 8> samplingRates{100, 200, 40, 80, 50, 20, 250, 500};
    for (const auto df : samplingRates) 
    {
        if (std::abs(df - samplingRate) < tolerance){return true;}
    }
    return false;
}

/// For more information see:
/// https://gitlab.com/aqms-swg/aqms-jiggle/-/blob/master/src/org/trinet/filters/WAFilter.java
std::tuple<double, double, double> 
    getWoodAndersonConstants(const double samplingRate,
                             const InputUnits inputUnits)
{
    // defaults : "vsp with delta-t = 0.01 sec" = short period velocity
    // h0   = 0.781; f0  = 1.29; g0 = 2963.0;
    double h0 =  0.78262E+00;
    double f0 =  0.12886E+01;
    double g0 =  0.29694E+04; // 7 Hz
    if (inputUnits == InputUnits::Velocity)
    {
        if (samplingRate == 100)
        {
            // Velocity : "vsp with delta-t = 0.01 sec" = short period velocity
            //h0   = 0.781; f0  = 1.29; g0 = 2963.0; // from paper, like 10 Hz
            //h0 =  0.78067E+00; f0 =  0.12884E+01; g0 =  0.29635E+04; // 10 Hz
            h0 =  0.78262E+00; f0 =  0.12886E+01; g0 =  0.29694E+04; // 7 Hz
        }
        else if (samplingRate == 200) // NEW added 09/14/2007 aww 
        {
            //h0 =  0.79245E+00; f0 =  0.12691E+01; g0 =  0.28861E+04; // 10 Hz
            h0 =  0.79251E+00; f0 =  0.12692E+01; g0 =  0.28865E+04; // 7 Hz
        }
        else if (samplingRate == 40)
        {
            //h0   = 0.743; f0  = 1.342; g0 = 3186.0; //  from Hiroo 08/18/06, like 7 HZ value
            //h0 =  0.72958E+00; f0 =  0.13398E+01; g0 =  0.31401E+04; // 10 Hz
            h0 =  0.74306E+00; f0 =  0.13422E+01; g0 =  0.31859E+04; // 7 Hz
        }
        else if (samplingRate == 250)
        {
            h0 =  0.79384E+00; f0 =  0.12656E+01; g0 =  0.28695E+04; // VEL fmax=7 250sps WA2800
        }
        else if (samplingRate == 500)
        {
            h0 =  0.79702E+00; f0 =  0.12578E+01; g0 =  0.28349E+04; // VEL fmax=7 500sps WA2800
        }
        else if (samplingRate == 50) // NEW added 09/14/2007 aww
        {
            //h0 =  0.74985E+00; f0 =  0.13238E+01; g0 =  0.30929E+04; // 10 Hz
            h0 =  0.75820E+00; f0 =  0.13251E+01; g0 =  0.31203E+04; // 7 Hz
        }
        else if (samplingRate == 80) // NEW added 09/14/2007 aww
        {
            //h0 =  0.77412E+00; f0 =  0.12976E+01; g0 =  0.29996E+04; // 10 Hz
            h0 =  0.77721E+00; f0 =  0.12980E+01; g0 =  0.30093E+04; // 7 Hz
        }
        else if (samplingRate == 20)
        {
            // Velocity : "vbb with delta-t = 0.05 sec" = broad-band velocity
            //h0   = 0.568; f0  = 1.39; g0 = 3110.0; // from paper, like 10 Hz
            //h0 =  0.57034E+00; f0 =  0.13924E+01; g0 =  0.31175E+04; // 10 Hz
            h0 =  0.63382E+00; f0 =  0.14094E+01; g0 =  0.33662E+04; // 7 Hz
        }
        else
        {
#ifndef NDEBUG
            assert(false);
#endif
            throw std::runtime_error("Unhandled sampling rate in velocity");
        }
    }
    else if (inputUnits == InputUnits::Acceleration) // accelerometer values
    {
        if (samplingRate == 100)
        {
            // Acceleration : "lg with delta-t = 0.01 sec" = "low gain" acceleration
            //h0   = 0.781; f0  = 1.29; g0 = 2963.0; // from paper, like 10 Hz?
            //h0 =  0.78401E+00; f0 =  0.12873E+01; g0 =  0.29683E+04; // 10 Hz
            h0 =  0.78427E+00; f0 =  0.12876E+01; g0 =  0.29700E+04; // 7 Hz
        }
        else if (samplingRate == 200)
        {
            //BUG ? g0 in line below appears to be 2080 value for 10Hz not the 2800 value 09/14/2007
            //h0 =  0.793; f0  = 1.27; g0 = 2144.0; // from Hiroo 9/1/04
            //h0 =  0.79245E+00; f0 =  0.12691E+01; g0 =  0.28861E+04; // 10 Hz
            h0 =  0.79251E+00; f0 =  0.12692E+01; g0 =  0.28865E+04; // 7 Hz
        }
        else if (samplingRate == 40) 
        {
            //h0   = 0.754; f0  = 1.336; g0 = 3191.0; // from Hiroo 08/18/06, like 7 HZ 
            //h0 =  0.75220E+00; f0 =  0.13331E+01; g0 =  0.31774E+04; // 10 Hz
            h0 =  0.75422E+00; f0 =  0.13357E+01; g0 =  0.31913E+04; // 7 Hz
        }
        else if (samplingRate == 250) // NEW added 09/14/2007 aww
        {
            h0 =  0.79409E+00; f0 =  0.12654E+01; g0 =  0.28695E+04; // ACC fmax=7 250sps WA2800
        }
        else if (samplingRate == 500) // NEW added 09/14/2007 aww
        {
            h0 =  0.79711E+00; f0 =  0.12578E+01; g0 =  0.28349E+04; // ACC fmax=7 500sps WA2800
        }
        else if (samplingRate == 50) // NEW added 09/14/2007 aww
        {
            //h0 =  0.76394E+00; f0 =  0.13194E+01; g0 =  0.31152E+04; // 10 Hz
            h0 =  0.76517E+00; f0 =  0.13210E+01; g0 =  0.31235E+04; // 7 Hz
        }
        else if (samplingRate == 80) 
        {
            //h0   = 0.781;  f0  = 1.29; g0 =  2963.0; // from paper
            //h0   = 0.774;  f0  = 1.30; g0 =  2999.0; // from RT code
            //h0 =  0.77937E+00; f0 =  0.12958E+01; g0 =  0.30073E+04; // 10 Hz
            h0 =  0.77980E+00; f0 =  0.12963E+01; g0 =  0.30101E+04; // 7 Hz
        }
        else if (samplingRate == 20) // NEW added 09/14/2007 aww
        {
            //h0 =  0.67044E+00; f0 =  0.13636E+01; g0 =  0.32957E+04; // 10 Hz
            h0 =  0.68405E+00; f0 =  0.13829E+01; g0 =  0.34011E+04; // 7 Hz
        }
        else
        {
#ifndef NDEBUG
            assert(false);
#endif
            throw std::runtime_error("Unhandled sampling rate in acceleration");
        }
    }
    else
    {
#ifndef NDEBUG
        assert(false);
#endif
        throw std::runtime_error("Unhandled units");
    }
    return std::tuple<double, double, double> (h0, f0, g0);
}
}

class TimeDomainWoodAndersonParameters::TimeDomainWoodAndersonParametersImpl
{
public:
    void setConstants()
    {
        if (mHaveInputUnits && mSamplingRate > 0)
        {
            std::tie(mH, mF, mG) = getWoodAndersonConstants(mSamplingRate,
                                                            mInputUnits);
        }
    }
    double mSamplingRate = 0;
    double mTaperPct = 5;
    double mSimpleResponse = 0;
    double mHighPassFilterQ = 0.998; // Default q in 1998 paper
    double mH = 0;
    double mG = 0;
    double mF = 0;
    InputUnits mInputUnits = InputUnits::Velocity;
    WoodAndersonGain mWAGain = WoodAndersonGain::WA_2800;
    WindowType mWindow = WindowType::Sine;
    DetrendType mDetrendType = DetrendType::RemoveMean;
    HighPassFilter mHighPassFilter = HighPassFilter::No;
    bool mHaveInputUnits = false;
};

/// C'tor
TimeDomainWoodAndersonParameters::TimeDomainWoodAndersonParameters() :
    pImpl(std::make_unique<TimeDomainWoodAndersonParametersImpl> ())
{
}

/// Copy c'tor 
TimeDomainWoodAndersonParameters::TimeDomainWoodAndersonParameters(
    const TimeDomainWoodAndersonParameters &parameters)
{
    *this = parameters;
}

/// Move c'tor
TimeDomainWoodAndersonParameters::TimeDomainWoodAndersonParameters(
    TimeDomainWoodAndersonParameters &&parameters) noexcept
{
    *this = parameters;
}

/// Copy assigment
TimeDomainWoodAndersonParameters& TimeDomainWoodAndersonParameters::operator=(
    const TimeDomainWoodAndersonParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<TimeDomainWoodAndersonParametersImpl>
            (*parameters.pImpl);
    return *this;
}

/// Move assignment
TimeDomainWoodAndersonParameters& TimeDomainWoodAndersonParameters::operator=(
    TimeDomainWoodAndersonParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Simple response
void TimeDomainWoodAndersonParameters::setSimpleResponse(const double gain)
{
    if (gain == 0){throw std::invalid_argument("Simple response is zero");}
    pImpl->mSimpleResponse = gain;
}

double TimeDomainWoodAndersonParameters::getSimpleResponse() const
{
    if (!haveSimpleResponse())
    {
        throw std::runtime_error("Simple response not set");
    }
    return pImpl->mSimpleResponse;
}

bool TimeDomainWoodAndersonParameters::haveSimpleResponse() const noexcept
{
    return (pImpl->mSimpleResponse != 0);
}

/// Input units
void TimeDomainWoodAndersonParameters::setInputUnits(
    const InputUnits inputUnits) noexcept
{
    pImpl->mInputUnits = inputUnits;
    pImpl->mHaveInputUnits = true;
    pImpl->setConstants();
}

InputUnits TimeDomainWoodAndersonParameters::getInputUnits() const
{
    if (!haveInputUnits()){throw std::runtime_error("Input units not set");}
    return pImpl->mInputUnits;
}

bool TimeDomainWoodAndersonParameters::haveInputUnits() const noexcept
{
    return pImpl->mHaveInputUnits;
}

/// Sampling rate
void TimeDomainWoodAndersonParameters::setSamplingRate(const double samplingRate)
{
    if (!this->isSamplingRateSupported(samplingRate))
    {
        throw std::invalid_argument("Sampling rate: "
                                  + std::to_string(samplingRate)
                                  + " (Hz) is not supported");
    }
    pImpl->mSamplingRate = samplingRate;
    pImpl->setConstants();
}

double TimeDomainWoodAndersonParameters::getSamplingRate() const
{
    if (!haveSamplingRate()){throw std::runtime_error("Sampling rate not set");}
    return pImpl->mSamplingRate;
}

bool TimeDomainWoodAndersonParameters::haveSamplingRate() const noexcept
{
    return (pImpl->mSamplingRate > 0);
}

bool TimeDomainWoodAndersonParameters::isSamplingRateSupported(
    const double samplingRate) const noexcept
{
    return isValidSamplingRateWork(samplingRate);
}

/// Taper pct
void TimeDomainWoodAndersonParameters::setTaper(
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

/// Taper pct
double TimeDomainWoodAndersonParameters::getTaperPercentage() const noexcept
{
    return pImpl->mTaperPct;
}

WindowType TimeDomainWoodAndersonParameters::getWindowType() const noexcept
{
    return pImpl->mWindow;
}

/// Detrend
void TimeDomainWoodAndersonParameters::setDetrendType(
    const DetrendType detrend)
{
    pImpl->mDetrendType = detrend;
}

DetrendType TimeDomainWoodAndersonParameters::getDetrendType() const noexcept
{
    return pImpl->mDetrendType;
}

/// WA gain
void TimeDomainWoodAndersonParameters::setWoodAndersonGain(
    const WoodAndersonGain waGain) noexcept
{
    pImpl->mWAGain = waGain;
}

WoodAndersonGain 
    TimeDomainWoodAndersonParameters::getWoodAndersonGain() const noexcept
{
    return pImpl->mWAGain;
}

/// High pass filter q
void TimeDomainWoodAndersonParameters::setHighPassRCFilter(
    const double q,
    HighPassFilter filter)
{
    if (std::abs(q) >= 1)
    {
        throw std::invalid_argument("|q| = " + std::to_string(q)
                                  + " must be less than 1");
    }
    pImpl->mHighPassFilterQ = q;
    pImpl->mHighPassFilter = filter;
}

double TimeDomainWoodAndersonParameters::getHighPassFilterQ() const noexcept
{
    return pImpl->mHighPassFilterQ;
}

HighPassFilter
    TimeDomainWoodAndersonParameters::getHighPassFilter() const noexcept
{
    return pImpl->mHighPassFilter;
}

/// Reset class
void TimeDomainWoodAndersonParameters::clear() noexcept
{
    pImpl = std::make_unique<TimeDomainWoodAndersonParametersImpl> ();
}

/// Destructor
TimeDomainWoodAndersonParameters::~TimeDomainWoodAndersonParameters() = default;

/// Optimized coefficients 
double TimeDomainWoodAndersonParameters::getOptimizedDampingConstant() const
{
    if (!haveSamplingRate()){throw std::runtime_error("Sampling rate not set");}
    if (!haveInputUnits()){throw std::runtime_error("Input units not set");}
    return pImpl->mH;
}

double
   TimeDomainWoodAndersonParameters::getOptimizedNaturalAngularFrequency() const
{
    if (!haveSamplingRate()){throw std::runtime_error("Sampling rate not set");}
    if (!haveInputUnits()){throw std::runtime_error("Input units not set");}
    return pImpl->mF;
}

double TimeDomainWoodAndersonParameters::getOptimizedGain() const
{
    if (!haveSamplingRate()){throw std::runtime_error("Sampling rate not set");}
    if (!haveInputUnits()){throw std::runtime_error("Input units not set");}
    return pImpl->mG;
}

