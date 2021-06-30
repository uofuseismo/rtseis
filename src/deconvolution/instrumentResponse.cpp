#include <cstdio>
#include <cstdlib>
#include <valarray>
#include <complex>
#include <cmath>
#include "rtseis/deconvolution/instrumentResponse.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterDesign/response.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/sos.hpp"

using namespace RTSeis::Deconvolution;

class InstrumentResponse::InstrumentResponseImpl
{
public:
    RTSeis::FilterRepresentations::BA mBA;
    double mSamplingRate = 0;
    bool mHaveTransferFunction = false;
    bool mIsAnalog = true;
};

InstrumentResponse::InstrumentResponse() :
    pImpl(std::make_unique<InstrumentResponseImpl> ())
{
}

/// Operators
InstrumentResponse&
InstrumentResponse::operator=(const InstrumentResponse &response)
{
    if (&response == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<InstrumentResponseImpl> (*response.pImpl);
    return *this;
}

InstrumentResponse&
InstrumentResponse::operator=(InstrumentResponse &&response) noexcept
{
    if (&response == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(response.pImpl);
    return *this;
}

/// Destructors
InstrumentResponse::~InstrumentResponse() = default;

/// Reset class
void InstrumentResponse::clear() noexcept
{
    pImpl->mBA.clear();
    pImpl->mSamplingRate = 0;
    pImpl->mHaveTransferFunction = false;
    pImpl->mIsAnalog = true;
}

/// Set sampling rate
void InstrumentResponse::setSamplingRate(const double df)
{
    if (df <= 0)
    {
        throw std::invalid_argument("Sampling rate = " + std::to_string(df)
                                  + " must be positive");
    }
    pImpl->mSamplingRate = df; 
}

/// Get sampling rate
double InstrumentResponse::getSamplingRate() const
{
    if (!haveSamplingRate())
    {   
        throw std::runtime_error("Sampling rate not yet set");
    }   
    return pImpl->mSamplingRate;
}

/// Have sampling rate
bool InstrumentResponse::haveSamplingRate() const noexcept
{
    return (pImpl->mSamplingRate > 0); 
}

/// Set an analog response
void InstrumentResponse::setAnalogTransferFunction(
    const RTSeis::FilterRepresentations::ZPK &zpk) noexcept
{
    auto ba = RTSeis::FilterDesign::IIR::zpk2tf(zpk);
    setAnalogTransferFunction(ba);
}

void InstrumentResponse::setAnalogTransferFunction(
     const RTSeis::FilterRepresentations::BA &ba)
{
    pImpl->mHaveTransferFunction = false;
    auto b = ba.getNumeratorCoefficients();
    auto a = ba.getDenominatorCoefficients();
    if (b.empty()){throw std::invalid_argument("No numerator coefficients");}
    if (a.empty()){throw std::invalid_argument("No demoninator coefficients");}
    bool okay = false;
    for (const auto &ai : a)
    {
        if (ai != 0){okay = true;}
    } 
    if (!okay){throw std::invalid_argument("All a coefficients are zero");} 
    pImpl->mIsAnalog = true;
    pImpl->mBA = ba;
    pImpl->mHaveTransferFunction = true;
}

/// Set a digital response
void InstrumentResponse::setDigitalTransferFunction(
    const RTSeis::FilterRepresentations::ZPK &zpk) noexcept
{
    auto ba = RTSeis::FilterDesign::IIR::zpk2tf(zpk);
    setDigitalTransferFunction(ba);
}

void InstrumentResponse::setDigitalTransferFunction(
     const RTSeis::FilterRepresentations::BA &ba)
{
    pImpl->mHaveTransferFunction = false;
    auto b = ba.getNumeratorCoefficients();
    auto a = ba.getDenominatorCoefficients();
    if (b.empty()){throw std::invalid_argument("No numerator coefficients");}
    if (a.empty()){throw std::invalid_argument("No demoninator coefficients");}
    bool okay = false;
    for (const auto &ai : a)
    {
        if (ai != 0){okay = true;}
    }
    if (!okay){throw std::invalid_argument("All a coefficients are zero");} 
    pImpl->mIsAnalog = false;
    pImpl->mBA = ba;
    pImpl->mHaveTransferFunction = true;
}

/*
void InstrumentResponse::setResponse(
    const RTSeis::FilterRepresentations::ZPK &zpk,
    const bool lanalog)
{
    pImpl->mIsAnalog = lanalog;
}
*/

std::vector<std::complex<double>>
InstrumentResponse::compute(const std::vector<double> &frequencies) const
{
    std::vector<std::complex<double>> response;
    if (frequencies.empty()){return response;}
    response.resize(frequencies.size());
    auto responsePtr = response.data();
    compute(frequencies.size(), frequencies.data(), &responsePtr); 
    return response;
}

void InstrumentResponse::compute(const int nFrequencies,
                                 const double frequencies[],
                                 std::complex<double> *responseIn[]) const
{
    if (nFrequencies < 1){return;}
    if (!haveTransferFunction())
    {
        throw std::runtime_error("Transfer function not yet set");
    }
    if (frequencies == nullptr)
    {
        throw std::invalid_argument("freqruencies is NULL");
    }
    auto response = *responseIn;
    if (response == nullptr){throw std::invalid_argument("response is NULL");}
    if (isAnalogTransferFunction())
    {
        constexpr double twopi = 2*M_PI;
        std::vector<double> omega(nFrequencies);
        std::transform(frequencies, frequencies + nFrequencies, omega.begin(),
                       [&twopi](auto &x){return twopi*x;});
        FilterDesign::Response::freqs(pImpl->mBA, omega);
    }
}
                        

/// Analog response?
bool InstrumentResponse::isAnalogTransferFunction() const
{
    if (!haveTransferFunction())
    {
        throw std::runtime_error("Transfer function not yet set");
    }
    return pImpl->mIsAnalog;
}

/// Have response?
bool InstrumentResponse::haveTransferFunction() const noexcept
{
    return pImpl->mHaveTransferFunction;
}

/*
double InstrumentResponse::getNyquistFrequency() const
{
    double fnyq = getSamplingRate()/2.0;
    return fnyq;
}
*/

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
