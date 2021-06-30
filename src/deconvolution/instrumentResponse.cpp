#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cmath>
#include "rtseis/deconvolution/instrumentResponse.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/sos.hpp"

using namespace RTSeis::Deconvolution;

class InstrumentResponse::InstrumentResponseImpl
{
public:
    //class RTSeis::FilterRepresentations::ZPK mZPK;
    RTSeis::FilterRepresentations::BA mBA;
    //class RTSeis::FilterRepresentations mSOS;
    double mSamplingRate = 0;
    bool mHaveResponse = false;
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
    pImpl->mHaveResponse = false;
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
void InstrumentResponse::setAnalogResponse(
    const RTSeis::FilterRepresentations::ZPK &zpk) noexcept
{
    pImpl->mIsAnalog = true;
    pImpl->mBA.clear();
    pImpl->mBA = RTSeis::FilterDesign::IIR::zpk2tf(zpk);
    pImpl->mHaveResponse = true;
}

/// Set a digital response
void InstrumentResponse::setDigitalResponse(
    const RTSeis::FilterRepresentations::ZPK &zpk) noexcept
{
    pImpl->mIsAnalog = false;
    pImpl->mBA.clear();
    pImpl->mBA = RTSeis::FilterDesign::IIR::zpk2tf(zpk);
    pImpl->mHaveResponse = true;
}

/*
void InstrumentResponse::setResponse(
    const RTSeis::FilterRepresentations::ZPK &zpk,
    const bool lanalog)
{
    pImpl->mIsAnalog = lanalog;
}
*/

/// Analog response?
bool InstrumentResponse::isAnalogResponse() const
{
    if (haveResponse()){throw std::runtime_error("Response not yet set");}
    return pImpl->mIsAnalog;
}

/// Have response?
bool InstrumentResponse::haveResponse() const noexcept
{
    return pImpl->mHaveResponse;
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
