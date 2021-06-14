#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cmath>
#include "private/throw.hpp"
#include "rtseis/utilities/deconvolution/instrumentResponse.hpp"
#include "rtseis/utilities/filterDesign/iir.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/sos.hpp"

using namespace RTSeis::Utilities::Deconvolution;

class InstrumentResponse::InstrumentResponseImpl
{
public:
    //class RTSeis::FilterRepresentations::ZPK mZPK;
    class RTSeis::FilterRepresentations::BA mBA;
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

void InstrumentResponse::clear() noexcept
{
    pImpl->mBA.clear();
    pImpl->mSamplingRate = 0;
    pImpl->mHaveResponse = false;
    pImpl->mIsAnalog = true;
}

void InstrumentResponse::setSamplingRate(const double df)
{
    if (df <= 0)
    {
        RTSEIS_THROW_RTE("Sampling rate = %lf must be positive", df);
    }
    pImpl->mSamplingRate = df; 
}

void InstrumentResponse::setAnalogResponse(
    const RTSeis::FilterRepresentations::ZPK &zpk) noexcept
{
    pImpl->mIsAnalog = true;
    pImpl->mBA.clear();
    pImpl->mBA = RTSeis::Utilities::FilterDesign::IIR::zpk2tf(zpk);
    pImpl->mHaveResponse = true;
}

void InstrumentResponse::setDigitalResponse(
    const RTSeis::FilterRepresentations::ZPK &zpk) noexcept
{
    pImpl->mIsAnalog = false;
    pImpl->mBA.clear();
    pImpl->mBA = RTSeis::Utilities::FilterDesign::IIR::zpk2tf(zpk);
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

bool InstrumentResponse::isAnalogResponse() const
{
    if (haveResponse())
    {
        RTSEIS_THROW_RTE("%s", "Response not yet set");
    }
    return pImpl->mIsAnalog;
}

bool InstrumentResponse::haveResponse() const noexcept
{
    return pImpl->mHaveResponse;
}

/*
double InstrumentResponse::getSamplingRate() const
{
    if (pImpl->mSamplingRate <= 0)
    {
        RTSEIS_THROW_RTE("%s", "Sampling rate not yet set");
    }
    return pImpl->mSamplingRate;
}

double InstrumentResponse::getNyquistFrequency() const
{
    double fnyq = getSamplingRate()/2.0;
    return fnyq;
}
*/
