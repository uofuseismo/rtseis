#include <cstdio>
#include <cstdlib>

//namespace FR = RTSeis::Utilities::FilterRepresentations;

class InstrumentResponse::InstrumentResponseImpl
{
public:
    class RTSeis::Utilities::FilterRepresentations mZPK;
    class RTSeis::Utilities::FilterRepresentations mBA;
    class RTSeis::Utilities::FilterRepresentations mSOS;
    double mSamplingRate = 0;
    bool mHaveResponse = false;
    bool mIsAnalog = true;
}

InstrumentResponse::InstrumentResponse :
    pImpl(std::make_unique<InstrumentResponse> ())
{
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
    const RTSeis::Utilities::FilterRepresentations &zpk)
{
    constexpr bool lanalog = true;
    setResponse(zpk, lanalog);
}

void InstrumentResponse::setDigitalResponse(
    const RTSeis::Utilities::FilterRepresentations &zpk)
{
    constexpr bool lanalog = false;
    setResponse(zpk, lanalog);
}

void InstrumentResponse::setResponse(
    const RTSeis::Utilities::FilterRepresentations::ZPK &zpk,
    const bool lanalog)
{
    pImpl->mIsAnalalog = lanalog;
}

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
