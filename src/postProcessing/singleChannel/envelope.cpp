#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include "rtseis/postProcessing/singleChannel/envelope.hpp"
#include "rtseis/utilities/windowFunctions.hpp"
#include "rtseis/private/throw.hpp"
#include "rtseis/enums.h"

using namespace RTSeis::PostProcessing::SingleChannel;

enum class Implementation
{
    FIR,     /*!< Use a Hilbert FIR filter to compute the envelope. */
    ANALYTIC /*!< Use the analytic signal as computed by the Fourier transform
                  to compute the envelope. */
};

class EnvelopeFIRParameters::EnvelopeFIRParametersImpl
{
public:
    /// The beta in the Kaiser window
    double beta = 8.0;
    /// The number of points in the FIR filter
    int nfir = 301;
    /// Precision 
    RTSeis::Precision precision = RTSeis::Precision::DOUBLE; 
    /// The processing mode
    const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING;
};

EnvelopeFIRParameters::EnvelopeFIRParameters(
    const int n,
    const double beta,
    const RTSeis::Precision precision) :
    pImpl(std::make_unique<EnvelopeFIRParametersImpl>())
{
    setFilterLength(n);
    setBeta(beta);
    setPrecision(precision);
}

EnvelopeFIRParameters::EnvelopeFIRParameters(
    const EnvelopeFIRParameters &envfir)
{
    *this = envfir;
}

EnvelopeFIRParameters::EnvelopeFIRParameters(
    EnvelopeFIRParameters &&envfir) noexcept
{
    *this = std::move(envfir);
}

EnvelopeFIRParameters& EnvelopeFIRParameters::operator=(
    const EnvelopeFIRParameters &envfir)
{
    if (&envfir == this){return *this;}
    pImpl = std::make_unique<EnvelopeFIRParametersImpl>();
    pImpl->nfir = envfir.pImpl->nfir;
    pImpl->beta = envfir.pImpl->beta;
    pImpl->precision = envfir.pImpl->precision;
    return *this;
}

EnvelopeFIRParameters& EnvelopeFIRParameters::operator=(
    EnvelopeFIRParameters &&envfir) noexcept
{
    if (&envfir == this){return *this;}
    pImpl = std::move(envfir.pImpl);
    return *this;
}

EnvelopeFIRParameters::~EnvelopeFIRParameters() = default;

void EnvelopeFIRParameters::clear() noexcept
{
    pImpl->nfir = 301;
    pImpl->beta = 8;
    pImpl->precision = RTSeis::Precision::DOUBLE;
}

void EnvelopeFIRParameters::setBeta(const double beta) noexcept
{
    pImpl->beta = beta;
}

double EnvelopeFIRParameters::getBeta() const noexcept
{
    return pImpl->beta;
}

void EnvelopeFIRParameters::setFilterLength(const int nfir)
{
    if (nfir < 1)
    {
        RTSEIS_THROW_IA("filter length = %d must be at least 1", nfir);
    }
    pImpl->nfir = nfir;
}

int EnvelopeFIRParameters::getFilterLength() const noexcept
{
    return pImpl->nfir;
}

void EnvelopeFIRParameters::setPrecision(
    const RTSeis::Precision precision) noexcept
{
    pImpl->precision = precision;
}

RTSeis::Precision EnvelopeFIRParameters::getPrecision() const noexcept
{
    return pImpl->precision;
}

bool EnvelopeFIRParameters::isValid() const noexcept
{
    if (pImpl->nfir < 1){return false;}
    return true;
}

//============================================================================//
//                          DFT Envelope Parameters                           //
//============================================================================//

class EnvelopeDFTParameters::EnvelopeDFTParametersImpl
{
public:
    /// Precision 
    RTSeis::Precision precision = RTSeis::Precision::DOUBLE;
};

EnvelopeDFTParameters::EnvelopeDFTParameters(
    const RTSeis::Precision precision) :
    pImpl(std::make_unique<EnvelopeDFTParametersImpl>())
{
    setPrecision(precision);
}

EnvelopeDFTParameters::EnvelopeDFTParameters(
    const EnvelopeDFTParameters &envdft)
{
    *this = envdft;
}

EnvelopeDFTParameters::EnvelopeDFTParameters(
    EnvelopeDFTParameters &&envdft) noexcept
{
    *this = std::move(envdft);
}

EnvelopeDFTParameters& EnvelopeDFTParameters::operator=(
    const EnvelopeDFTParameters &envdft)
{
    if (&envdft == this){return *this;}
    pImpl = std::make_unique<EnvelopeDFTParametersImpl>();
    pImpl->precision = envdft.pImpl->precision;
    return *this;
}

EnvelopeDFTParameters& EnvelopeDFTParameters::operator=(
    EnvelopeDFTParameters &&envdft) noexcept
{
    if (&envdft == this){return *this;}
    pImpl = std::move(envdft.pImpl);
    return *this;
}

EnvelopeDFTParameters::~EnvelopeDFTParameters() = default;

void EnvelopeDFTParameters::setPrecision(
    const RTSeis::Precision precision) noexcept
{
    pImpl->precision = precision;
    return;
}

RTSeis::Precision EnvelopeDFTParameters::getPrecision() const noexcept
{
    return pImpl->precision;
}

bool EnvelopeDFTParameters::isValid() const noexcept
{
    return true;
}

//============================================================================//
//                         Envelope Implementation                            //
//============================================================================//

class Envelope::EnvelopeImpl
{
public:
    EnvelopeFIRParameters envfir;
    EnvelopeDFTParameters envdft;
    Implementation implementation = Implementation::ANALYTIC;
    double mean = 0;
    bool linit = true; // Default DFT implementation will work
};

Envelope::Envelope() :
    pImpl(std::make_unique<EnvelopeImpl>())
{
}

Envelope::~Envelope() = default;

Envelope::Envelope(const EnvelopeFIRParameters &envfir) :
    pImpl(std::make_unique<EnvelopeImpl> ())
{
    if (!envfir.isValid())
    {
        pImpl->linit = false;
        RTSEIS_THROW_IA("%s", "Parameters are invalid");
    }
    pImpl->envfir = envfir;
    pImpl->linit = true;
    return;  
}

Envelope::Envelope(const EnvelopeDFTParameters &envdft) :
    pImpl(std::make_unique<EnvelopeImpl> ())
{
    if (!envdft.isValid())
    {
        pImpl->linit = false;
        RTSEIS_THROW_IA("%s", "Parameters are invalid");
    }
    pImpl->envdft = envdft;
    pImpl->linit = true;
    return;  
}
/*

void Envelope::apply(const int nx, const double yupper[], const double ylower[])
{
    if (nx <= 0){return;}
    if (yupper == nullptr){RTSEIS_THROW_IA("%s", "yupper is null");}
    if (ylower == nullptr){RTSEIS_THROW_IA("%s", "ylower is null");}
    apply(nx, x, y);
    // Compute the lower envelope from the upper envelope

    return; 
}

void Envelope::apply(const int nx, const double x[], double y[])
{
    if (nx <= 0){return;}
    if (!isInitialized())
    {
        RTSEIS_THROW_IA("%s", "Envelope not initialized");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_THROW_IA("%s", "y is NULL");}
        RTSEIS_THROW_IA("%s", "Invalid arrays");
    }
    // Remove the mean
    double pMean;
    ippsMean_64f(x, nx, &pMean);
    ippsSubC_64f(x, pMean, y); 
    // Compute the absolute value of the analytic signal (hilbert transform)
    if (pImpl->implementation == ANALYTIC)
    {
        // Compute the Hilbert transform
        Utilities::Transforms::Hilbert hilbert;
        hilbert(nx, RTSeis::Precision::DOUBLE);
        auto *yhilb = static_cast<std::complex<double> *> (IppsMalloc_64fc(nx));
        hilbert.apply(nx, y, yhilb); 
        // Take the absolute value of the analytic signal then restore mean
        auto *yhilbIPP = static_cast<Ipp64fc *> (yhilb);
        ippsPowerSpectr_64fc(yhilbIPP, y, nx); 
        // Now add in the mean
        ippsAddC_64f_I(y, pMean, nx);
        ippsFree(yhilb);
    }
    // Perform FIR filtering
    else
    {
        Utilities::FilterImplementations::FIRFilter firReal;
        Utilities::FilterImplementations::FIRFilter firImag;
    }
    return;
}
*/

bool Envelope::isInitialized(void) const noexcept
{
    return pImpl->linit;
}
