#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <memory>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/postProcessing/singleChannel/taper.hpp"
#include "rtseis/utilities/windowFunctions.hpp"
#include "rtseis/enums.h"

using namespace RTSeis::PostProcessing::SingleChannel;

class TaperParameters::TaperParametersImpl
{
public:
    /// Percentage of signal to apply taper to
    double pct = 5;
    /// Taper type.
    Type type = Type::HAMMING;
    /// Precision 
    RTSeis::Precision precision = RTSeis::Precision::DOUBLE; 
    /// The processing mode can't be toggled
    const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING;
};

class Taper::TaperImpl
{
public:
    TaperParameters parms; 
    std::vector<double> w8;
    std::vector<float>  w4;
    int winLen0 =-1;
    bool linit = true;
};

TaperParameters::TaperParameters(const double pct,
                                 const Type type,
                                 const RTSeis::Precision precision) :
    pImpl(std::make_unique<TaperParametersImpl>())
{
    setPercentage(pct);
    setTaperType(type);
    setPrecision(precision);
    return;
}

TaperParameters::TaperParameters(const TaperParameters &parms)
{
    *this = parms;
    return;
}

TaperParameters::TaperParameters(TaperParameters &&parms)
{
    *this = std::move(parms);
    return;
}

TaperParameters& TaperParameters::operator=(const TaperParameters &parms)
{
    if (&parms == this){return *this;}
    pImpl = std::make_unique<TaperParametersImpl> (*parms.pImpl);
    return *this; 
}

TaperParameters& TaperParameters::operator=(TaperParameters &&parms)
{
    if (&parms == this){return *this;}
    pImpl = std::move(parms.pImpl);
    return *this;
}

TaperParameters::~TaperParameters(void) = default;

void TaperParameters::setPrecision(const RTSeis::Precision precision)
{
    pImpl->precision = precision;
    return;
}

RTSeis::Precision TaperParameters::getPrecision(void) const
{
    return pImpl->precision;
}

void TaperParameters::setTaperType(const Type type)
{
    pImpl->type = type;
    return;
}

TaperParameters::Type TaperParameters::getTaperType(void) const
{
    return pImpl->type;
}

void TaperParameters::setPercentage(const double pct)
{
    if (pct < 0 || pct > 100)
    {
        RTSEIS_THROW_IA("%s", "Percentage must be in range [0,100]");
    }
    pImpl->pct = pct;
    return;
}

double TaperParameters::getPercentage(void) const
{
    return pImpl->pct;
}

void TaperParameters::clear(void)
{
    if (pImpl)
    {
        pImpl->pct = 5;
        pImpl->type = Type::HAMMING;
        pImpl->precision = RTSeis::Precision::DOUBLE;
    }
    return;
}

bool TaperParameters::isValid(void) const
{
    if (pImpl->pct < 0 || pImpl->pct > 100){return false;}
    return true;
}

//============================================================================//
//                                    Tapering                                //
//============================================================================//

Taper::Taper(void) :
    pImpl(std::make_unique<TaperImpl>())
{
    return;
}

Taper::Taper(const TaperParameters &parameters) :
    pImpl(std::make_unique<TaperImpl>())
{
    setParameters(parameters);
    return;
}

Taper::~Taper(void) = default;

void Taper::clear(void)
{
    pImpl->parms.clear();
    pImpl->w8.clear();
    pImpl->w4.clear();
    pImpl->winLen0 =-1;
    pImpl->linit = true;
    return;
}

Taper& Taper::operator=(const Taper &taper)
{
    if (&taper == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<TaperImpl> (*taper.pImpl);
    return *this;
}

void Taper::setParameters(const TaperParameters &parameters)
{
    clear(); // Sets winLen0 to -1
    if (!parameters.isValid())
    {
        RTSEIS_THROW_IA("%s", "Taper parameters are invalid");
    }
    pImpl->parms = parameters;
    return;
}

void Taper::apply(const int nx, const double x[], double y[])
{
    if (nx <= 0){return;}
    if (!pImpl->linit)
    {
        RTSEIS_THROW_IA("%s", "Taper never initialized");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_THROW_IA("%s", "y is NULL");}
        RTSEIS_THROW_IA("%s", "Invalid arguments");
    }
    // Deal with an edge case
    if (nx < 3)
    {
        y[0] = 0;
        if (nx == 2){y[1] = 0;}
        return;
    }
    // Compute taper length
    double pct = pImpl->parms.getPercentage();
    int npct = static_cast<int> (static_cast<double> (nx)*pct/100 + 0.5) + 1;
    int m = std::max(2, std::min(nx, npct));
    // Redesign the window?  If the parameters were (re)set then winLen0 is -1.
    // Otherwise, if the same length signal is coming at us then the precision
    // of the module can't change so we can just use the old window.
    if (pImpl->winLen0 != m)
    {
        pImpl->w8.resize(m); // Prevent reallocations in a bit
        TaperParameters::Type type = pImpl->parms.getTaperType();
        double *w8data = pImpl->w8.data();
        if (type == TaperParameters::Type::HAMMING)
        {
            RTSeis::Utilities::WindowFunctions::hamming(m, &w8data);
        }
        else if (type == TaperParameters::Type::BLACKMAN)
        {
            RTSeis::Utilities::WindowFunctions::blackman(m, &w8data);
        }
        else if (type == TaperParameters::Type::HANN)
        {
            RTSeis::Utilities::WindowFunctions::hann(m, &w8data);
        }
        else if (type == TaperParameters::Type::BARTLETT)
        {
            RTSeis::Utilities::WindowFunctions::bartlett(m, &w8data);
        }
        else if (type == TaperParameters::SINE)
        {
            RTSeis::Utilities::WindowFunctions::sine(m, &w8data);
        }
        else
        {
#ifdef DEBUG
            assert(false);
#endif
            RTSEIS_THROW_IA("%s", "Unsupported window");
        }
    }
    // Taper first (m+1)/2 points
    int mp12 = m/2;
    const double *w = pImpl->w8.data();
    ippsMul_64f(w, x, y, mp12);
    // Copy the intermediate portion of the signal
    int ncopy = nx - mp12 - mp12; // Subtract out two window lengths
    if (ncopy > 0)
    {
        ippsCopy_64f(&x[mp12], &y[mp12], ncopy);
    } 
    // Taper last (m+1)/2 points 
    ippsMul_64f(&w[m-mp12], &x[nx-mp12], &y[nx-mp12], mp12);
    return;
}

bool Taper::isInitialized(void) const
{
    return pImpl->linit;
}
