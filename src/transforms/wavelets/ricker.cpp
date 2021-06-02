#include <iostream>
#include <complex>
#include <cmath>
#include <cassert>
#include "rtseis/utilities/transforms/wavelets/ricker.hpp"
#include "rtseis/utilities/transforms/wavelets/derivativeOfGaussian.hpp"

using namespace RTSeis::Utilities::Transforms::Wavelets;

class Ricker::RickerImpl
{
public:
    DerivativeOfGaussian mDOG; // Will specialize to second derivative
};

/// C'tor
Ricker::Ricker() :
    pImpl(std::make_unique<RickerImpl> ())
{
    pImpl->mDOG.setOrder(2);
}

/// Copy c'tor
Ricker::Ricker(const Ricker &ricker)
{
    *this = ricker;
}

/// Move c'tor
[[maybe_unused]]
Ricker::Ricker(Ricker &&ricker) noexcept
{
    *this = std::move(ricker);
}

/// Clone
std::unique_ptr<IContinuousWavelet>
Ricker::clone() const
{
    return std::make_unique<Ricker> (*this);
}

/// Copy assignment
Ricker& Ricker::operator=(const Ricker &ricker)
{
    if (&ricker == this){return *this;}
    pImpl = std::make_unique<RickerImpl> (*ricker.pImpl);
    return *this;
}

/// Move assignment
Ricker& Ricker::operator=(Ricker &&ricker) noexcept
{
    if (&ricker == this){return *this;}
    pImpl = std::move(ricker.pImpl);
    return *this;
}

/// Destructor
Ricker::~Ricker() = default;

/// Sampling period
void Ricker::setSamplingPeriod(const double dt)
{
    pImpl->mDOG.setSamplingPeriod(dt);
}

double Ricker::getSamplingPeriod() const noexcept
{
    return pImpl->mDOG.getSamplingPeriod();
}

/// Evaluate
void Ricker::evaluate(const int n, const double scale,
                      const double k[], std::complex<double> *w[]) const
{
    pImpl->mDOG.evaluate(n, scale, k, w);
}

void Ricker::evaluate(const int n, const float scale,
                      const float k[], std::complex<float> *w[]) const
{
    pImpl->mDOG.evaluate(n, scale, k, w);
}

/// Clear
void Ricker::clear() noexcept
{
    pImpl->mDOG.clear();
    pImpl->mDOG.setOrder(2);
} 

double Ricker::computeConeOfInfluenceScalar() const noexcept
{
    return pImpl->mDOG.computeConeOfInfluenceScalar();
}
