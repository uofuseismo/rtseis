#include <iostream>
#include <complex>
#include <cmath>
#include <cassert>
#include "rtseis/utilities/transforms/wavelets/morlet.hpp"

using namespace RTSeis::Utilities::Transforms::Wavelets;

namespace
{

/*
template<typename T>
void evaluateFrequencyDomain(const int n,
                             const T waveNumber,
                             const T dt,
                             const T scale,
                             const T k[],
                             std::complex<T> *daughter)
{
    const std::complex<T> czero(0, 0);
    const T half = 0.5;
    const T zero = 0;
    const T pi14 = 0.7511255444649425; // std::pow(M_PI, -0.25);
    T arg = 2*M_PI*scale/dt;
    T norm = std::sqrt(arg)*pi14;
    int nhalf = n/2;
    std::fill(daughter, daughter + n, czero);
    #pragma omp simd
    for (int i=0; i<nhalf + 1; ++i)
    {
        T res = scale*k[i] - waveNumber;
        T expnt = -half*(res*res);
        daughter[i] = std::complex<T> (norm*expnt, zero);
    }
}
*/

template<typename T>
void evaluateTimeDomain(const int n, 
                        const T omega0, // Dimensionless wavenumber: omega_0
                        const T scale,  // Frequency (Hz)
                        std::complex<T> *daughter,
                        const bool normalize)
{
    const T one = 1;
    const T half = 0.5;
    //auto norm = static_cast<T> (1./(std::sqrt(scale*std::sqrt(M_PI))));
    T norm = 0.7511255444649425;
    if (normalize){norm = one/std::sqrt(std::abs(scale))*norm;}
    T xhalf = half*static_cast<T> (n - 1);
    #pragma omp simd
    for (int i=0; i<n; ++i)
    {
        T xs = (static_cast<T> (i) - xhalf)/scale;
        T damp = norm*std::exp(-half*(xs*xs));
        T carg = omega0*xs;
        // \pi^{-1/4} \exp(i \omega_0 \eta) 
        auto dr = damp*std::cos(carg);
        auto di = damp*std::sin(carg);
        daughter[i] = std::complex<T> (dr, di);
    }
}

}

class Morlet::MorletImpl
{
public:
    //double mSamplingPeriod = 1; // Sampling period (seconds)
    double mOmega0 = 6; // Morlet wavelet parameter
    bool mNormalize = false; // Divide by 1/sqrt(s)
};

/// C'tor
Morlet::Morlet() :
    pImpl(std::make_unique<MorletImpl> ())
{
}

/// Copy c'tor
Morlet::Morlet(const Morlet &morlet)
{
    *this = morlet;
}

/// Move c'tor
[[maybe_unused]]
Morlet::Morlet(Morlet &&morlet) noexcept
{
    *this = std::move(morlet);
}

/// Clone
std::unique_ptr<IContinuousWavelet>
Morlet::clone() const
{
    return std::make_unique<Morlet> (*this);
}

/// Copy assignment
Morlet& Morlet::operator=(const Morlet &morlet)
{
    if (&morlet == this){return *this;}
    pImpl = std::make_unique<MorletImpl> (*morlet.pImpl);
    return *this;
}

/// Move assignment
Morlet& Morlet::operator=(Morlet &&morlet) noexcept
{
    if (&morlet == this){return *this;}
    pImpl = std::move(morlet.pImpl);
    return *this;
}

/// Destructor
Morlet::~Morlet() = default;

/// Wavelet parameter 
void Morlet::setParameter(double omega0)
{
    if (omega0 <= 0)
    {
        throw std::invalid_argument("omega0 = " + std::to_string(omega0)
                                  + " must be positive");
    }
    pImpl->mOmega0 = omega0;
}

double Morlet::getParameter() const noexcept
{
    return pImpl->mOmega0;
}

/// Normalization
void Morlet::enableNormalization() noexcept
{
    pImpl->mNormalize = true;
}

void Morlet::disableNormalization() noexcept
{
    pImpl->mNormalize = false;
}

bool Morlet::normalize() const noexcept
{
    return pImpl->mNormalize;
}

/// Sampling period
/*
void Morlet::setSamplingPeriod(const double dt)
{
    if (dt <= 0)
    {
        throw std::invalid_argument("dt = " + std::to_string(dt)
                                  + " must be positive");
    }
    pImpl->mSamplingPeriod = dt;
}

double Morlet::getSamplingPeriod() const noexcept
{
    return pImpl->mSamplingPeriod;
}
*/

/*
/// Evaluate
void Morlet::evaluate(const int n, const double scale,
                      const double k[],
                      std::complex<double> *daughterIn[]) const
{
    if (n < 1){return;}
    auto daughter = *daughterIn;
    if (scale < 0){throw std::invalid_argument("The scale must be positive");}
    if (k == nullptr){throw std::invalid_argument("k is NULL");}
    if (daughter == nullptr){throw std::invalid_argument("daughter is NULL");}
    evaluateFrequencyDomain(n,
                   getWaveNumber(),
                   getSamplingPeriod(),
                   scale,
                   k, daughter);
}

void Morlet::evaluate(const int n, const float scale,
                      const float k[],
                      std::complex<float> *daughterIn[]) const
{
    if (n < 1){return;}
    auto daughter = *daughterIn;
    if (scale < 0){throw std::invalid_argument("The scale must be positive");}
    if (k == nullptr){throw std::invalid_argument("k is NULL");}
    if (daughter == nullptr){throw std::invalid_argument("daughter is NULL");}
    evaluateFrequencyDomain(n,
                            static_cast<float> (getWaveNumber()),
                            static_cast<float> (getSamplingPeriod()),
                            scale,
                            k, daughter);
}

/// Cone of influence
double Morlet::computeConeOfInfluenceScalar() const noexcept
{
    const double sqrt2i = 1./std::sqrt(2.);
    const double fourpi = 4*M_PI; 
    auto param = getWaveNumber();
    auto fourierFactor = fourpi/(param + std::sqrt(2 + param*param));
    auto coiScalar = fourierFactor*sqrt2i;
    return coiScalar;
}
*/

void Morlet::evaluate(const int n, const double scale,
                      std::complex<double> *daughterIn[]) const
{
    if (n < 1){return;}
    if (scale <= 0){throw std::invalid_argument("scale must be positive");}
    auto daughter = *daughterIn;
    if (daughter == nullptr){throw std::invalid_argument("daughter is NULL");}
    auto omega0 = getParameter();
    auto lnorm = normalize();
    evaluateTimeDomain(n, omega0, scale, daughter, lnorm);
}

void Morlet::evaluate(const int n, const float scale,
                      std::complex<float> *daughterIn[]) const
{
    if (n < 1){return;}
    if (scale <= 0){throw std::invalid_argument("scale must be positive");}
    auto daughter = *daughterIn;
    if (daughter == nullptr){throw std::invalid_argument("daughter is NULL");}
    auto omega0 = static_cast<float> (getParameter());
    auto lnorm = normalize();
    evaluateTimeDomain(n, omega0, scale, daughter, lnorm);
}

/// Clear
void Morlet::clear() noexcept
{
    pImpl->mOmega0 = 6;
    pImpl->mNormalize = false;
} 
