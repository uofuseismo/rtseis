#include <iostream>
#include <complex>
#include <cmath>
#include <cassert>
#include "rtseis/utilities/transforms/wavelets/derivativeOfGaussian.hpp"

using namespace RTSeis::Utilities::Transforms::Wavelets;

namespace
{
template <typename T>
void evaluateDOG(const int n,
                 const bool isReal,
                 const int m,
                 const int derivativeSign,
                 const T dt,
                 const T gammaNorm,
                 const T scale,
                 const T k[],
                 std::complex<T> *daughter)
{
    const T half = 0.5;
    T arg = 2*M_PI*scale/dt;
    T norm =-derivativeSign*std::sqrt(arg)*gammaNorm;
    // Evaluate ala Torrence and Compo 1998, Table 1
    if (isReal)
    {
        #pragma omp simd
        for (int i=0; i<n; ++i)
        {
            T sk = scale*k[i];
            T halfsk2 = half*(sk*sk);
            T re = norm*std::pow(sk, m)*std::exp(-halfsk2);
            daughter[i] = std::complex<T> (re, 0);
        }
    }
    else
    {
        #pragma omp simd
        for (int i=0; i<n; ++i)
        {
            T sk = scale*k[i];
            T halfsk2 = half*(sk*sk);
            T im = norm*std::pow(sk, m)*std::exp(-halfsk2);
            daughter[i] = std::complex<T> (0, im);
        }
    }
}
}

class DerivativeOfGaussian::DerivativeOfGaussianImpl
{
public:
    double mNorm = 0; // Normalization $ \frac{1}{\sqrt{ \Gamma(m + 1/2) }}
    double mSamplingPeriod = 1; // Sampling period (seconds)
    int mDerivativeSign =-1; // Sign on derivative (2nd derivative)
    int mOrder = 2; // Ricker (Same as Torrence's software)
    bool mIsReal = true; // Result is real or complex (2nd derivative real)
};

/// C'tor
DerivativeOfGaussian::DerivativeOfGaussian() :
    pImpl(std::make_unique<DerivativeOfGaussianImpl> ())
{
}

/// Copy c'tor
DerivativeOfGaussian::DerivativeOfGaussian(const DerivativeOfGaussian &dog)
{
    *this = dog;
}

/// Move c'tor
DerivativeOfGaussian::DerivativeOfGaussian(DerivativeOfGaussian &&dog) noexcept
{
    *this = std::move(dog);
}

/// Copy assignment
DerivativeOfGaussian& 
DerivativeOfGaussian::operator=(const DerivativeOfGaussian &dog)
{
    if (&dog == this){return *this;}
    pImpl = std::make_unique<DerivativeOfGaussianImpl> (*dog.pImpl);
    return *this;
}

/// Move assignment
DerivativeOfGaussian&
DerivativeOfGaussian::operator=(DerivativeOfGaussian &&dog) noexcept
{
    if (&dog == this){return *this;}
    pImpl = std::move(dog.pImpl);
    return *this;
}

/// Clone
std::unique_ptr<IContinuousWavelet>
DerivativeOfGaussian::clone() const
{
    return std::make_unique<DerivativeOfGaussian> (*this);
}

/// Destructor
DerivativeOfGaussian::~DerivativeOfGaussian() = default;

/// Clear
void DerivativeOfGaussian::clear() noexcept
{
    pImpl->mNorm = 1;
    pImpl->mSamplingPeriod = 1;
    pImpl->mDerivativeSign =-1;
    pImpl->mOrder = 2;
    pImpl->mIsReal = true;
}

/// Set the order
void DerivativeOfGaussian::setOrder(int order)
{
    if (order < 0)
    {
        throw std::invalid_argument("Order = " + std::to_string(order)
                                  + " must be positive");
    }
    pImpl->mOrder = order;
    pImpl->mNorm = 1/std::sqrt( std::tgamma(order + 0.5) );
    // Derivative in Fourier domain is (i)^m
    if (order%4 == 0) // i**0 = 1
    {
        pImpl->mDerivativeSign = 1;
        pImpl->mIsReal = true;
    }
    else if (order%4 == 1) // i**1 = i
    {
        pImpl->mDerivativeSign = 1;
        pImpl->mIsReal = false;
    }
    else if (order%4 == 2) // i**2 =-1
    {
        pImpl->mDerivativeSign =-1;
        pImpl->mIsReal = true;
    }
    else //if (order%4 == 3) // i**3 =-i
    {
#ifndef DNDEBUG
        assert(order%4 == 3);
#endif
        pImpl->mDerivativeSign =-1;
        pImpl->mIsReal = false;
    }
}

int DerivativeOfGaussian::getOrder() const noexcept
{
    return pImpl->mOrder;
}

/// Set the sampling period
void DerivativeOfGaussian::setSamplingPeriod(const double dt)
{
    if (dt <= 0)
    {
        throw std::invalid_argument("dt = " + std::to_string(dt)
                                  + " must be positive");
    }
    pImpl->mSamplingPeriod = dt;
}

double DerivativeOfGaussian::getSamplingPeriod() const noexcept
{
    return pImpl->mSamplingPeriod;
}

/*
/// Initialize the class
void DerivativeOfGaussian::initialize(const int order, const double dt)
{
    clear();
    if (order < 0){throw std::invalid_argument("Order must be positive");}
    if (dt <= 0){throw std::invalid_argument("dt must be positive");}
    // Normalization factor taken right from Table 1 of Torrence and Compo, 1998
    pImpl->mNorm = 1/std::sqrt( std::tgamma(order + 0.5) );
    // Derivative in Fourier domain is (i)^m 
    if (order%4 == 0) // i**0 = 1
    {
        pImpl->mDerivativeSign = 1;
        pImpl->mIsReal = true;
    }
    else if (order%4 == 1) // i**1 = i
    {
        pImpl->mDerivativeSign = 1;
        pImpl->mIsReal = false;
    }
    else if (order%4 == 2) // i**2 =-1
    {
        pImpl->mDerivativeSign =-1;
        pImpl->mIsReal = true;
    }
    else //if (order%4 == 3) // i**3 =-i
    {
#ifndef DNDEBUG
        assert(order%4 == 3);
#endif
        pImpl->mDerivativeSign =-1;
        pImpl->mIsReal = false;
    } 
    pImpl->mSamplingPeriod = dt;
    pImpl->mOrder = order;
} 

/// Initialized?
bool DerivativeOfGaussian::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}
*/

/*
/// Get the e-folding time
double DerivativeOfGaussian::getEFoldingTime(double s) const
{
    return M_SQRT2*s;
}

/// Get the wavelength
double DerivativeOfGaussian::getWavelength(double s) const
{
    double lambda = (2*M_PI*s)/std::sqrt(pImpl->mOrder + 0.5);
    return lambda; 
}
*/

/// Evaluate the wavelet
void DerivativeOfGaussian::evaluate(const int n,
                                    const double scale,
                                    const double k[],
                                    std::complex<double> *daughterIn[]) const
{
    if (n < 1){return;}
    auto daughter = *daughterIn;
    if (scale < 0){throw std::invalid_argument("The scale must be positive");}
    if (k == nullptr){throw std::invalid_argument("k is NULL");}
    if (daughter == nullptr){throw std::invalid_argument("daughter is NULL");}
    evaluateDOG(n,
                pImpl->mIsReal, getOrder(), pImpl->mDerivativeSign,
                getSamplingPeriod(), pImpl->mNorm, scale, k, daughter); 
/*
    auto m = getOrder();
    auto dt = getSamplingPeriod();
    double norm =-pImpl->mDerivativeSign*std::sqrt(2*M_PI*scale/dt)*pImpl->mNorm;
    // Evaluate ala Torrence and Compo 1998, Table 1
    if (pImpl->mIsReal)
    {
        #pragma omp simd
        for (int i=0; i<n; ++i)
        {
            auto sk = scale*k[i];
            auto halfsk2 = 0.5*(sk*sk);
            auto re = norm*std::pow(sk, m)*std::exp(-halfsk2); 
            daughter[i] = std::complex<double> (re, 0);
        }
    } 
    else
    {
        #pragma omp simd
        for (int i=0; i<n; ++i)
        {
            auto sk = scale*k[i];
            auto halfsk2 = 0.5*(sk*sk);
            auto im = norm*std::pow(sk, m)*std::exp(-halfsk2);
            daughter[i] = std::complex<double> (0, im);
        }
    }
*/
}

void DerivativeOfGaussian::evaluate(const int n,
                                    const float scale,
                                    const float k[],
                                    std::complex<float> *daughterIn[]) const
{
    if (n < 1){return;}
    auto daughter = *daughterIn;
    if (scale < 0){throw std::invalid_argument("The scale must be positive");}
    if (k == nullptr){throw std::invalid_argument("k is NULL");}
    if (daughter == nullptr){throw std::invalid_argument("daughter is NULL");}
    evaluateDOG(n,
                pImpl->mIsReal,
                getOrder(),
                pImpl->mDerivativeSign,
                static_cast<float> (getSamplingPeriod()),
                static_cast<float> (pImpl->mNorm),
                scale, k, daughter);
}

double DerivativeOfGaussian::computeConeOfInfluenceScalar() const noexcept
{
    const double sqrt2i = 1./std::sqrt(2.);
    auto m = getOrder();
    auto fourierFactor = 2*M_PI/std::sqrt(m + 0.5);
    auto coiScalar = fourierFactor*sqrt2i;
    return coiScalar;
}
