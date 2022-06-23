#include <stdexcept>
#include <string>
#include <cmath>
#include <ipps.h>
#include "rtseis/filterImplementations/detrend.hpp"

using namespace RTSeis::FilterImplementations;

template<class T>
class Detrend<T>::DetrendImpl
{
public:
//private:
    double mMean = 0;
    double mSlope = 0;
    double mIntercept = 0;
    DetrendType mType = DetrendType::LINEAR;
    bool mInitialized = false;
};

/// Constructor
template<class T>
Detrend<T>::Detrend() :
    pImpl(std::make_unique<DetrendImpl> ())
{
}

/// Copy constructor
template<class T>
Detrend<T>::Detrend(const Detrend &detrend)
{
    *this = detrend;
}

/// Move c'tor
template<class T>
Detrend<T>::Detrend(Detrend &&detrend) noexcept
{
    *this = std::move(detrend);
}

/// Copy constructor
template<class T>
Detrend<T>& Detrend<T>::operator=(const Detrend &detrend)
{
    if (&detrend == this){return *this;}
    pImpl = std::make_unique<DetrendImpl> (*detrend.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Detrend<T>& Detrend<T>::operator=(Detrend &&detrend) noexcept
{
    if (&detrend == this){return *this;}
    pImpl = std::move(detrend.pImpl);
    return *this;
}


/// Destructor
template<class T>
Detrend<T>::~Detrend() = default;

/// Clear/reset module
template<class T>
void Detrend<T>::clear() noexcept
{
    pImpl->mMean = 0;
    pImpl->mSlope = 0;
    pImpl->mIntercept = 0;
    pImpl->mType = DetrendType::LINEAR;
    pImpl->mInitialized = false;
}

/// Initialize
template<>
void Detrend<double>::initialize(const DetrendType type)
{
    clear();
    pImpl->mType = type;
    pImpl->mInitialized = true;
}

template<>
void Detrend<float>::initialize(const DetrendType type)
{
    clear();
    pImpl->mType = type;
    pImpl->mInitialized = true;
}

/// Cehck if inititalized
template<class T>
bool Detrend<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Apply
template<>
void Detrend<double>::apply(const int nx, const double x[], double *yin[])
{
    double *y = *yin;
    pImpl->mMean = 0;
    pImpl->mSlope = 0;
    pImpl->mIntercept = 0;
    if (nx <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    if (pImpl->mType == DetrendType::LINEAR)
    {
        removeTrend(nx, x, &y, &pImpl->mIntercept, &pImpl->mSlope);
    }
    else
    {
        removeMean(nx, x, &y, &pImpl->mMean);
    }
}

/// Apply
template<>
void Detrend<float>::apply(const int nx, const float x[], float *yin[])
{
    float *y = *yin;
    pImpl->mMean = 0;
    pImpl->mSlope = 0;
    pImpl->mIntercept = 0;
    if (nx <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    if (pImpl->mType == DetrendType::LINEAR)
    {
        float intercept;
        float slope;
        removeTrend(nx, x, &y, &intercept, &slope);
        pImpl->mIntercept = static_cast<double> (intercept);
        pImpl->mSlope = static_cast<float> (slope);
    }
    else
    {
        float mean;
        removeMean(nx, x, &y, &mean);
        pImpl->mMean = mean;
    }
}

/// Remove trend (float)
void RTSeis::FilterImplementations::removeTrend(
    const int length, const float x[], float *yin[],
    float *intercept, float *slope)
{
    if (length <= 0){return;}
    if (x == nullptr || *yin == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Handle an edge case - this is actually under-determined
    if (length < 2)
    {
        *yin[0] = 0.f;
        if (intercept != nullptr){*intercept = x[0];}
        if (slope != nullptr){*slope = 0.0f;}
        return;
    }
    // Mean of x - analytic formula for evenly spaced samples starting
    // at indx 0. This is computed by simplifying Gauss's formula.
    auto len64 = static_cast<uint64_t> (length);
    auto mean_x = 0.5*static_cast<double> (len64 - 1); 
    // Note, the numerator is the sum of consecutive squared numbers.
    // In addition we simplify.
    auto var_x = static_cast<double> ( ((len64 - 1))*(2*(len64 - 1) + 1) )/6.
               - mean_x*mean_x;
    float mean_y;
    ippsMean_32f(x, length, &mean_y, ippAlgHintAccurate);
    auto cov_xy = 0.0;
    #pragma omp simd reduction(+:cov_xy)
    for (auto i=0; i<length; ++i)
    {
        cov_xy = cov_xy
               + static_cast<double> (i)*static_cast<double> (x[i]);
    }
    // This is computed by expanding (x_i - bar(x))*(y_i - bar(y)),
    // using the definition of the mean, and simplifying
    cov_xy = (cov_xy/static_cast<double> (length)) - mean_x*mean_y;
    auto b1 = cov_xy/var_x;
    auto b0 = mean_y - b1*mean_x;
    auto b0f = static_cast<float> (b0);
    auto b1f = static_cast<float> (b1);
    // Remove the mean
    float *y = *yin;
    #pragma omp simd
    for (auto i=0; i<length; ++i)
    {   
        y[i] = x[i] - (b0f + b1f*static_cast<float> (i));
    }
    if (intercept){*intercept = b0;}
    if (slope){*slope = b1;}
}

/// Remove trend (double)
void RTSeis::FilterImplementations::removeTrend(
    const int length, const double x[], double *yin[],
    double *intercept, double *slope)
{
    if (length <= 0){return;}
    if (x == nullptr || *yin == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Handle an edge case - this is actually underdetermined
    if (length < 2)
    {
        *yin[0] = 0;
        if (intercept != nullptr){*intercept = x[0];}
        if (slope != nullptr){*slope = 0;}
        return;
    }
    // Mean of x - analytic formula for evenly spaced samples starting
    // at indx 0. This is computed by simplifying Gauss's formula.
    auto len64 = static_cast<uint64_t> (length);
    auto mean_x = 0.5*static_cast<double> (len64 - 1); 
    // Note, the numerator is the sum of consecutive squared numbers.
    // In addition we simplify.
    auto var_x = static_cast<double> ( ((len64 - 1))*(2*(len64 - 1) + 1) )/6.
               - mean_x*mean_x;
    double mean_y;
    ippsMean_64f(x, length, &mean_y);
    auto cov_xy = 0.0;
    #pragma omp simd reduction(+:cov_xy)
    for (auto i=0; i<length; ++i)
    {   
        cov_xy = cov_xy + static_cast<double> (i)*x[i];
    }   
    // This is computed by expanding (x_i - bar(x))*(y_i - bar(y)),
    // using the definition of the mean, and simplifying
    cov_xy = (cov_xy/static_cast<double> (length)) - mean_x*mean_y;
    auto b1 = cov_xy/var_x;
    auto b0  = mean_y - b1*mean_x;
    // Remove the mean
    double *y = *yin;
    #pragma omp simd
    for (auto i=0; i<length; ++i)
    {   
        y[i] = x[i] - (b0 + b1*static_cast<double> (i));
    }
    if (intercept){*intercept = b0;}
    if (slope){*slope = b1;}
}

/// Remove mean (double)
void RTSeis::FilterImplementations::removeMean(
    const int nx, const double x[], double *y[], double *mean)
{
    if (nx <= 0){return;} 
    if (x == nullptr || *y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    double pMean;
    ippsMean_64f(x, nx, &pMean); // Compute mean of input
    ippsSubC_64f(x, pMean, *y, nx); // y - mean(x)
    if (mean){*mean = pMean;}
}

/// Remove mean (float)
void RTSeis::FilterImplementations::removeMean(
    const int nx, const float x[], float *y[], float *mean)
{
    if (nx <= 0){return;} 
    if (x == nullptr || *y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    float pMean;
    ippsMean_32f(x, nx, &pMean, ippAlgHintAccurate);
    ippsSubC_32f(x, pMean, *y, nx); // y - mean(x)
    if (mean){*mean = pMean;}
}

///--------------------------------------------------------------------------///
///                        Template instantiation                            ///
///--------------------------------------------------------------------------///
template class RTSeis::FilterImplementations::Detrend<double>;
template class RTSeis::FilterImplementations::Detrend<float>;
