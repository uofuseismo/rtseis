#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ipps.h>
#include <ippversion.h>
#include <ippcore.h>
#define RTSEIS_LOGGING 1
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/filterImplementations/detrend.hpp"

using namespace RTSeis::Utilities::FilterImplementations;

class Detrend::DetrendImpl
{
public:
//private:
    double mMean = 0;
    double mSlope = 0;
    double mIntercept = 0;
    RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
    DetrendType mType=  DetrendType::LINEAR;
    bool mInitialized = false;
};

/// Constructor
Detrend::Detrend() :
    pImpl(std::make_unique<DetrendImpl> ())
{
}

/// Copy constructor
Detrend::Detrend(const Detrend &detrend)
{
    *this = detrend;
}

/*
Detrend::Detrend(const Detrend &&detrend)
{
    *this = std::move(detrend);
}
*/

/// Copy constructor
Detrend& Detrend::operator=(const Detrend &detrend)
{
    if (&detrend == this){return *this;}
    pImpl = std::make_unique<DetrendImpl> (*detrend.pImpl);
    return *this;
}

/*
Detrend& Detrend::operator=(const Detrend &&detrend)
{
    if (&detrend == this){return *this;}
    pImpl = std::move(detrend.pImpl);
}
*/

/// Destructor
Detrend::~Detrend() = default;

/// Clear/reset module
void Detrend::clear() noexcept
{
    pImpl->mMean = 0;
    pImpl->mSlope = 0;
    pImpl->mIntercept = 0;
    pImpl->mPrecision = RTSeis::Precision::DOUBLE;
    pImpl->mType = DetrendType::LINEAR;
    pImpl->mInitialized = false;
}

/// Initialize
void Detrend::initialize(const DetrendType type,
                         const RTSeis::Precision precision)
{
    clear();
    pImpl->mType = type;
    pImpl->mPrecision = precision;
    pImpl->mInitialized = true;
}

/// Cehck if inititalized
bool Detrend::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Apply
void Detrend::apply(const int nx, const double x[], double *yin[])
{
    double *y = *yin;
    pImpl->mMean = 0;
    pImpl->mSlope = 0;
    pImpl->mIntercept = 0;
    if (nx <= 0){return;}
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not inititialized");
    }
    if (pImpl->mType == DetrendType::LINEAR && nx < 2)
    {   
        RTSEIS_THROW_IA("%s", "At least 2 data points for linear detrend");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    if (pImpl->mPrecision == RTSeis::Precision::FLOAT)
    {
        float *x32 = ippsMalloc_32f(nx);
        ippsConvert_64f32f(x, x32, nx);
        float *y32 = ippsMalloc_32f(nx);
        apply(nx, x32, &y32);
        ippsConvert_32f64f(y32, *yin, nx); 
        ippsFree(x32);
        ippsFree(y32);
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
void Detrend::apply(const int nx, const float x[], float *yin[])
{
    float *y = *yin;
    pImpl->mMean = 0;
    pImpl->mSlope = 0;
    pImpl->mIntercept = 0;
    if (nx <= 0){return;}
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not inititialized");
    }
    if (pImpl->mType == DetrendType::LINEAR && nx < 2)
    {
        RTSEIS_THROW_IA("%s", "At least 2 data points for linear detrend");
    } 
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    if (pImpl->mPrecision == RTSeis::Precision::DOUBLE)
    {
        double *x64 = ippsMalloc_64f(nx);
        ippsConvert_32f64f(x, x64, nx);
        double *y64 = ippsMalloc_64f(nx);
        apply(nx, x64, &y64);
        ippsConvert_64f32f(y64, *yin, nx); 
        ippsFree(x64);
        ippsFree(y64);
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

void RTSeis::Utilities::FilterImplementations::removeTrend(
    const int length, const float x[], float *yin[],
    float *intercept, float *slope)
{
    if (length <= 0){return;}
    if (length < 2 || x == nullptr || *yin == nullptr)
    {
        if (length < 2)
        {
            RTSEIS_THROW_IA("length = %d must be at least 2", length);
        }
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
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

void RTSeis::Utilities::FilterImplementations::removeTrend(
    const int length, const double x[], double *yin[],
    double *intercept, double *slope)
{
    if (length <= 0){return;}
    if (length < 2 || x == nullptr || *yin == nullptr)
    {
        if (length < 2)
        {
            RTSEIS_THROW_IA("length = %d must be at least 2", length);
        }
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
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

void RTSeis::Utilities::FilterImplementations::removeMean(
    const int nx, const double x[], double *y[], double *mean)
{
    if (nx <= 0){return;} 
    if (x == nullptr || *y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    double pMean;
    ippsMean_64f(x, nx, &pMean); // Compute mean of input
    ippsSubC_64f(x, pMean, *y, nx); // y - mean(x)
    if (mean){*mean = pMean;}
}
void RTSeis::Utilities::FilterImplementations::removeMean(
    const int nx, const float x[], float *y[], float *mean)
{
    if (nx <= 0){return;} 
    if (x == nullptr || *y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    float pMean;
    ippsMean_32f(x, nx, &pMean, ippAlgHintAccurate);
    ippsSubC_32f(x, pMean, *y, nx); // y - mean(x)
    if (mean){*mean = pMean;}
}