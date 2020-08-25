#include <cstdio>
#include <cstdlib>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/normalization/minMax.hpp"


using namespace RTSeis::Utilities::Normalization;

class MinMax::MinMaxImpl
{
public:
    double mDataMin64f = 0;
    double mDataMax64f = 0;
    double mTargetMin64f = 0;
    double mTargetMax64f = 0;
    double mDen64f = 0;
    bool mRescaleToUnitInterval = true; 
    bool mInitialized = false;
};

/// Constructors
MinMax::MinMax() :
    pImpl(std::make_unique<MinMaxImpl> ())
{
}

[[maybe_unused]] MinMax::MinMax(const MinMax &minMax)
{
    *this = minMax;
}

MinMax::MinMax(MinMax &&minMax) noexcept
{
    *this = std::move(minMax);
}

/// Operators
MinMax& MinMax::operator=(const MinMax &minMax)
{
    if (&minMax == this){return *this;}
    pImpl = std::make_unique<MinMaxImpl> (*minMax.pImpl);
    return *this;
}

MinMax& MinMax::operator=(MinMax &&minMax) noexcept
{
    if (&minMax == this){return *this;}
    pImpl = std::move(minMax.pImpl);
    return *this;
}

/// Destructors
MinMax::~MinMax() = default;

void MinMax::clear() noexcept
{
    pImpl->mDataMin64f = 0;
    pImpl->mDataMax64f = 0;
    pImpl->mTargetMin64f = 0;
    pImpl->mTargetMax64f = 0;
    pImpl->mDen64f = 0;
    pImpl->mRescaleToUnitInterval = true;
    pImpl->mInitialized = false;
}

/// Initializers
void MinMax::initialize(const std::pair<double, double> dataMinMax,
                        const std::pair<double, double> targetMinMax)
{
    clear();
    if (dataMinMax.first == dataMinMax.second)
    {
        RTSEIS_THROW_IA("dataMin = %lf cannot equal dataMax == %lf", 
                        dataMinMax.first, dataMinMax.second);
    }
    pImpl->mDataMin64f = dataMinMax.first;
    pImpl->mDataMax64f = dataMinMax.second;
    pImpl->mTargetMin64f = targetMinMax.first;
    pImpl->mTargetMax64f = targetMinMax.second;;
    pImpl->mDen64f = (pImpl->mDataMax64f - pImpl->mDataMin64f);
    pImpl->mRescaleToUnitInterval = true;
    if (pImpl->mTargetMin64f != 0 || pImpl->mTargetMax64f != 1)
    {
        pImpl->mRescaleToUnitInterval = false;
    }
    pImpl->mInitialized = true;
}

void MinMax::initialize(const int npts, const double x[],
                        const std::pair<double, double> targetMinMax)
{
    clear();
    if (npts < 2 || x == nullptr)
    {
        if (npts < 2){RTSEIS_THROW_IA("npts = %d must be at least 2", npts);}
        RTSEIS_THROW_IA("%s", "x cannot be NULL");
    }
    // Find the min and max
    double dataMin, dataMax;
    ippsMinMax_64f(x, npts, &dataMin, &dataMax);
    if (dataMin == dataMax)
    {
        RTSEIS_THROW_IA("%s", "All elements of x are identical");
    }
    std::pair<double, double> dataMinMax(dataMin, dataMax);
    initialize(dataMinMax, targetMinMax);
}        

/// Apply min max normalization
void MinMax::apply(const int npts, const double x[], double *yIn[])
{
    if (npts < 1){return;}
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    double *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    // Rescale
    if (pImpl->mRescaleToUnitInterval)
    {
        // y = 0 + (x - xmin)/(xmax - xmin)*(1 - 0)
        ippsNormalize_64f(x, y, npts, pImpl->mDataMin64f, pImpl->mDen64f);
    }
    else
    {
        // y = a + (x - xmin)/(xmax - xmin)*(b - a)
        //   = a - xmin/(xmax - xmin)*(b - a) + x*(b - a)/(xmax - xmin)
        //   = \alpha + \beta x 
        auto bma = pImpl->mTargetMax64f - pImpl->mTargetMin64f;
        auto alpha = pImpl->mTargetMin64f
                   - (pImpl->mDataMin64f*bma)/pImpl->mDen64f;
        auto beta  = bma/pImpl->mDen64f;
        #pragma omp simd
        for (auto i=0; i<npts; ++i)
        {
            y[i] = alpha + beta*x[i];
        }
    }
}

void MinMax::apply(const int npts, const float *x, float *yIn[])
{
    if (npts < 1){return;}
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    float *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    // Rescale
    if (pImpl->mRescaleToUnitInterval)
    {
        // y = 0 + (x - xmin)/(xmax - xmin)*(1 - 0)
        ippsNormalize_32f(x, y, npts,
                          static_cast<float> (pImpl->mDataMin64f),
                          static_cast<float> (pImpl->mDen64f));
    }
    else
    {
        auto bma = pImpl->mTargetMax64f - pImpl->mTargetMin64f;
        auto alpha = static_cast<float> (
            pImpl->mTargetMin64f
                - (pImpl->mDataMin64f * bma) / pImpl->mDen64f);
        auto beta = static_cast<float> (bma / pImpl->mDen64f);
        #pragma omp simd
        for (auto i = 0; i < npts; ++i)
        {
            y[i] = alpha + beta * x[i];
        }
    }
}

bool MinMax::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}