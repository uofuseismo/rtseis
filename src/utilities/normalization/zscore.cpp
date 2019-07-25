#include <cstdio>
#include <cstdlib>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/normalization/zscore.hpp"

using namespace RTSeis::Utilities::Normalization;

class ZScore::ZScoreImpl
{
public:
    double mMean64f = 0;
    double mStd64f = 0;
    float mMean32f = 0;
    float mStd32f = 0;
    bool mInitialized = false;
};

/// Constructors
ZScore::ZScore() :
    pImpl(std::make_unique<ZScoreImpl> ())
{
}

ZScore::ZScore(const ZScore &zscore)
{
    *this = zscore;
}

ZScore::ZScore(ZScore &&zscore) noexcept
{
    *this = std::move(zscore);
}

/// Operators
ZScore& ZScore::operator=(const ZScore &zscore)
{
    if (&zscore == this){return *this;}
    pImpl = std::make_unique<ZScoreImpl> (*zscore.pImpl);
    return *this;
}

ZScore& ZScore::operator=(ZScore &&zscore) noexcept
{
    if (&zscore == this){return *this;}
    pImpl = std::move(zscore.pImpl);
    return *this;
}

/// Destructors
ZScore::~ZScore() = default;

void ZScore::clear() noexcept
{
    pImpl->mMean64f = 0;
    pImpl->mStd64f = 0;
    pImpl->mMean32f = 0;
    pImpl->mStd32f = 0;
    pImpl->mInitialized = false;
}

/// Initializes the class
void ZScore::initialize(const double mean,
                        const double std)
{
    clear();
    if (std <= 0)
    {
        RTSEIS_THROW_IA("standard deviation = %lf must be positive", std);
    }
    pImpl->mMean64f = mean;
    pImpl->mStd64f = std;
    pImpl->mMean32f = static_cast<float> (mean);
    pImpl->mStd32f = static_cast<float> (mean);
    pImpl->mInitialized = true;
}

void ZScore::initialize(const int nx, const double x[])
{
    clear();
    if (nx < 2 || x == nullptr)
    {
        if (nx < 2){RTSEIS_THROW_IA("nx = %d must be at least 2", nx);}
        RTSEIS_THROW_IA("%s", "x is NULL");
    }
    double pMean, pStdDev;
    ippsMeanStdDev_64f(x, nx, &pMean, &pStdDev);
    if (pStdDev == 0)
    {
        RTSEIS_THROW_IA("%s", "x cannot be all the same values");
    }
    initialize(pMean, pStdDev); 
}

bool ZScore::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Standardizes the data
void ZScore::apply(const int npts, const double x[], double *yIn[])
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
    ippsNormalize_64f(x, y, npts, pImpl->mMean64f, pImpl->mStd64f);
}

void ZScore::apply(const int npts, const float x[], float *yIn[])
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
    ippsNormalize_32f(x, y, npts, pImpl->mMean32f, pImpl->mStd32f);
}
