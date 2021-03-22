#include <cstdio>
#include <cstdlib>
#include <string>
#include <ipps.h>
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

/// Copy c'tor
ZScore::ZScore(const ZScore &zscore)
{
    *this = zscore;
}

/// Move c'tor
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

/// Clears the class
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
        throw std::invalid_argument("standard deviation = "
                                  + std::to_string(std) + " must be positive");
    }
    pImpl->mMean64f = mean;
    pImpl->mStd64f = std;
    pImpl->mMean32f = static_cast<float> (mean);
    pImpl->mStd32f = static_cast<float> (mean);
    pImpl->mInitialized = true;
}

/// Initialize the class
void ZScore::initialize(const int nx, const double x[])
{
    clear();
    if (nx < 2 || x == nullptr)
    {
        if (nx < 2)
        {
            throw std::invalid_argument("nx = " + std::to_string(nx)
                                     + " must be at least 2");
        }
        throw std::invalid_argument("x is NULL");
    }
    double pMean, pStdDev;
    ippsMeanStdDev_64f(x, nx, &pMean, &pStdDev);
    if (pStdDev == 0)
    {
        throw std::invalid_argument("x cannot be all the same values");
    }
    initialize(pMean, pStdDev); 
}

/// Initialized?
bool ZScore::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Standardizes the data
void ZScore::apply(const int npts, const double x[], double *yIn[])
{
    if (npts < 1){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    double *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    ippsNormalize_64f(x, y, npts, pImpl->mMean64f, pImpl->mStd64f);
}

void ZScore::apply(const int npts, const float x[], float *yIn[])
{
    if (npts < 1){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    float *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    ippsNormalize_32f(x, y, npts, pImpl->mMean32f, pImpl->mStd32f);
}
