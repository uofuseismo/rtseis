#include <cstdio>
#include <cstdlib>
#include <string>
#include <ipps.h>
#include "rtseis/utilities/normalization/zscore.hpp"

using namespace RTSeis::Utilities::Normalization;

namespace
{

void ippsMeanStdDev(const double *x, const int nx,
                    double *pMean, double *pStdDev)
{
    ippsMeanStdDev_64f(x, nx, pMean, pStdDev);
}

void ippsMeanStdDev(const float *x, const int nx,
                    float *pMean, float *pStdDev)
{
    ippsMeanStdDev_32f(x, nx, pMean, pStdDev, ippAlgHintAccurate);
}

void ippsNormalize(const double *x, double *y, const int nx,
                   const double pMean, const double pStd)
{
    ippsNormalize_64f(x, y, nx, pMean, pStd);
}

void ippsNormalize(const float *x, float *y, const int nx,
                   const float pMean, const float pStd)
{
    ippsNormalize_32f(x, y, nx, pMean, pStd);
}

}

template<class T>
class ZScore<T>::ZScoreImpl
{
public:
    T mMean = 0;
    T mStd = 0;
    bool mInitialized = false;
};

/// Constructors
template<class T>
ZScore<T>::ZScore() :
    pImpl(std::make_unique<ZScoreImpl> ())
{
}

/// Copy c'tor
template<class T>
ZScore<T>::ZScore(const ZScore &zscore)
{
    *this = zscore;
}

/// Move c'tor
template<class T>
ZScore<T>::ZScore(ZScore &&zscore) noexcept
{
    *this = std::move(zscore);
}

/// Operators
template<class T>
ZScore<T>& ZScore<T>::operator=(const ZScore<T> &zscore)
{
    if (&zscore == this){return *this;}
    pImpl = std::make_unique<ZScoreImpl> (*zscore.pImpl);
    return *this;
}

template<class T>
ZScore<T>& ZScore<T>::operator=(ZScore<T> &&zscore) noexcept
{
    if (&zscore == this){return *this;}
    pImpl = std::move(zscore.pImpl);
    return *this;
}

/// Destructors
template<class T>
ZScore<T>::~ZScore() = default;

/// Clears the class
template<class T>
void ZScore<T>::clear() noexcept
{
    pImpl->mMean = 0;
    pImpl->mStd = 0;
    pImpl->mInitialized = false;
}

/// Initializes the class
template<class T>
void ZScore<T>::initialize(const double mean,
                           const double std)
{
    clear();
    if (std <= 0)
    {
        throw std::invalid_argument("standard deviation = "
                                  + std::to_string(std) + " must be positive");
    }
    pImpl->mMean = static_cast<T> (mean);
    pImpl->mStd = static_cast<T> (std);
    pImpl->mInitialized = true;
}

/// Initialize the class
template<class T>
void ZScore<T>::initialize(const int nx, const T x[])
{
    clear();
    auto [pMean, pStdDev] = computeMeanAndStandardDeviation(nx, x);
    if (pStdDev == 0)
    {
        throw std::invalid_argument("x cannot be all the same values");
    }
    initialize(static_cast<double> (pMean), static_cast<double> (pStdDev));
}

/// Initialized?
template<class T>
bool ZScore<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Standardizes the data
template<class T>
void ZScore<T>::apply(const int npts, const T x[], T *yIn[])
{
    if (npts < 1){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    T *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    ippsNormalize(x, y, npts, pImpl->mMean, pImpl->mStd);
}

/// Computes the mean and standard deviation
std::pair<double, double>
RTSeis::Utilities::Normalization::computeMeanAndStandardDeviation(
    const int nx, const double x[])
{
    if (nx < 2)
    {
        throw std::invalid_argument("x must contain at least 2 samples");
    }
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    double pMean, pStdDev;
    ippsMeanStdDev(x, nx, &pMean, &pStdDev);
    return std::pair(pMean, pStdDev);
}

std::pair<float, float>
RTSeis::Utilities::Normalization::computeMeanAndStandardDeviation(
    const int nx, const float x[])
{   
    if (nx < 2)
    {
        throw std::invalid_argument("x must contain at least 2 samples");
    }
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    float pMean, pStdDev; 
    ippsMeanStdDev(x, nx, &pMean, &pStdDev);
    return std::pair(pMean, pStdDev);
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///

template class RTSeis::Utilities::Normalization::ZScore<double>;
template class RTSeis::Utilities::Normalization::ZScore<float>;
