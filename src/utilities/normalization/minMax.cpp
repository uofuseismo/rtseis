#include <stdexcept>
#include <string>
#include <algorithm>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/utilities/normalization/minMax.hpp"


using namespace RTSeis::Utilities::Normalization;

template<class T>
class MinMax<T>::MinMaxImpl
{
public:
    T mDataMin = 0;
    T mDataMax = 0;
    T mTargetMin = 0;
    T mTargetMax = 0;
    T mDen = 0;
    bool mRescaleToUnitInterval = true; 
    bool mInitialized = false;
};

/// Constructors
template<class T>
MinMax<T>::MinMax() :
    pImpl(std::make_unique<MinMaxImpl> ())
{
}

template<class T>
[[maybe_unused]] MinMax<T>::MinMax(const MinMax<T> &minMax)
{
    *this = minMax;
}

template<class T>
MinMax<T>::MinMax(MinMax &&minMax) noexcept
{
    *this = std::move(minMax);
}

/// Operators
template<class T>
MinMax<T>& MinMax<T>::operator=(const MinMax &minMax)
{
    if (&minMax == this){return *this;}
    pImpl = std::make_unique<MinMaxImpl> (*minMax.pImpl);
    return *this;
}

template<class T>
MinMax<T>& MinMax<T>::operator=(MinMax &&minMax) noexcept
{
    if (&minMax == this){return *this;}
    pImpl = std::move(minMax.pImpl);
    return *this;
}

/// Destructors
template<class T>
MinMax<T>::~MinMax() = default;

template<class T>
void MinMax<T>::clear() noexcept
{
    pImpl->mDataMin = 0;
    pImpl->mDataMax = 0;
    pImpl->mTargetMin = 0;
    pImpl->mTargetMax = 0;
    pImpl->mDen = 0;
    pImpl->mRescaleToUnitInterval = true;
    pImpl->mInitialized = false;
}

/// Initializers
template<class T>
void MinMax<T>::initialize(const std::pair<T, T> dataMinMax,
                           const std::pair<T, T> targetMinMax)
{
    clear();
    if (dataMinMax.first == dataMinMax.second)
    {
        auto errmsg = "dataMin = " + std::to_string(dataMinMax.first)
                    + " cannot equal dataMax = "
                    + std::to_string(dataMinMax.second);
        throw std::invalid_argument(errmsg);
    }
    pImpl->mDataMin = dataMinMax.first;
    pImpl->mDataMax = dataMinMax.second;
    pImpl->mTargetMin = targetMinMax.first;
    pImpl->mTargetMax = targetMinMax.second;;
    pImpl->mDen = (pImpl->mDataMax - pImpl->mDataMin);
    pImpl->mRescaleToUnitInterval = true;
    if (pImpl->mTargetMin != 0 || pImpl->mTargetMax != 1)
    {
        pImpl->mRescaleToUnitInterval = false;
    }
    pImpl->mInitialized = true;
}

template<class T>
void MinMax<T>::initialize(const int npts, const T x[],
                           const std::pair<T, T> targetMinMax)
{
    clear();
    if (npts < 2 || x == nullptr)
    {
        if (npts < 2)
        {
            throw std::invalid_argument("npts = " + std::to_string(npts)
                                      + " must be at least 2");
        }
        throw std::invalid_argument("x cannot be NULL");
    }
    // Find the min and max
    const auto [dataMin, dataMax] = std::minmax_element(x, x+npts);
    /*
    double dataMin, dataMax;
    ippsMinMax_64f(x, npts, &dataMin, &dataMax);
    */
    if (*dataMin == *dataMax)
    {
        throw std::invalid_argument("All elements of x are identical");
    }
    std::pair<T, T> dataMinMax(*dataMin, *dataMax);
    initialize(dataMinMax, targetMinMax);
}        

/// Apply min max normalization
template<>
void MinMax<double>::apply(const int npts,
                           const double x[], double *yIn[]) const
{
    if (npts < 1){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    double *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Rescale
    if (pImpl->mRescaleToUnitInterval)
    {
        // y = 0 + (x - xmin)/(xmax - xmin)*(1 - 0)
        ippsNormalize_64f(x, y, npts, pImpl->mDataMin, pImpl->mDen);
    }
    else
    {
        // y = a + (x - xmin)/(xmax - xmin)*(b - a)
        //   = a - xmin/(xmax - xmin)*(b - a) + x*(b - a)/(xmax - xmin)
        //   = \alpha + \beta x 
        auto bma = pImpl->mTargetMax - pImpl->mTargetMin;
        auto alpha = pImpl->mTargetMin - (pImpl->mDataMin*bma)/pImpl->mDen;
        auto beta  = bma/pImpl->mDen;
        #pragma omp simd
        for (auto i=0; i<npts; ++i)
        {
            y[i] = alpha + beta*x[i];
        }
    }
}

template<>
void MinMax<float>::apply(const int npts, const float *x, float *yIn[]) const
{
    if (npts < 1){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    float *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Rescale
    if (pImpl->mRescaleToUnitInterval)
    {
        // y = 0 + (x - xmin)/(xmax - xmin)*(1 - 0)
        ippsNormalize_32f(x, y, npts, pImpl->mDataMin, pImpl->mDen);
    }
    else
    {
        auto bma = pImpl->mTargetMax - pImpl->mTargetMin;
        auto alpha = pImpl->mTargetMin - (pImpl->mDataMin*bma)/pImpl->mDen;
        auto beta = bma/pImpl->mDen;
        #pragma omp simd
        for (auto i = 0; i < npts; ++i)
        {
            y[i] = alpha + beta * x[i];
        }
    }
}

template<class T>
bool MinMax<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Template instantiation
template class RTSeis::Utilities::Normalization::MinMax<double>;
template class RTSeis::Utilities::Normalization::MinMax<float>;
