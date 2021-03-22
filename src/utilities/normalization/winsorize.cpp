#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <ipps.h>
#include "rtseis/utilities/normalization/winsorize.hpp"

using namespace RTSeis::Utilities::Normalization;

namespace
{
void ippsThreshold_LTValGTVal(const double *pSrc,
                              double *pDst,
                              const int len,
                              const double levelLT,
                              const double valueLT,
                              const double levelGT,
                              const double valueGT)
{
#ifndef NDEBUG
    auto status = 
#endif
    ippsThreshold_LTValGTVal_64f(pSrc, pDst, len,
                                 levelLT, valueLT,
                                 levelGT, valueGT);
#ifndef NDEBUG
    assert(status == ippStsNoErr);
#endif
}

void ippsThreshold_LTValGTVal(const float *pSrc,
                              float *pDst,
                              const int len,
                              const float levelLT,
                              const float valueLT,
                              const float levelGT,
                              const float valueGT)
{
#ifndef NDEBUG
    auto status = 
#endif
    ippsThreshold_LTValGTVal_32f(pSrc, pDst, len,
                                 levelLT, valueLT,
                                 levelGT, valueGT);
#ifndef NDEBUG
    assert(status == ippStsNoErr);
#endif
}

}

template<class T>
class Winsorize<T>::WinsorizeImpl
{
public:
    void createMemory(const int n)
    {
        auto ns = static_cast<size_t> (n);
        if (mLowMemory)
        {
            mWork.resize(ns, 0);
        }
        else
        {
           if (ns > mWork.size())
           {
               mWork.resize(ns, 0);
           }
        }
    }
    std::vector<T> mWork;
    std::pair<double, double> mLimits{0.05, 0.95};
    bool mInclusive = true;
    bool mLowMemory = false;
    bool mInitialized = false;
};

/// C'tor
template<class T>
Winsorize<T>::Winsorize() :
    pImpl(std::make_unique<WinsorizeImpl> ())
{
}

/// Copy c'tor
template<class T>
Winsorize<T>::Winsorize(const Winsorize &winsorize)
{
    *this = winsorize;
}

/// Move c'tor
template<class T>
Winsorize<T>::Winsorize(Winsorize &&winsorize) noexcept
{
    *this = std::move(winsorize);
}

/// Copy assignment
template<class T>
Winsorize<T>& Winsorize<T>::operator=(const Winsorize<T> &winsorize)
{
    if (&winsorize == this){return *this;}
    pImpl = std::make_unique<WinsorizeImpl> (*winsorize.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Winsorize<T>& Winsorize<T>::operator=(Winsorize<T> &&winsorize) noexcept
{
    if (&winsorize == this){return *this;}
    pImpl = std::move(winsorize.pImpl);
    return *this;
}

/// Destructor
template<class T>
Winsorize<T>::~Winsorize() = default;

/// Clear class
template<class T>
void Winsorize<T>::clear() noexcept
{
    pImpl->mWork.clear();
    pImpl->mLimits = std::make_pair<double, double> (0.05, 0.95);
    pImpl->mInclusive = true;
    pImpl->mLowMemory = false;
    pImpl->mInitialized = false;
}

/// Initialize
template<class T>
void Winsorize<T>::initialize(const std::pair<double, double> &limits,
                              const bool inclusive,
                              const bool lowMemory)
{
    clear();
    if (limits.first < 0)
    {
        throw std::invalid_argument("Lower limit = "
                                  + std::to_string(limits.first)
                                  + " must be at least 0");
    }
    if (limits.second > 100)
    {
        throw std::invalid_argument("Upper limit = "
                                  + std::to_string(limits.second)
                                  + " cannot exceed 100");
    }
    if (limits.first > limits.second)
    {
        throw std::invalid_argument("Lower limit = "
                                  + std::to_string(limits.first)
                                  + " must be less than upper limit = "
                                  + std::to_string(limits.second));
    } 
    pImpl->mLimits = std::make_pair<double, double> (limits.first/100,
                                                     limits.second/100);
    pImpl->mInclusive = inclusive;
    pImpl->mLowMemory = lowMemory;
    pImpl->mInitialized = true;
}

/// Initialized?
template<class T>
bool Winsorize<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Winsorize the array 
template<class T>
void Winsorize<T>::apply(const int n, const T *__restrict__ x, T *yIn[])
{
    if (n < 1){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    T *__restrict__ y =  *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Handle an edge case
    if (n < 3)
    {
        std::copy(x, x + n, y);
        return;
    }
    // Figure out the lower/upper indices
    int lowerIndex = 0;
    int upperIndex = n - 1;
    if (pImpl->mInclusive)
    {
        lowerIndex = static_cast<int> (pImpl->mLimits.first*n);
        upperIndex = static_cast<int> (pImpl->mLimits.second*n);
    }
    else
    {
        lowerIndex = static_cast<int> (std::round(pImpl->mLimits.first*n));
        upperIndex = static_cast<int> (std::round(pImpl->mLimits.second*n)) - 1;
    }
    // Some sanity checks
    lowerIndex = std::max(0, lowerIndex);
    upperIndex = std::max(lowerIndex, std::min(n - 1, upperIndex)); 
    // Make sure memory is there to handle calculation
    pImpl->createMemory(n);
    // Copy the input signal and sort
    T *__restrict__ workPtr = pImpl->mWork.data();
    std::copy(x, x + n, workPtr);
    std::sort(workPtr, workPtr + n);
    // Threshold
    T lowerLimit = pImpl->mWork.at(lowerIndex);
    T upperLimit = pImpl->mWork.at(upperIndex);
    ippsThreshold_LTValGTVal(x, y, n,
                             lowerLimit, lowerLimit,
                             upperLimit, upperLimit);
}

/// Template instantiation
template class RTSeis::Utilities::Normalization::Winsorize<double>;
template class RTSeis::Utilities::Normalization::Winsorize<float>;
