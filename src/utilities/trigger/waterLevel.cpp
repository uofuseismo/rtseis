#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <ipps.h>
#include "private/throw.hpp"
#include "rtseis/utilities/trigger/waterLevel.hpp"

using namespace RTSeis::Utilities::Trigger::PostProcessing;

template<class T>
class WaterLevel<T>::WaterLevelImpl
{
public:
    std::vector<std::pair<int, int>> mWindows;
    //std::array<int8_t, 2049> mBitMask;
    double mOnTolerance = 0;
    double mOffTolerance = 0;
    bool mInitialized = false;
};

/// C'tor
template<class T>
WaterLevel<T>::WaterLevel() :
    pImpl(std::make_unique<WaterLevelImpl> ())
{
}

/// Copy c'tor
template<class T>
WaterLevel<T>::WaterLevel(const WaterLevel &trigger)
{
    *this = trigger;
}

/// Move c'tor
template<class T>
WaterLevel<T>::WaterLevel(WaterLevel &&trigger) noexcept
{
    *this = std::move(trigger);
}

/// Copy assignment operator
template<class T>
WaterLevel<T>& WaterLevel<T>::operator=(const WaterLevel &trigger)
{
    if (&trigger == this){return *this;}
    pImpl = std::make_unique<WaterLevelImpl> (*trigger.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
WaterLevel<T>& WaterLevel<T>::operator=(WaterLevel &&trigger) noexcept
{
    if (&trigger == this){return *this;}
    pImpl = std::move(trigger.pImpl);
    return *this;
}

/// Destructor
template<class T>
WaterLevel<T>::~WaterLevel() = default;

/// Clears the class
template<class T>
void WaterLevel<T>::clear() noexcept
{
    pImpl->mWindows.clear();
    pImpl->mOnTolerance = 0;
    pImpl->mOffTolerance = 0;
    pImpl->mInitialized = false;
}

/// Gets the number of triggers
template<class T>
int WaterLevel<T>::getNumberOfWindows() const noexcept
{
    return static_cast<int> (pImpl->mWindows.size());
}

/// Initialize the class
template<class T>
void WaterLevel<T>::initialize(const double onTolerance,
                               const double offTolerance)
{
    clear();
    pImpl->mOnTolerance = onTolerance;
    pImpl->mOffTolerance = offTolerance;
    pImpl->mInitialized = true;
}

/// Is the class initialized?
template<class T>
bool WaterLevel<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Gets the triggers
template<class T>
std::vector<std::pair<int, int>> WaterLevel<T>::getWindows() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mWindows;
}

/// Gets the triggers
template<class T>
void WaterLevel<T>::getWindows(const int nWindows,
                               std::pair<int, int> *windowsIn[]) const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    auto nwinRef = getNumberOfWindows();
    if (nWindows != nwinRef)
    {
        RTSEIS_THROW_IA("nWindows = %d must equal %d", nWindows, nwinRef);
    }
    auto windows = *windowsIn;
    std::copy(pImpl->mWindows.begin(), pImpl->mWindows.end(), windows);
}

/// Applies
template<class T>
void WaterLevel<T>::apply(const int nSamples, const T x[])
{
    pImpl->mWindows.clear();
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (nSamples < 1){return;} // Nothing to do
    if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL\n");}
    // Look for all points where the value of x is above the tolerances
    pImpl->mWindows.reserve(2048);
    T on  = static_cast<T> (pImpl->mOnTolerance);
    T off = static_cast<T> (pImpl->mOffTolerance);
    int isOn =-1;
    if (x[0] > on){isOn = 0;}
    for (int i=1; i<nSamples; ++i)
    {
        // Searching for end of window
        if (isOn > -1)
        {
            if (x[i-1] >= off && x[i] < off)
            {
                pImpl->mWindows.push_back(std::pair(isOn, i));
                isOn = -1;
            }
        }
        // Searching for start of window
        else
        {
            if (x[i-1] < on && x[i] >= on){isOn = i;}
        }
    }
/*
    bool isOn = false;
    bool isOff = false;
    bool mStart = true;
    auto mBitMask = pImpl->mBitMask.data();
    std::pair<int, int> onset(0, 0);
    size_t k = 0;
    int chunkSize = static_cast<int> (pImpl->mBitMask.size()) - 1;
    mBitMask[0] =-1;
    for (int i=0; i<nSamples; i=i+chunkSize)
    {
        int j1 = i;
        int j2 = std::min(nSamples, j1 + chunkSize);
        int nloc = j2 - j1;
        auto first = mBitMask[0];
        int check = 0;
        #pragma omp simd reduction(+:check)
        for (int j=0; j<nloc; ++j)
        {
            isOn = false;
            isOff = false;
            if (x[i+j] > on)  isOn  = true;
            if (x[i+j] > off) isOff = true;
            mBitMask[j+1] = isOn && isOff ? +1 : -1;
            if (mBitMask[j+1] != first){check = 1;}
        }
        if (check == 0){continue;} // Not interested in this window
        // Now we look for all the sign switches in this window
        for (int j=1; j<nloc; ++j)
        {
            if (mBitMask[j-1]*mBitMask[j] ==-1)
            {
                if (mStart)
                {
                    onset.first = i + j - 1;
                }
                else
                {
                    onset.second = i + j - 1;
                    pImpl->mWindows.push_back(onset);
                }
                mStart = !mStart;
            }
        }
        // Update the final conditions
        mBitMask[0] = mBitMask[nloc];
    }
*/
}

/// Template instantiation
template class RTSeis::Utilities::Trigger::PostProcessing::WaterLevel<double>;
template class RTSeis::Utilities::Trigger::PostProcessing::WaterLevel<float>;
