#include <string>
#include <vector>
#include "private/throw.hpp"
#include "rtseis/trigger/waterLevel.hpp"

using namespace RTSeis::Trigger;

template<RTSeis::ProcessingMode E, class T>
class WaterLevel<E, T>::WaterLevelImpl
{
public:
    std::vector<std::pair<int, int>> mWindows;
    //std::array<int8_t, 2049> mBitMask;
    double mOnTolerance = 0;
    double mOffTolerance = 0;
    bool mInitialized = false;
};

/// C'tor
template<RTSeis::ProcessingMode E, class T>
WaterLevel<E, T>::WaterLevel() :
    pImpl(std::make_unique<WaterLevelImpl> ())
{
}

/// Copy c'tor
template<RTSeis::ProcessingMode E, class T>
WaterLevel<E, T>::WaterLevel(const WaterLevel &trigger)
{
    *this = trigger;
}

/// Move c'tor
template<RTSeis::ProcessingMode E, class T>
WaterLevel<E, T>::WaterLevel(WaterLevel &&trigger) noexcept
{
    *this = std::move(trigger);
}

/// Copy assignment operator
template<RTSeis::ProcessingMode E, class T>
WaterLevel<E, T>& WaterLevel<E, T>::operator=(const WaterLevel &trigger)
{
    if (&trigger == this){return *this;}
    pImpl = std::make_unique<WaterLevelImpl> (*trigger.pImpl);
    return *this;
}

/// Move assignment operator
template<RTSeis::ProcessingMode E, class T>
WaterLevel<E, T>& WaterLevel<E, T>::operator=(WaterLevel &&trigger) noexcept
{
    if (&trigger == this){return *this;}
    pImpl = std::move(trigger.pImpl);
    return *this;
}

/// Destructor
template<RTSeis::ProcessingMode E, class T>
WaterLevel<E, T>::~WaterLevel() = default;

/// Clears the class
template<RTSeis::ProcessingMode E, class T>
void WaterLevel<E, T>::clear() noexcept
{
    pImpl->mWindows.clear();
    pImpl->mOnTolerance = 0;
    pImpl->mOffTolerance = 0;
    pImpl->mInitialized = false;
}

/// Gets the number of triggers
template<RTSeis::ProcessingMode E, class T>
int WaterLevel<E, T>::getNumberOfWindows() const noexcept
{
    return static_cast<int> (pImpl->mWindows.size());
}

/// Initialize the class
template<RTSeis::ProcessingMode E, class T>
void WaterLevel<E, T>::initialize(const double onTolerance,
                               const double offTolerance)
{
    clear();
    pImpl->mOnTolerance = onTolerance;
    pImpl->mOffTolerance = offTolerance;
    pImpl->mInitialized = true;
}

/// Is the class initialized?
template<RTSeis::ProcessingMode E, class T>
bool WaterLevel<E, T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Gets the triggers
template<RTSeis::ProcessingMode E, class T>
std::vector<std::pair<int, int>> WaterLevel<E, T>::getWindows() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mWindows;
}

/// Gets the triggers
template<RTSeis::ProcessingMode E, class T>
void WaterLevel<E, T>::getWindows(
    const int nWindows,
    std::pair<int, int> *windowsIn[]) const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    auto nwinRef = getNumberOfWindows();
    if (nWindows != nwinRef)
    {
        throw std::invalid_argument("nWindows = " + std::to_string(nWindows)
                                 + " must equal " + std::to_string(nwinRef));
    }
    auto windows = *windowsIn;
    std::copy(pImpl->mWindows.begin(), pImpl->mWindows.end(), windows);
}

/// Applies
template<RTSeis::ProcessingMode E, class T>
void WaterLevel<E, T>::apply(const int nSamples, const T x[])
{
    pImpl->mWindows.clear();
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (nSamples < 1){return;} // Nothing to do
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    // Look for all points where the value of x is above the tolerances
    pImpl->mWindows.reserve(2048);
    T on  = static_cast<T> (pImpl->mOnTolerance);
    T off = static_cast<T> (pImpl->mOffTolerance);
    int isOn =-1;
    if (x[0] > on){isOn = 0;}
    for (int i = 1; i < nSamples; ++i)
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

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class RTSeis::Trigger::WaterLevel<RTSeis::ProcessingMode::POST_PROCESSING, double>;
template class RTSeis::Trigger::WaterLevel<RTSeis::ProcessingMode::POST_PROCESSING, float>;
