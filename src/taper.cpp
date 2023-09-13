#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#include "rtseis/taper.hpp"
#include "rtseis/window.hpp"
#include "rtseis/vector.hpp"

using namespace RTSeis;

namespace
{
template<typename T>
void apply(const RTSeis::Vector<T> &window, RTSeis::Vector<T> &y)
{
    auto signalLength = y.size();
    if (signalLength == 0){return;}
    if (signalLength == 1)
    {
        if (!window.empty())
        {
            y[0] = window.at(0)*y.at(0);
        }
    }
    else if (signalLength == 2)
    {
        if (window.size() == 2)
        {
            y[0] = window.at(0)*y.at(0);
            y[1] = window.at(1)*y.at(1);
        }
    }
    else
    {
        // Taper first (m+1)/2 points
        auto halfWindow = static_cast<int> (window.size()/2);
        std::transform(window.cbegin(), window.cbegin() + halfWindow,
                       y.begin(),
                       y.begin(),
                       std::multiplies<T> ());  
        auto i0 = static_cast<int> (y.size() - halfWindow);
        std::transform(window.cend() - halfWindow,
                       window.cend(),
                       y.begin() + i0,
                       y.begin() + i0,
                       std::multiplies<T> ());
    }
}
}

template<class T>
class Taper<T>::TaperImpl
{
public:
    RTSeis::Window<T> mWindow;
    Taper<T>::Window mWindowType;
    double mPercentage{5};
    double mBeta{0.5};
    int mWindowLength{0};
    int mSignalLength{0};
    bool mDoDesign{true};
    bool mInitialized{false};
};

/// Constructor
template<class T>
Taper<T>::Taper() :
    pImpl(std::make_unique<TaperImpl> ()) 
{
}

/// Copy constructor
template<class T>
Taper<T>::Taper(const Taper &taper) :
    System::ISystem<T, T> (taper) 
{
    *this = taper;
}

/// Move constructor
template<class T>
Taper<T>::Taper(Taper &&taper) noexcept :
    System::ISystem<T, T> (taper) 
{
    *this = std::move(taper);
}

/// Copy assignment 
template<class T>
Taper<T>& Taper<T>::operator=(const Taper<T> &taper)
{
    if (&taper == this){return *this;}
    pImpl = std::make_unique<TaperImpl> (*taper.pImpl);
    return *this;
}

/// Move assignment 
template<class T>
Taper<T>& Taper<T>::operator=(Taper<T> &&taper) noexcept
{
    if (&taper == this){return *this;}
    pImpl = std::move(taper.pImpl);
    return *this;
}

/// Reset class
template<class T>
void Taper<T>::clear() noexcept
{
    System::ISystem<T, T>::clear();
    pImpl = std::make_unique<TaperImpl> ();
}

/// Destructor
template<class T>
Taper<T>::~Taper() = default;

/// Initialized
template<class T>
bool Taper<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Taper type
template<class T>
Taper<T>::Window Taper<T>::getWindow() const
{
    if (!isInitialized()){throw std::runtime_error("Taper not initialized");}
    return pImpl->mWindowType;
}

/// Percentage
template<class T>
double Taper<T>::getPercentage() const
{
    if (!isInitialized()){throw std::runtime_error("Taper not initialized");}
    return pImpl->mPercentage;
}

/// Apply
template<class T>
void Taper<T>::apply()
{
    if (!isInitialized()){throw std::runtime_error("Taper not initialized");}
    auto x = System::ISystem<T, T>::getInputReference();
    if (pImpl->mSignalLength == 0 || getPercentage() == 0)
    {
        System::ISystem<T, T>::setOutput(x);
        return;
    }
#ifndef NDEBUG
    assert(!x.empty());
    assert(!pimpl->mDoDesign);
#endif
    // Copy constructor - all that's left to do now is multiply by window
    // the window function
    Vector<T> y(x);
    if (pImpl->mWindowType == Window::Boxcar)
    {
        RTSeis::Vector<T> zeros;
        zeros.resize(pImpl->mWindowLength, 0);
        ::apply<T>(zeros, y);
    }
    else
    {
        const auto &window = pImpl->mWindow.getWindowReference();
        ::apply<T>(window, y);
    }
    System::ISystem<T, T>::setOutput(y);
}

/// Initialize class
template<class T>
void Taper<T>::initialize(const Window windowType,
                          const double percentage,
                          const double beta)
{
    if (percentage < 0 || percentage > 100)
    {
        throw std::invalid_argument("Percentage = "
                                  + std::to_string(percentage)
                                  + "  must be in range [0,100]");
    }
    if (windowType == Window::Kaiser)
    {
        if (beta < 0)
        {
            throw std::invalid_argument("Beta must be postiive");
        }
    } 
    clear();
    pImpl->mWindowType = windowType;
    pImpl->mPercentage = percentage;
    pImpl->mInitialized = true;
} 

/// Set data
template<class T>
void Taper<T>::setInput(const Vector<T> &signal)
{
    if (!isInitialized()){throw std::runtime_error("Taper not initialized");}
    if (signal.empty()){throw std::invalid_argument("Signal is empty");}
    auto signalLength = static_cast<int> (signal.size());
    if (signalLength != pImpl->mSignalLength || pImpl->mDoDesign)
    {
        pImpl->mDoDesign = true;
        double percentage = getPercentage(); 
        if (percentage == 0)
        {
            pImpl->mWindowLength = 0;
        }
        else if (std::abs(percentage - 100) <
                 std::numeric_limits<double>::epsilon())
        {
            pImpl->mWindowLength = signal.size();
        }
        else
        {
            auto nPercent
                = static_cast<int> (std::round((signalLength*percentage)/100.))
                + 1;
            pImpl->mWindowLength
                = std::max(2, std::min(signalLength, nPercent));
        }
        if (pImpl->mWindowType != Window::Boxcar)
        {
            if (pImpl->mWindowType == Window::Hamming)
            {
                pImpl->mWindow.initialize(pImpl->mWindowLength,
                                          WindowType::Hamming,
                                          pImpl->mBeta);
            }
            else if (pImpl->mWindowType == Window::Blackman)
            {
                pImpl->mWindow.initialize(pImpl->mWindowLength,
                                          WindowType::Blackman,
                                          pImpl->mBeta);
            }
            else if (pImpl->mWindowType == Window::Hanning)
            {
                pImpl->mWindow.initialize(pImpl->mWindowLength,
                                          WindowType::Hanning,
                                          pImpl->mBeta);
            }
            else if (pImpl->mWindowType == Window::Sine)
            {
                pImpl->mWindow.initialize(pImpl->mWindowLength,
                                          WindowType::Sine,
                                          pImpl->mBeta);
            }
            else if (pImpl->mWindowType == Window::Bartlett)
            {
                pImpl->mWindow.initialize(pImpl->mWindowLength,
                                          WindowType::Bartlett,
                                          pImpl->mBeta);
            }
            else if (pImpl->mWindowType == Window::Kaiser)
            {
                pImpl->mWindow.initialize(pImpl->mWindowLength,
                                          WindowType::Kaiser,
                                          pImpl->mBeta);
            }
            else
            {
#ifndef NDEBUG
                assert(false);
#endif
                throw std::runtime_error("Unhandled window function in taper");
            }
#ifndef NDEBUG
            assert(pImpl->mWindow.isInitialized());
#endif
        }
        else
        {
            pImpl->mWindow.clear();
        }
        pImpl->mDoDesign = false;
    }
    System::ISystem<T, T>::setInput(signal);
    pImpl->mSignalLength = System::ISystem<T, T>::getInputReference().size(); 
}

/// Template instantiation
template class RTSeis::Taper<double>;
template class RTSeis::Taper<float>;
