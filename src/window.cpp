#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#include "rtseis/window.hpp"
#include "rtseis/vector.hpp"
#include "src/filterDesign/windowFunctions.hpp"

using namespace RTSeis;

template<class T>
class Window<T>::WindowImpl
{
public:
    void design()
    {
        if (mLength > 0)
        {
            mWindow.resize(mLength, 0);
            T *window = mWindow.data();
            if (mType == Window<T>::Type::Hanning)
            {
                ::hann(mLength, &window);
            }
            else if (mType == Window<T>::Type::Blackman)
            {
                ::blackman(mLength, &window);
            }
            else if (mType == Window<T>::Type::Hamming)
            {
                ::hamming(mLength, &window);
            }
            else if (mType == Window<T>::Type::Kaiser)
            {
                ::kaiser(mLength, &window, static_cast<T> (mBeta));
            }
            else if (mType == Window<T>::Type::Bartlett)
            {
                ::bartlett(mLength, &window);
            }
            else if (mType == Window<T>::Type::Sine)
            {
                ::sine(mLength, &window);
            }
            else
            {   
#ifndef NDEBUG
                assert(false);
#endif
                throw std::runtime_error("Unhandled window type"); 
            }
        }
        else
        {
            mWindow.clear();
        }
    }
    Vector<T> mWindow;
    Window::Type mType{Window::Type::Hanning};
    double mBeta{0.5};
    int mLength{0};
    bool mInitialized{false};
};

/// Constructor
template<class T>
Window<T>::Window() :
    pImpl(std::make_unique<WindowImpl> ()) 
{
}

/// Copy constructor
template<class T>
Window<T>::Window(const Window &window)
{
    *this = window;
}

/// Move constructor
template<class T>
Window<T>::Window(Window &&window) noexcept
{
    *this = std::move(window);
}

/// Copy assignment 
template<class T>
Window<T>& Window<T>::operator=(const Window<T> &window)
{
    if (&window == this){return *this;}
    pImpl = std::make_unique<WindowImpl> (*window.pImpl);
    return *this;
}

/// Move assignment 
template<class T>
Window<T>& Window<T>::operator=(Window<T> &&window) noexcept
{
    if (&window == this){return *this;}
    pImpl = std::move(window.pImpl);
    return *this;
}


/// Initialized
template<class T>
bool Window<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Window type
template<class T>
Window<T>::Type Window<T>::getType() const
{
    if (!isInitialized()){throw std::runtime_error("Window not initialized");}
    return pImpl->mType;
}

/// Initialize class
template<class T>
void Window<T>::initialize(const int length,
                           const Window<T>::Type type,
                           const double beta)
{
    if (length < 1)
    {
        throw std::invalid_argument("Window length = "
                                  + std::to_string(length)
                                  + "  must be positive");
    }
    if (type == Window<T>::Type::Kaiser)
    {
        if (beta < 0)
        {
            throw std::invalid_argument("beta must be positive");
        }
    }
    clear();
    pImpl->mLength = length;
    pImpl->mType = type;
    if (beta >= 0){pImpl->mBeta = beta;}
    pImpl->design();
    pImpl->mInitialized = true;
} 


/// The window
template<class T>
Vector<T> Window<T>::getWindow() const
{
    return Vector<T> {getWindowReference()};
}

template<class T>
const Vector<T> &Window<T>::getWindowReference() const
{
    if (!isInitialized()){throw std::runtime_error("Window not initialized");}
    return *&pImpl->mWindow;
}

/// Reset class
template<class T>
void Window<T>::clear() noexcept
{
    pImpl = std::make_unique<WindowImpl> ();
}

/// Destructor
template<class T>
Window<T>::~Window() = default;

/// Template instantiation
template class RTSeis::Window<double>;
template class RTSeis::Window<float>;
