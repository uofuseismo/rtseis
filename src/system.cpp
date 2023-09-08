#include <complex>
#include "rtseis/system/system.hpp"
#include "rtseis/vector.hpp"

using namespace RTSeis::System;

template<class U, class T>
class ISystem<U, T>::ISystemImpl
{
public:
    RTSeis::Vector<U> mX;
    RTSeis::Vector<T> mY;
};

/// Constructor
template<class U, class T>
ISystem<U, T>::ISystem() :
    pImpl(std::make_unique<ISystemImpl> ())
{
}

/// Copy constructor
template<class U, class T>
ISystem<U, T>::ISystem(const ISystem<U, T> &system)
{
    *this = system;
}

/// Move constructor
template<class U, class T>
ISystem<U, T>::ISystem(ISystem<U, T> &&system) noexcept
{
    *this = std::move(system);
}

/// Copy assignment
template<class U, class T>
ISystem<U, T>& ISystem<U, T>::operator=(const ISystem<U, T> &system)
{
    if (&system == this){return *this;}
    pImpl = std::make_unique<ISystemImpl> (*system.pImpl);
    return *this;
}

/// Move assignment
template<class U, class T>
ISystem<U, T>& ISystem<U, T>::operator=(ISystem<U, T> &&system) noexcept
{
    if (&system == this){return *this;}
    pImpl = std::move(system.pImpl);
    return *this;
}

/// Set the input signal
template<class U, class T>
void ISystem<U, T>::setInput(const RTSeis::Vector<U> &x)
{
    pImpl->mX = x;
}

template<class U, class T>
void ISystem<U, T>::setInput(RTSeis::Vector<U> &&x) noexcept
{
    pImpl->mX = std::move(x);
}

/// Get the input signal
template<class U, class T>
RTSeis::Vector<U> ISystem<U, T>::getInput() const
{
    return pImpl->mX;
}

template<class U, class T>
const RTSeis::Vector<U>& ISystem<U, T>::getInputReference() const noexcept
{
    return pImpl->mX;
}

/// Reset class
template<class U, class T> void ISystem<U, T>::clear() noexcept
{
    pImpl = std::make_unique<ISystemImpl> ();
}

/// Destructor
template<class U, class T> ISystem<U, T>::~ISystem() = default;

/// Set the output signal
template<class U, class T>
void ISystem<U, T>::setOutput(const RTSeis::Vector<T> &y)
{
    pImpl->mY = y;
}

template<class U, class T>
void ISystem<U, T>::setOutput(RTSeis::Vector<T> &&y) noexcept
{
    pImpl->mY = std::move(y);
}

/// Get the output signal
template<class U, class T>
RTSeis::Vector<T> ISystem<U, T>::getOutput() const
{
    return pImpl->mY;
}

template<class U, class T>
const RTSeis::Vector<T>& ISystem<U, T>::getOutputReference() const noexcept
{
    return pImpl->mY;
} 

///--------------------------------------------------------------------------///
///                             Template Instantiation                       ///
///--------------------------------------------------------------------------///
template class RTSeis::System::ISystem<double, double>;
template class RTSeis::System::ISystem<float,  float>;
template class RTSeis::System::ISystem<float,  double>;
template class RTSeis::System::ISystem<double, float>;
template class RTSeis::System::ISystem<double, std::complex<double>>;
template class RTSeis::System::ISystem<std::complex<double>, double>;
template class RTSeis::System::ISystem<std::complex<float>, float>;
