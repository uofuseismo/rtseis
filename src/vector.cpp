#include <iostream>
#include <complex>
#include <vector>
#include <boost/align.hpp>
#include "rtseis/vector.hpp"

using namespace RTSeis;

template<class T>
class Vector<T>::VectorImpl
{
public:
    std::vector<T, boost::alignment::aligned_allocator<T, 64>> mX;
    //std::vector<double, boost::alignment::aligned_allocator<double, 64> > vector;
};

/// Constructor
template<class T>
Vector<T>::Vector() :
    pImpl(std::make_unique<VectorImpl> ())
{
}

/// Copy constructor
template<class T>
Vector<T>::Vector(const Vector<T> &v)
{
    *this = v;
}

/// Move constructor
template<class T>
Vector<T>::Vector(Vector<T> &&v) noexcept
{
    *this = std::move(v);
}

/// Construct from vector
template<class T>
Vector<T>::Vector(const std::vector<T> &v) :
    pImpl(std::make_unique<VectorImpl> ())
{
    pImpl->mX.resize(v.size());
    std::copy(v.begin(), v.end(), pImpl->mX.begin()); 
}

/// Copy assignment
template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T> &v)
{
    if (&v == this){return *this;}
    pImpl = std::make_unique<VectorImpl> (*v.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Vector<T>& Vector<T>::operator=(Vector<T> &&v) noexcept
{
    if (&v == this){return *this;}
    pImpl = std::move(v.pImpl);
    return *this;
}

/// Resize
template<class T>
void Vector<T>::resize(const size_t n)
{
    pImpl->mX.resize(n);
}

template<class T>
void Vector<T>::resize(const size_t n, const T value)
{
    pImpl->mX.resize(n, value);
}

/// Reserve
template<class T>
void Vector<T>::reserve(const size_t n)
{
    pImpl->mX.reserve(n);
}

/// Data
template<class T>
T* Vector<T>::data() noexcept
{
    return pImpl->mX.data();
}

template<class T>
const T* Vector<T>::data() const noexcept
{
    return pImpl->mX.data();
}

/// Size
template<class T>
size_t Vector<T>::size() const noexcept
{
    return pImpl->mX.size();
}

/// Release memory
template<class T>
void Vector<T>::clear() noexcept
{
    pImpl->mX.clear();
}

/// Empty?
template<class T>
bool Vector<T>::empty() const noexcept
{
    return pImpl->mX.empty();
}

/// Destructor
template<class T>
Vector<T>::~Vector() = default;

/// Iterators
template<class T>
Vector<T>::iterator Vector<T>::begin()
{
    return pImpl->mX.begin();
}

template<class T>
Vector<T>::iterator Vector<T>::end()
{
    return pImpl->mX.end();
}

template<class T>
Vector<T>::const_iterator Vector<T>::cbegin() const
{
    return pImpl->mX.cbegin();
}

template<class T>
Vector<T>::const_iterator Vector<T>::cend() const
{
    return pImpl->mX.cend();
}

template<class T>
T& Vector<T>::operator[](const size_t index)
{
    return pImpl->mX[index];
}

template<class T>
T& Vector<T>::operator[](const size_t index) const
{
    return pImpl->mX[index];
}

template<class T>
T& Vector<T>::at(size_t index)
{
    return pImpl->mX.at(index);
}

template<class T>
T& Vector<T>::at(const size_t index) const
{
    return pImpl->mX.at(index);
}

template class RTSeis::Vector<double>;
template class RTSeis::Vector<float>;
template class RTSeis::Vector<int>;
template class RTSeis::Vector<std::complex<double>>;
template class RTSeis::Vector<std::complex<float>>;
