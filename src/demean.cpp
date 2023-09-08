#include <numeric>
#include "rtseis/demean.hpp"
#include "rtseis/vector.hpp"

using namespace RTSeis;

template<class T>
class Demean<T>::DemeanImpl
{
public:
    double mMean{0};
};

/// Constructor
template<class T>
Demean<T>::Demean() :
    pImpl(std::make_unique<DemeanImpl> ())
{
}

/// Copy constructor
template<class T>
Demean<T>::Demean(const Demean &d) :
    System::ISystem<T, T> (d)
{
    *this = d;
}

/// Move constructor
template<class T>
Demean<T>::Demean(Demean &&d) noexcept :
    System::ISystem<T, T> (d)
{
    *this = std::move(d);
}

/// Copy assignment 
template<class T>
Demean<T>& Demean<T>::operator=(const Demean<T> &demean)
{
    if (&demean == this){return *this;}
    pImpl = std::make_unique<DemeanImpl> (*demean.pImpl);
    return *this;
}

/// Move assignment 
template<class T>
Demean<T>& Demean<T>::operator=(Demean<T> &&demean) noexcept
{
    if (&demean == this){return *this;}
    pImpl = std::move(demean.pImpl);
    return *this;
}

/// Reset class
template<class T>
void Demean<T>::clear() noexcept
{
    System::ISystem<T, T>::clear();
    pImpl = std::make_unique<DemeanImpl> ();
}

/// Destructor
template<class T>
Demean<T>::~Demean() = default;

/// Initialized?
template<class T>
bool Demean<T>::isInitialized() const noexcept
{
    return true;
}

/// Apply
template<class T>
void Demean<T>::apply()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    // Get handle on input
    pImpl->mMean = 0;
    auto x = this->getInputReference();
    if (x.empty()){throw std::runtime_error("No samples in input signal");}
    auto n = x.size();
    // Tabulate the mean 
    const T *__attribute__((aligned(64))) __restrict__ xPtr = x.data();
    double mean = std::accumulate(xPtr, xPtr + n, 0.0);
    mean = mean/n;
    // Create result space
    RTSeis::Vector<T> y;
    y.resize(n);
    // Remove mean: y = x - \bar{x}
    T *__attribute__((aligned(64))) __restrict__ yPtr = y.data();
    std::transform(xPtr, xPtr + n, yPtr,
                   std::bind2nd(std::minus<T> (), mean));
    this->setOutput(std::move(y));
    pImpl->mMean = mean; 
}

/// Get the mean
template<class T>
double Demean<T>::getMean() const noexcept
{
    return pImpl->mMean;
}

/// Template instantiation
template class RTSeis::Demean<double>;
template class RTSeis::Demean<float>;
