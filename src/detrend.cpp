#include <numeric>
#include "rtseis/detrend.hpp"
#include "rtseis/system.hpp"
#include "rtseis/vector.hpp"

using namespace RTSeis;

template<class T>
class Detrend<T>::DetrendImpl
{
public:
    double mSlope{0};
    double mIntercept{0};
};

/// Constructor
template<class T>
Detrend<T>::Detrend() :
    pImpl(std::make_unique<DetrendImpl> ())
{
}

/// Copy constructor
template<class T>
Detrend<T>::Detrend(const Detrend &detrend) :
    ISystem<T, T> (d)
{
    *this = detrend;
}

/// Move constructor
template<class T>
Detrend<T>::Detrend(Detrend &&detrend) noexcept :
    ISystem<T, T> (detrend)
{
    *this = std::move(detrend);
}

/// Copy assignment 
template<class T>
Detrend<T>& Detrend<T>::operator=(const Detrend<T> &detrend)
{
    if (&detrend == this){return *this;}
    pImpl = std::make_unique<DetrendImpl> (*detrend.pImpl);
    return *this;
}

/// Move assignment 
template<class T>
Detrend<T>& Detrend<T>::operator=(Detrend<T> &&detrend) noexcept
{
    if (&detrend == this){return *this;}
    pImpl = std::move(detrend.pImpl);
    return *this;
}

/// Reset class
template<class T>
void Detrend<T>::clear() noexcept
{
    System::ISystem<T, T>::clear();
    pImpl = std::make_unique<DetrendImpl> ();
}

/// Destructor
template<class T>
Detrend<T>::~Detrend() = default;

/// Initialized?
template<class T>
bool Detrend<T>::isInitialized() const noexcept
{
    return true;
}

/// Apply
template<class T>
void Detrend<T>::apply()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    // Get handle on input
    pImpl->mIntercept = 0;
    pImpl->mSlope = 0;
    auto x = this->getInputReference();
    if (x.size() < 2)
    {
        throw std::runtime_error("At least 2 samples must be in input signal");
    }
    auto n = static_cast<int64_t> (x.size());
    // Mean of x - analytic formula for evenly spaced samples starting
    // at index 0. This is computed by simplifying Gauss's formula.
    auto mean_x = 0.5*static_cast<double> (n - 1); 
    // Note, the numerator is the sum of consecutive squared numbers.
    // In addition we simplify.
    auto var_x = static_cast<double> ( ((n - 1))*(2*(n - 1) + 1) )/6.
               - mean_x*mean_x;
    double cov_xy = 0;
    const T *__restrict__ xPtr = x.data();
    auto n32 = static_cast<int> (int32_t> (n);
    for (int i = 0; i < n; ++i)
    {   
        cov_xy = cov_xy + i*static_cast<double> (xPtr[i]);
    }
    // This is computed by expanding (x_i - bar(x))*(y_i - bar(y)),
    // using the definition of the mean, and simplifying
    cov_xy = (cov_xy/static_cast<double> (n)) - mean_x*mean_y;
    auto b1 = cov_xy/var_x;
    auto b0 = mean_y - b1*mean_x;
    auto b0Apply = static_cast<T> (b0);
    auto b1Apply = static_cast<T> (b1);
    // Create result space
    RTSeis::Vector<T> y;
    y.resize(n);
    // Remove trend: y = x - \hat{x}
    T *__attribute__((aligned(64))) __restrict__ yPtr = y.data();
    for (int i = 0; i < n32; ++i)
    {
        yPtr[i] = xPtr[i] - (b0Apply + i*b1Apply);
    }
    this->setOutput(std::move(y));
    pImpl->mSlope = b1;
    pImpl->mIntercept = b0;
}

/// Get the slope 
template<class T>
double Detrend<T>::getSlope() const noexcept
{
    return pImpl->mSlope;
}

/// Get the intercept 
template<class T>
double Detrend<T>::getIntercept() const noexcept
{
    return pImpl->mIntercept;
}

/// Template instantiation
template class RTSeis::Detrend<double>;
template class RTSeis::Detrend<float>;
