#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <array>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/characteristicFunction/iirKurtosis.hpp"

namespace RealTime = RTSeis::Utilities::CharacteristicFunction::RealTime;
namespace PostProcessing
    = RTSeis::Utilities::CharacteristicFunction::PostProcessing;

template<class T>
class RTSeis::Utilities::CharacteristicFunction::IIRKurtosisImpl
{
public:
    /// Default constructor
    IIRKurtosisImpl(bool lRealTime) :
        mRealTime(lRealTime)
    {
    }
    /// Initializes the class
    void initialize(const T c1)
    {
        mC1 = c1;
        mInitialized = true;
    }
    /// Applies the filter
    void apply(const int n, const T x[], T y[]) noexcept
    {
        // Compute some constants 
        const T two = 2;
        const T one = 1;
        const T half = one/two;
        auto c1 = mC1;
        auto a1 = one - c1; 
        auto c2 = half*(one - a1*a1);
        auto one_p_c1 = one + c1; 
        auto two_c1 = two*c1;
        auto bias =-3.0*c1 - 3.0;
        // Get the delay lines
        auto mu1Delay = mDelay[0];
        auto mu2Delay = mDelay[1];
        auto k4barDelay = mDelay[2];
        // Apply the filter
        for (int i=0; i<n; i++)
        {
             auto mu1 = a1*mu1Delay + c1*x[i]; // Update IIR averaging
             auto dx  = x[i] - mu1Delay;
             auto dx2 = dx*dx;
             auto mu2 = a1*mu2Delay + c2*dx2;  // Update IIR averaging
             dx2 = dx2/mu2Delay;
             // one_p_c1 - two_c1*dx2
             auto xscal = std::fma(-two_c1, dx2, one_p_c1); // one_p_c1 - two_c1*dx2;
             auto dx22 = dx2*dx2;
             auto c1_dx22 = c1*dx22 + bias;
             // xscal*k4barDelay + c1_dx22 + bias
             y[i] = std::fma(xscal, k4barDelay, c1_dx22);
             // Update delay lines 
             mu1Delay = mu1;
             mu2Delay = mu2;
             k4barDelay = y[i]; // This is k4bar
        }
        // Update the final conditions
        if (mRealTime)
        {
            mDelay[0] = mu1Delay;
            mDelay[1] = mu2Delay;
            mDelay[2] = k4barDelay;
        }
        else
        {
            resetInitialConditions();
        }
    }
    /// Reset the initial conditions
    void resetInitialConditions() noexcept
    {
        mDelay = mZi;
    }
    /// Set the initial conditions
    void setInitialConditions(const T mu1, const T mu2, const T k4bar)
    {
        mZi[0] = mu1;
        mZi[1] = mu2;
        mZi[2] = k4bar;
        mDelay = mZi;
    }
    /// Resets the class
    void clear() noexcept
    {
        mDelay = {0, 1, 0};
        mZi = {0, 1, 0};
        mC1 = 0.9;
        mRealTime = false;
        mInitialized = false;
    }
//private:
    /// The delay line.  This is ordered the mean of the first order moment,
    /// the mean of the second order moment, and the initial unbiased 
    /// kurtosis value.
    std::array<T, 3> mDelay = {0, 1, 0};
    /// A saved copy of the initial conditions
    std::array<T, 3> mZi = {0, 1, 0};
    /// The pole in the filter
    T mC1 = 0.9;
    /// Indicates whether this is for real-time or post-processing
    bool mRealTime = false;
    /// Indicates whether this is initialized or not
    bool mInitialized = false;
};

//----------------------------------------------------------------------------//
//                                Post-Processing                             //
//----------------------------------------------------------------------------//
/// Default constructor 
template<class T>
PostProcessing::IIRKurtosis<T>::IIRKurtosis() :
    pImpl(std::make_unique<IIRKurtosisImpl<T>> (false))
{
}

/// Copy c'tor
template<class T>
PostProcessing::IIRKurtosis<T>::IIRKurtosis(const IIRKurtosis &iirKurtosis)
{
    *this = iirKurtosis;
}

/// Move c'tor
template<class T>
PostProcessing::IIRKurtosis<T>::IIRKurtosis(IIRKurtosis &&iirKurtosis) noexcept
{
    *this = std::move(iirKurtosis);
}


/// Copy assignment operator
template<class T>
PostProcessing::IIRKurtosis<T>&
PostProcessing::IIRKurtosis<T>::operator=(const IIRKurtosis &iirKurtosis)
{
    if (&iirKurtosis == this){return *this;}
    pImpl = std::make_unique<IIRKurtosisImpl<T>> (*iirKurtosis.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
PostProcessing::IIRKurtosis<T>&
PostProcessing::IIRKurtosis<T>::operator=(IIRKurtosis &&iirKurtosis) noexcept
{
    if (&iirKurtosis == this){return *this;}
    pImpl= std::move(iirKurtosis.pImpl);
    return *this;
}

/// Destructor
template<class T>
PostProcessing::IIRKurtosis<T>::~IIRKurtosis() = default;

/// Resets/clears the class
template<class T>
void PostProcessing::IIRKurtosis<T>::clear() noexcept
{
    pImpl->clear();
}

/// Initializes the class
template<class T>
void PostProcessing::IIRKurtosis<T>::initialize(const T c1)
{
    if (std::abs(c1) >= 1)
    {
        RTSEIS_THROW_IA("|c1| = %e must be less than 1", c1);
    }
    pImpl->initialize(c1);
}

/// Determines if this is initialized
template<class T>
bool PostProcessing::IIRKurtosis<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Applies the filter
template<class T>
void PostProcessing::IIRKurtosis<T>::apply(
    const int nx, const T x[], T *yIn[])
{
    if (nx < 1){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    auto y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    pImpl->apply(nx, x, y);
}

/// Sets the initial conditions
template<class T>
void PostProcessing::IIRKurtosis<T>::setInitialConditions(
    const T mu1, const T mu2, const T k4bar)
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    pImpl->setInitialConditions(mu1, mu2, k4bar);
}

//-----------------------------------------------------------------------------//
//                                 Real Time                                   //
//-----------------------------------------------------------------------------//
/// Default constructor 
template<class T>
RealTime::IIRKurtosis<T>::IIRKurtosis() :
    pImpl(std::make_unique<IIRKurtosisImpl<T>> (true))
{
}

/// Copy c'tor
template<class T>
RealTime::IIRKurtosis<T>::IIRKurtosis(const IIRKurtosis &iirKurtosis)
{
    *this = iirKurtosis;
}

/// Move c'tor
template<class T>
RealTime::IIRKurtosis<T>::IIRKurtosis(IIRKurtosis &&iirKurtosis) noexcept
{
    *this = std::move(iirKurtosis);
}


/// Copy assignment operator
template<class T>
RealTime::IIRKurtosis<T>&
RealTime::IIRKurtosis<T>::operator=(const IIRKurtosis &iirKurtosis)
{
    if (&iirKurtosis == this){return *this;}
    pImpl = std::make_unique<IIRKurtosisImpl<T>> (*iirKurtosis.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
RealTime::IIRKurtosis<T>&
RealTime::IIRKurtosis<T>::operator=(IIRKurtosis &&iirKurtosis) noexcept
{
    if (&iirKurtosis == this){return *this;}
    pImpl= std::move(iirKurtosis.pImpl);
    return *this;
}

/// Destructor
template<class T>
RealTime::IIRKurtosis<T>::~IIRKurtosis() = default;

/// Resets/clears the class
template<class T>
void RealTime::IIRKurtosis<T>::clear() noexcept
{
    pImpl->clear();
} 

/// Initializes the class
template<class T>
void RealTime::IIRKurtosis<T>::initialize(const T c1)
{
    if (std::abs(c1) >= 1)
    {
        RTSEIS_THROW_IA("|c1| = %e must be less than 1", c1);
    }
    pImpl->initialize(c1);
}

/// Determines if this is initialized
template<class T>
bool RealTime::IIRKurtosis<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Applies the filter
template<class T>
void RealTime::IIRKurtosis<T>::apply(
    const int nx, const T x[], T *yIn[])
{
    if (nx < 1){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    auto y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    pImpl->apply(nx, x, y); 
}

/// Sets the initial conditions
template<class T>
void RealTime::IIRKurtosis<T>::setInitialConditions(
    const T mu1, const T mu2, const T k4bar)
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    pImpl->setInitialConditions(mu1, mu2, k4bar);
}

/// Resets the initial conditions
template<class T>
void RealTime::IIRKurtosis<T>::resetInitialConditions()
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    pImpl->resetInitialConditions();
}

//-----------------------------------------------------------------------------//
//                            Template instantiation                           //
//-----------------------------------------------------------------------------//
template class
RTSeis::Utilities::CharacteristicFunction::PostProcessing::IIRKurtosis<double>;
template class
RTSeis::Utilities::CharacteristicFunction::PostProcessing::IIRKurtosis<float>;
template class
RTSeis::Utilities::CharacteristicFunction::RealTime::IIRKurtosis<double>;
template class
RTSeis::Utilities::CharacteristicFunction::RealTime::IIRKurtosis<float>;
