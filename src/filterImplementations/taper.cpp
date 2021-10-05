#include <iostream>
#include <string>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <memory>
#include <ipps.h>
#include "rtseis/postProcessing/singleChannel/taper.hpp"
#include "rtseis/filterImplementations/taper.hpp"
#include "rtseis/utilities/windowFunctions.hpp"

/*
using namespace RTSeis::PostProcessing::SingleChannel;

class TaperParameters::TaperParametersImpl
{
public:
    /// Percentage of signal to apply taper to
    double pct = 5;
    /// Taper type.
    Type type = Type::HAMMING;
};

/// Default c'tor
TaperParameters::TaperParameters(const double pct,
                                 const Type type) :
    pImpl(std::make_unique<TaperParametersImpl> ())
{
    setPercentage(pct);
    setTaperType(type);
}

/// Copy c'tor
TaperParameters::TaperParameters(const TaperParameters &parms)
{
    *this = parms;
}

/// Move c'tor
TaperParameters::TaperParameters(TaperParameters &&parms) noexcept
{
    *this = std::move(parms);
}

/// Copy assignment
TaperParameters& TaperParameters::operator=(const TaperParameters &parms)
{
    if (&parms == this){return *this;}
    pImpl = std::make_unique<TaperParametersImpl> (*parms.pImpl);
    return *this;
}

/// Move assignment
TaperParameters& TaperParameters::operator=(TaperParameters &&parms) noexcept
{
    if (&parms == this){return *this;}
    pImpl = std::move(parms.pImpl);
    return *this;
}

/// Destructor
TaperParameters::~TaperParameters() = default;

/// Taper type
void TaperParameters::setTaperType(const Type type) noexcept
{
    pImpl->type = type;
}

TaperParameters::Type TaperParameters::getTaperType() const
{
    return pImpl->type;
}

/// Percentage
void TaperParameters::setPercentage(const double pct)
{
    if (pct < 0 || pct > 100)
    {
        throw std::invalid_argument("Percentage = " + std::to_string(pct)
                                  + " must be in range [0,100]");
    }
    pImpl->pct = pct;
}

double TaperParameters::getPercentage() const
{
    return pImpl->pct;
}

/// Reset class
void TaperParameters::clear() noexcept
{
    if (pImpl)
    {
        pImpl->pct = 5;
        pImpl->type = Type::HAMMING;
    }
}

/// Is valid?
bool TaperParameters::isValid() const noexcept
{
    if (pImpl->pct < 0 || pImpl->pct > 100){return false;}
    return true;
}
*/

//============================================================================//
//                                    Tapering                                //
//============================================================================//

using namespace RTSeis::FilterImplementations;

template<class T>
class Taper<T>::TaperImpl
{
public:
    //TaperParameters parms; 
    std::vector<double> mW8; 
    std::vector<float>  mW4; 
    double mBeta = 8;
    double mPercentage = 5;
    int mWindowLength0 =-1;
    TaperWindowType mWindowType = TaperWindowType::Hamming;
    bool mInitialized = false;
};

/// C'tor
template<class T>
Taper<T>::Taper() :
    pImpl(std::make_unique<TaperImpl>())
{
}

/// Copy c'tor
template<class T>
Taper<T>::Taper(const Taper<T> &taper)
{
    *this = taper;
}

/// Move c'tor
template<class T>
Taper<T>::Taper(Taper<T> &&taper) noexcept
{
    *this = std::move(taper);
}

/*
/// C'tor from parameters
template<class T>
Taper<T>::Taper(const TaperParameters &parameters) :
    pImpl(std::make_unique<TaperImpl>())
{
    initialize(parameters);
}
*/

/// Destructor
template<class T>
Taper<T>::~Taper() = default;

/// Reset class
template<class T>
void Taper<T>::clear() noexcept
{
    pImpl->mW8.clear();
    pImpl->mW4.clear();
    pImpl->mBeta = 8;
    pImpl->mWindowLength0 =-1;
    pImpl->mPercentage = 5;
    pImpl->mWindowType = TaperWindowType::Hamming;
    pImpl->mInitialized = false;
}

/// Copy assignment
template<class T>
Taper<T>& Taper<T>::operator=(const Taper &taper)
{
    if (&taper == this){return *this;}
    pImpl = std::make_unique<TaperImpl> (*taper.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Taper<T>& Taper<T>::operator=(Taper &&taper) noexcept
{
    if (&taper == this){return *this;}
    pImpl = std::move(taper.pImpl);
    return *this;
}

/// Initialize
template<class T>
void Taper<T>::initialize(const double percentage,
                          const TaperWindowType type,
                          const double beta)
{
    clear(); // Sets mWindowLength0 to -1
    if (percentage < 0 || percentage > 100)
    {
        throw std::invalid_argument("Percentage = " + std::to_string(percentage)
                                  + " must be in range [0,100]");
    }
    if (type == TaperWindowType::Kaiser && beta < 0)
    {
        throw std::invalid_argument("Beta must be positive for Kaiser windows");
    }
    pImpl->mPercentage = percentage;
    pImpl->mBeta = beta;
    pImpl->mWindowType = type;
    pImpl->mInitialized = true;
}

/// Apply double
template<>
void Taper<double>::apply(const int nx, const double x[], double *yIn[])
{
    if (nx <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Taper not initialized");}
    double *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Deal with an edge case
    if (nx < 3)
    {
        y[0] = 0;
        if (nx == 2){y[1] = 0;}
        return;
    }
    // Compute taper length
    double pct = pImpl->mPercentage;
    int npct = static_cast<int> (std::round((nx*pct)/100.)) + 1;
    int m = std::max(2, std::min(nx, npct));
    // Redesign the window?  If the parameters were (re)set then mWindowLength0
    // is -1.
    // Otherwise, if the same length signal is coming at us then the precision
    // of the module can't change so we can just use the old window.
    if (pImpl->mWindowLength0 != m)
    {
        pImpl->mW8.resize(m); // Prevent reallocations in a bit
        auto type = pImpl->mWindowType;
        double beta = pImpl->mBeta;
        double *w8data = pImpl->mW8.data();
        if (type == TaperWindowType::Hamming)
        {
            RTSeis::Utilities::WindowFunctions::hamming(m, &w8data);
        }
        else if (type == TaperWindowType::Blackman)
        {
            RTSeis::Utilities::WindowFunctions::blackman(m, &w8data);
        }
        else if (type == TaperWindowType::Hann)
        {
            RTSeis::Utilities::WindowFunctions::hann(m, &w8data);
        }
        else if (type == TaperWindowType::Bartlett)
        {
            RTSeis::Utilities::WindowFunctions::bartlett(m, &w8data);
        }
        else if (type == TaperWindowType::Sine)
        {
            RTSeis::Utilities::WindowFunctions::sine(m, &w8data);
        }
        else if (type == TaperWindowType::Kaiser)
        {
            RTSeis::Utilities::WindowFunctions::kaiser(m, &w8data, beta);
        }
        else
        {
#ifndef NDEBUG
            assert(false);
#endif
            throw std::invalid_argument("Unsupported window");
        }
        pImpl->mWindowLength0 = m;
    }
    // Taper first (m+1)/2 points
    int mp12 = m/2;
    const double *w = pImpl->mW8.data();
    ippsMul_64f(w, x, y, mp12);
    // Copy the intermediate portion of the signal
    int ncopy = nx - mp12 - mp12; // Subtract out two window lengths
    if (ncopy > 0)
    {
        ippsCopy_64f(&x[mp12], &y[mp12], ncopy);
    } 
    // Taper last (m+1)/2 points 
    ippsMul_64f(&w[m-mp12], &x[nx-mp12], &y[nx-mp12], mp12);
}

/// Apply float
template<>
void Taper<float>::apply(const int nx, const float x[], float *yIn[])
{
    if (nx <= 0){return;}
    if (!isInitialized()){throw std::runtime_error("Taper never initialized");}
    float *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Deal with an edge case
    if (nx < 3)
    {
        y[0] = 0;
        if (nx == 2){y[1] = 0;}
        return;
    }
    // Compute taper length
    double pct = pImpl->mPercentage;
    int npct = static_cast<int> (std::round( (nx*pct)/100. )) + 1;
    int m = std::max(2, std::min(nx, npct));
    // Redesign the window?  If the parameters were (re)set then mWindowLength0
    // is -1.
    // Otherwise, if the same length signal is coming at us then the precision
    // of the module can't change so we can just use the old window.
    if (pImpl->mWindowLength0 != m)
    {
        pImpl->mW4.resize(m); // Prevent reallocations in a bit
        auto type = pImpl->mWindowType;
        auto beta = static_cast<float> (pImpl->mBeta);
        float *w4data = pImpl->mW4.data();
        if (type == TaperWindowType::Hamming)
        {
            RTSeis::Utilities::WindowFunctions::hamming(m, &w4data);
        }
        else if (type == TaperWindowType::Blackman)
        {
            RTSeis::Utilities::WindowFunctions::blackman(m, &w4data);
        }
        else if (type == TaperWindowType::Hann)
        {
            RTSeis::Utilities::WindowFunctions::hann(m, &w4data);
        }
        else if (type == TaperWindowType::Bartlett)
        {
            RTSeis::Utilities::WindowFunctions::bartlett(m, &w4data);
        }
        else if (type == TaperWindowType::Sine)
        {
            RTSeis::Utilities::WindowFunctions::sine(m, &w4data);
        }
        else if (type == TaperWindowType::Kaiser)
        {
            RTSeis::Utilities::WindowFunctions::kaiser(m, &w4data, beta);
        }
        else
        {
#ifndef NDEBUG
            assert(false);
#endif
            throw std::invalid_argument("Unsupported window");
        }
        pImpl->mWindowLength0 = m;
    }
    // Taper first (m+1)/2 points
    int mp12 = m/2;
    const float *w = pImpl->mW4.data();
    ippsMul_32f(w, x, y, mp12);
    // Copy the intermediate portion of the signal
    int ncopy = nx - mp12 - mp12; // Subtract out two window lengths
    if (ncopy > 0)
    {   
        ippsCopy_32f(&x[mp12], &y[mp12], ncopy);
    }   
    // Taper last (m+1)/2 points 
    ippsMul_32f(&w[m-mp12], &x[nx-mp12], &y[nx-mp12], mp12);
}

/// Initialized?
template<class T>
bool Taper<T>::isInitialized() const
{
    return pImpl->mInitialized;
}

///--------------------------------------------------------------------------///
///                             Template Instantiation                       ///
///--------------------------------------------------------------------------///
template class RTSeis::FilterImplementations::Taper<double>;
template class RTSeis::FilterImplementations::Taper<float>;
