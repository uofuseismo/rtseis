#include <string>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <memory>
#include <ipps.h>
#include "rtseis/postProcessing/singleChannel/taper.hpp"
#include "rtseis/utilities/windowFunctions.hpp"

using namespace RTSeis::PostProcessing::SingleChannel;

class TaperParameters::TaperParametersImpl
{
public:
    /// Percentage of signal to apply taper to
    double pct = 5;
    /// Taper type.
    Type type = Type::HAMMING;
    /// Precision 
    //RTSeis::Precision precision = RTSeis::Precision::DOUBLE; 
    /// The processing mode can't be toggled
    //const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING;
};

template<class T>
class Taper<T>::TaperImpl
{
public:
    TaperParameters parms; 
    std::vector<double> w8;
    std::vector<float>  w4;
    int winLen0 =-1;
    bool linit = false;
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

//============================================================================//
//                                    Tapering                                //
//============================================================================//
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

/// C'tor from parameters
template<class T>
Taper<T>::Taper(const TaperParameters &parameters) :
    pImpl(std::make_unique<TaperImpl>())
{
    setParameters(parameters);
}

/// Destructor
template<class T>
Taper<T>::~Taper() = default;

template<class T>
void Taper<T>::clear() noexcept
{
    pImpl->parms.clear();
    pImpl->w8.clear();
    pImpl->w4.clear();
    pImpl->winLen0 =-1;
    pImpl->linit = false;
}

template<class T>
Taper<T>& Taper<T>::operator=(const Taper &taper)
{
    if (&taper == this){return *this;}
    pImpl = std::make_unique<TaperImpl> (*taper.pImpl);
    return *this;
}

template<class T>
Taper<T>& Taper<T>::operator=(Taper &&taper) noexcept
{
    if (&taper == this){return *this;}
    pImpl = std::move(taper.pImpl);
    return *this;
}

template<class T>
void Taper<T>::setParameters(const TaperParameters &parameters)
{
    clear(); // Sets winLen0 to -1
    if (!parameters.isValid())
    {
        throw std::invalid_argument("Taper parameters are invalid");
    }
    pImpl->parms = parameters;
    pImpl->linit = true;
}

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
    double pct = pImpl->parms.getPercentage();
    int npct = static_cast<int> (static_cast<double> (nx)*pct/100 + 0.5) + 1;
    int m = std::max(2, std::min(nx, npct));
    // Redesign the window?  If the parameters were (re)set then winLen0 is -1.
    // Otherwise, if the same length signal is coming at us then the precision
    // of the module can't change so we can just use the old window.
    if (pImpl->winLen0 != m)
    {
        pImpl->w8.resize(m); // Prevent reallocations in a bit
        TaperParameters::Type type = pImpl->parms.getTaperType();
        double *w8data = pImpl->w8.data();
        if (type == TaperParameters::Type::HAMMING)
        {
            RTSeis::Utilities::WindowFunctions::hamming(m, &w8data);
        }
        else if (type == TaperParameters::Type::BLACKMAN)
        {
            RTSeis::Utilities::WindowFunctions::blackman(m, &w8data);
        }
        else if (type == TaperParameters::Type::HANN)
        {
            RTSeis::Utilities::WindowFunctions::hann(m, &w8data);
        }
        else if (type == TaperParameters::Type::BARTLETT)
        {
            RTSeis::Utilities::WindowFunctions::bartlett(m, &w8data);
        }
        else if (type == TaperParameters::SINE)
        {
            RTSeis::Utilities::WindowFunctions::sine(m, &w8data);
        }
        else
        {
#ifndef NDEBUG
            assert(false);
#endif
            throw std::invalid_argument("Unsupported window");
        }
    }
    // Taper first (m+1)/2 points
    int mp12 = m/2;
    const double *w = pImpl->w8.data();
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
    double pct = pImpl->parms.getPercentage();
    int npct = static_cast<int> (static_cast<double> (nx)*pct/100 + 0.5) + 1;
    int m = std::max(2, std::min(nx, npct));
    // Redesign the window?  If the parameters were (re)set then winLen0 is -1.
    // Otherwise, if the same length signal is coming at us then the precision
    // of the module can't change so we can just use the old window.
    if (pImpl->winLen0 != m)
    {
        pImpl->w4.resize(m); // Prevent reallocations in a bit
        TaperParameters::Type type = pImpl->parms.getTaperType();
        float *w4data = pImpl->w4.data();
        if (type == TaperParameters::Type::HAMMING)
        {
            RTSeis::Utilities::WindowFunctions::hamming(m, &w4data);
        }
        else if (type == TaperParameters::Type::BLACKMAN)
        {
            RTSeis::Utilities::WindowFunctions::blackman(m, &w4data);
        }
        else if (type == TaperParameters::Type::HANN)
        {
            RTSeis::Utilities::WindowFunctions::hann(m, &w4data);
        }
        else if (type == TaperParameters::Type::BARTLETT)
        {
            RTSeis::Utilities::WindowFunctions::bartlett(m, &w4data);
        }
        else if (type == TaperParameters::SINE)
        {
            RTSeis::Utilities::WindowFunctions::sine(m, &w4data);
        }
        else
        {
#ifndef NDEBUG
            assert(false);
#endif
            throw std::invalid_argument("Unsupported window");
        }
    }
    // Taper first (m+1)/2 points
    int mp12 = m/2;
    const float *w = pImpl->w4.data();
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

template<class T>
bool Taper<T>::isInitialized() const
{
    return pImpl->linit;
}

///--------------------------------------------------------------------------///
///                             Template Instantiation                       ///
///--------------------------------------------------------------------------///
template class RTSeis::PostProcessing::SingleChannel::Taper<double>;
template class RTSeis::PostProcessing::SingleChannel::Taper<float>;
