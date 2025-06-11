#include <complex>
#include "rtseis/filterRepresentations/finiteImpulseResponse.hpp"
#include "rtseis/vector.hpp"

using namespace RTSeis::FilterRepresentations;

template<class T>
class FiniteImpulseResponse<T>::FiniteImpulseResponseImpl
{
public:
    class RTSeis::Vector<T> mFilterCoefficients;
};

/// Constructor
template<class T>
FiniteImpulseResponse<T>::FiniteImpulseResponse(
    const RTSeis::Vector<T> &filterCoefficients)
{
    constexpr T zero{0};
    if (filterCoefficients.empty())
    {
        throw std::invalid_argument("The filter coefficients are empty");
    }
    pImpl = std::make_unique<FiniteImpulseResponseImpl> ();
    pImpl->mFilterCoefficients = filterCoefficients;
    int nCoefficients = static_cast<int> (pImpl->mFilterCoefficients.size());
    for (int i = nCoefficients - 1; i >= 0; --i)
    {
        if (pImpl->mFilterCoefficients.at(i) != zero){break;}
        pImpl->mFilterCoefficients.pop_back(); 
    }
    if (pImpl->mFilterCoefficients.empty())
    {
        pImpl = nullptr;
        throw std::runtime_error("The filter coefficients are all zero");
    }
}

/// Copy constructor
template<class T>
FiniteImpulseResponse<T>::FiniteImpulseResponse(
    const FiniteImpulseResponse<T> &firFilter)
{
    *this = firFilter;
}

/// Move constructor
template<class T>
FiniteImpulseResponse<T>::FiniteImpulseResponse(
    FiniteImpulseResponse<T> &&firFilter) noexcept
{
    *this = std::move(firFilter);
}

/// Copy assignment
template<class T>
FiniteImpulseResponse<T> &
FiniteImpulseResponse<T>::operator=(const FiniteImpulseResponse<T> &firFilter)
{
    if (&firFilter == this){return *this;}
    pImpl = std::make_unique<FiniteImpulseResponseImpl> (*firFilter.pImpl);
    return *this;
}

/// Move assignment
template<class T>
FiniteImpulseResponse<T> &
FiniteImpulseResponse<T>::operator=(
    FiniteImpulseResponse<T> &&firFilter) noexcept
{
    if (&firFilter == this){return *this;}
    pImpl = std::move(firFilter.pImpl);
    return *this;
}

/// Filter coefficients
template<class T>
const RTSeis::Vector<T> &
FiniteImpulseResponse<T>::getFilterCoefficientsReference() const noexcept
{
    return *&pImpl->mFilterCoefficients;
}

/// Filter coefficients
template<class T>
RTSeis::Vector<T> 
FiniteImpulseResponse<T>::getFilterCoefficients() const noexcept
{
    return pImpl->mFilterCoefficients;
}

/// Filter order
template<class T>
int FiniteImpulseResponse<T>::getOrder() const noexcept
{
    return static_cast<int> (pImpl->mFilterCoefficients.size() - 1);
}

/// Destructor
template<class T>
FiniteImpulseResponse<T>::~FiniteImpulseResponse() = default;

template class RTSeis::FilterRepresentations::FiniteImpulseResponse<double>;
template class RTSeis::FilterRepresentations::FiniteImpulseResponse<float>;
template class RTSeis::FilterRepresentations::FiniteImpulseResponse<int>;
template class RTSeis::FilterRepresentations::FiniteImpulseResponse<std::complex<double>>;
template class RTSeis::FilterRepresentations::FiniteImpulseResponse<std::complex<float>>;

