#include <stdexcept>
#include <memory>
#include "rtseis/utilities/normalization/signBit.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"

using namespace RTSeis::Utilities::Normalization;

class SignBit::SignBitImpl
{
public:
    bool mInitialized = false;
};

/// C'tor
SignBit::SignBit() :
    pImpl(new SignBitImpl())
{
}

/// Copy c'tor
SignBit::SignBit(const SignBit &signBit)
{
    *this = signBit;
}

/// Move c'tor
[[maybe_unused]]
SignBit::SignBit(SignBit &&signBit) noexcept
{
    *this = signBit;
}

/// Destructor
SignBit::~SignBit() = default;

/// Copy assignment
SignBit& SignBit::operator=(const SignBit &signBit)
{
    if (&signBit == this){return *this;}
    pImpl = std::make_unique<SignBitImpl> (*signBit.pImpl);
    return *this;
}

/// Move assignment
SignBit& SignBit::operator=(SignBit &&signBit) noexcept
{
    if (&signBit == this){return *this;}
    pImpl = std::move(signBit.pImpl);
    return *this;
}

/// Clears the class
void SignBit::clear() noexcept
{
    pImpl->mInitialized = false;
}

/// Initializes the class
void SignBit::initialize() noexcept
{
    clear();
    pImpl->mInitialized = true;
}

/// Check if class is initialized
bool SignBit::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Sets the initial conditions
void SignBit::setInitialConditions()
{
    if (!isInitialized())
    {
        throw std::runtime_error("signBit not initialized");
    }
}

/// Resets the initial conditions
void SignBit::resetInitialConditions()
{
    if (!isInitialized())
    {
        throw std::runtime_error("signBit not initialized");
    }
}

/// Applies the sign bit normalization
template<typename U>
void SignBit::apply(const int nx, const U x[], U *yIn[])
{
    if (nx <= 0){return;}
    if (!isInitialized())
    {
        throw std::runtime_error("signBit not initialized");
    }
    U *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
#ifdef __INTEL_COMPILER
    #pragma forceinline
#endif
    RTSeis::Utilities::Math::VectorMath::copysign(nx, x, y);
}

/*
void SignBit::apply(const int nx, const float x[], float *yIn[])
{
    if (nx <= 0){return;} 
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "signBit not initialized");
    }
    float *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_RTE("%s", "x is NULL");}
        RTSEIS_THROW_RTE("%s", "y is NULL");
    }
#ifdef __INTEL_COMPILER
    #pragma forceinline
#endif
    RTSeis::Utilities::Math::VectorMath::copysign(nx, x, y);
}
*/

///--------------------------------------------------------------------------///
///                       Template Instantiation                             ///
///--------------------------------------------------------------------------///
template
void RTSeis::Utilities::Normalization::SignBit::apply(int nx, const double x[],
                                                      double *yIn[]);
template
void RTSeis::Utilities::Normalization::SignBit::apply(int nx, const float x[],
                                                      float *yIn[]);
