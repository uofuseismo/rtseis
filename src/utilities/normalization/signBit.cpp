#include <cstdio>
#include <cstdlib>
#include <memory>
#define RTSEIS_LOGGING 1
#include "rtseis/private/throw.hpp"
#include "rtseis/log.h"
#include "rtseis/utilities/normalization/signBit.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"

using namespace RTSeis::Utilities::Normalization;

class SignBit::SignBitImpl
{
    public:
        SignBitImpl() = default;
        ~SignBitImpl() = default;
        SignBitImpl& operator=(const SignBitImpl &sb)
        {
            if (&sb == this){return *this;}
            linit = sb.linit;
            return *this;
        }
        bool linit = false;
};


SignBit::SignBit() :
    pSignBit_(new SignBitImpl())
{
}

SignBit::~SignBit() = default;

SignBit::SignBit(const SignBit &signBit)
{
    *this = signBit;
    return;
}

SignBit& SignBit::operator=(const SignBit &signBit)
{
    if (&signBit == this){return *this;}
    pSignBit_ = std::unique_ptr<SignBitImpl> (new SignBitImpl(*signBit.pSignBit_));
    pSignBit_->linit = signBit.pSignBit_->linit;
    return *this;
}

void SignBit::clear() noexcept
{
    pSignBit_->linit = false;
}

void SignBit::initialize() noexcept
{
    clear();
    pSignBit_->linit = true;
}

bool SignBit::isInitialized() const noexcept
{
    return pSignBit_->linit;
}

void SignBit::setInitialConditions()
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "signBit not initialized");
    }
}

void SignBit::resetInitialConditions()
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "signBit not initialized");
    }
}

void SignBit::apply(const int nx, const double x[], double y[])
{
    if (nx <= 0){return;}
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "signBit not initialized");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
#ifdef __INTEL_COMPILER
    #pragma forceinline
#endif
    RTSeis::Utilities::Math::VectorMath::copysign(nx, x, y);
}

void SignBit::apply(const int nx, const float x[], float y[])
{
    if (nx <= 0){return;} 
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "signBit not initialized");
    }
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
