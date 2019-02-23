#include <cstdio>
#include <cstdlib>
#include <memory>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/utilities/normalization/signBit.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"

using namespace RTSeis::Utilities::Normalization;

class SignBit::SignBitImpl
{
    public:
        SignBitImpl(void){return;}
        ~SignBitImpl(void){linit = false;}
        SignBitImpl& operator=(const SignBitImpl &sb)
        {
            if (&sb == this){return *this;}
            linit = sb.linit;
            return *this;
        }
        bool linit = false;
};


SignBit::SignBit(void) :
    pSignBit_(new SignBitImpl())
{
    return;
}

SignBit::~SignBit(void)
{
    clear();
    return;
}

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

void SignBit::clear(void)
{
    pSignBit_->linit = false;
    return;
}

int SignBit::initialize(void)
{
    clear();
    pSignBit_->linit = true;
    return 0;
}

bool SignBit::isInitialized(void) const
{
    return pSignBit_->linit;
}

int SignBit::setInitialConditions(void)
{
    if (!isInitialized())
    {
        RTSEIS_WARNMSG("%s", "signBit not initialized");
    }
    return 0;
}

int SignBit::resetInitialConditions(void)
{
    if (!isInitialized())
    {
        RTSEIS_WARNMSG("%s", "signBit not initialized");
    }
    return 0;
}

int SignBit::apply(const int nx, const double x[], double y[])
{
    if (nx <= 0){return 0;}
    if (!isInitialized())
    {
        RTSEIS_WARNMSG("%s", "signBit not initialized");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1;
    }
    RTSeis::Utilities::Math::VectorMath::copysign(nx, x, y);
    return 0;
}

int SignBit::apply(const int nx, const float x[], float y[])
{
    if (nx <= 0){return 0;} 
    if (!isInitialized())
    {
        RTSEIS_WARNMSG("%s", "signBit not initialized");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1; 
    }
    RTSeis::Utilities::Math::VectorMath::copysign(nx, x, y);
    return 0;
}
