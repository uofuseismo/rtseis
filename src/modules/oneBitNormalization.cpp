#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/oneBitNormalization.hpp"

using namespace RTSeis::Modules;

OneBitNormalizationParameters::OneBitNormalizationParameters(
    const bool lrt,
    const RTSeis::Precision prec) :
    precision_(prec),
    isRealTime_(lrt),
    isValid_(true)
{
    return;
}

OneBitNormalizationParameters::OneBitNormalizationParameters(
    const OneBitNormalizationParameters &parameters)
{
    *this = parameters;
    return;
}

OneBitNormalizationParameters&
OneBitNormalizationParameters::operator=(
    const OneBitNormalizationParameters &parameters)
{
    if (&parameters == this){return *this;}
    clear();
    precision_ = parameters.precision_;
    isRealTime_ = parameters.isRealTime_;
    isValid_ = parameters.isValid_;
    return *this;
}

OneBitNormalizationParameters::~OneBitNormalizationParameters(void)
{
    clear();
    return;
}

void OneBitNormalizationParameters::clear(void)
{
    precision_ = defaultPrecision_;
    isRealTime_ = false;
    isValid_ = true; // This is still a valid processing class 
    return;
}

bool OneBitNormalizationParameters::isValid(void) const
{
    return isValid_;
}

void OneBitNormalizationParameters::setRealTime(const bool lrt)
{
    isRealTime_ = lrt;
    validate_();
    return;
} 

bool OneBitNormalizationParameters::getRealTime(void) const
{
    return isRealTime_;
}

RTSeis::Precision OneBitNormalizationParameters::getPrecision(void) const
{
    return precision_;
}

void OneBitNormalizationParameters::validate_(void)
{
    isValid_ = false;
    if (getPrecision() != RTSeis::Precision::DOUBLE &&
        getPrecision() != RTSeis::Precision::FLOAT){return;}
    isValid_ = true;
    return;
}
//============================================================================//
//                                   End Parameters                           //
//============================================================================//

OneBitNormalization::OneBitNormalization(void) :
    isInitialized_(false)
{
    clear();
    return;
}

OneBitNormalization::OneBitNormalization(const OneBitNormalization &onebit)
{
    *this = onebit;
    return;
}

OneBitNormalization::OneBitNormalization(
     const OneBitNormalizationParameters &parameters)
{
    clear();
    if (!parameters.isValid())
    {
        RTSEIS_ERRMSG("%s", "Input parameters are not yet valid");
        return;
    }
    parms_ = parameters;
    isInitialized_ = true;
    return;
}
 
OneBitNormalization&
OneBitNormalization::operator=(const OneBitNormalization &onebit)
{
    if (&onebit == this){return *this;}
    clear();
    parms_ = onebit.parms_;
    isInitialized_ = onebit.isInitialized_;
    return *this;
}

int OneBitNormalization::initialize(const OneBitNormalizationParameters &parameters)
{
    clear();
    if (!parameters.isValid())
    {
        RTSEIS_ERRMSG("%s", "Parameters are not valid");
        return -1;
    }
    *this = OneBitNormalization(parameters);
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize class");
        return -1;
    }
    return 0;
}

int OneBitNormalization::setInitialConditions(void)
{
    return 0;
}

int OneBitNormalization::resetInitialConditions(void)
{
    return 0;
}

OneBitNormalization::~OneBitNormalization(void)
{
    clear();
    return;
}

void OneBitNormalization::clear(void)
{
    parms_.clear();
    isInitialized_ = false;
    return;
}

int OneBitNormalization::apply(const int nx, const double x[], double y[])
{
    if (nx <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class is not initialied");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is null");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is null");}
        return -1;
    }
    #pragma omp simd
    for (int i=0; i<nx; i++)
    {
        y[i] = std::copysign(1.0, x[i]);
    }
    return 0;
}

int OneBitNormalization::apply(const int nx, const float x[], float y[])
{
    if (nx <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class is not initialied");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is null");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is null");}
        return -1; 
    }
    #pragma omp simd
    for (int i=0; i<nx; i++)
    {
        y[i] = std::copysign(1.0f, x[i]);
    }
    return 0;
}

bool OneBitNormalization::isInitialized(void) const
{
    return isInitialized_;
}
