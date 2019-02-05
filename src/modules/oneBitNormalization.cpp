#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/oneBitNormalization.hpp"

using namespace RTSeis::Modules;

class OneBitNormalizationParameters::OneBitVars
{
    public:
        OneBitVars(void){return;}
        OneBitVars(const RTSeis::ProcessingMode mode,
                   const RTSeis::Precision prec) :
                   precision_(prec),
                   processingMode_(mode),
                   isValid_(true)
        {
            return;
        }
        OneBitVars& operator=(const OneBitVars &pParms)
        {
            if (&pParms == this){return *this;}
            precision_ = pParms.precision_;
            processingMode_ = pParms.processingMode_;
            isValid_ = pParms.isValid_;
            return *this;
        }
        OneBitVars(const OneBitVars &pParms)
        {
            *this = pParms;
            return;
        }
        ~OneBitVars(void)
        {
            clear();
        }  
        void clear(void)
        {
            processingMode_ = RTSeis::ProcessingMode::POST_PROCESSING;
            precision_ = defaultPrecision_;
            isValid_ = true;
            return;
        }
        bool isValid(void) const
        {
            return isValid_;
        }
        RTSeis::ProcessingMode getProcessingMode(void) const
        {
            return processingMode_;
        }
        int setProcessingMode(const RTSeis::ProcessingMode mode)
        {
            processingMode_ = mode;
            validate();
            return 0;
        }
        int setPrecision(const RTSeis::Precision prec)
        {
            precision_ = prec;
            validate();
            return 0;
        } 
        RTSeis::Precision getPrecision(void) const
        {
            return precision_;
        }
        void validate(void)
        {
            isValid_ = false;
            if (getPrecision() != RTSeis::Precision::DOUBLE &&
                getPrecision() != RTSeis::Precision::FLOAT){return;}
            isValid_ = true;
            return;
        }
        const RTSeis::Precision defaultPrecision_ = RTSeis::Precision::DOUBLE;
        RTSeis::Precision precision_ = defaultPrecision_;
        RTSeis::ProcessingMode processingMode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        bool isValid_ = true;
};

OneBitNormalizationParameters::OneBitNormalizationParameters(
    const RTSeis::ProcessingMode mode,
    const RTSeis::Precision prec) :
//    precision_(prec),
//    processingMode_(mode),
//    isValid_(true),
    pImpl_(new OneBitVars(mode, prec))
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
    pImpl_ = std::unique_ptr<OneBitVars> (new OneBitVars(*parameters.pImpl_));
/*
    precision_ = parameters.precision_;
    processingMode_ = parameters.processingMode_;
    isValid_ = parameters.isValid_;
*/
    return *this;
}

OneBitNormalizationParameters::~OneBitNormalizationParameters(void)
{
    clear();
    return;
}

void OneBitNormalizationParameters::clear(void)
{
    pImpl_->clear();
/*
    precision_ = defaultPrecision_;
    processingMode_ = RTSeis::ProcessingMode::POST_PROCESSING;
    isValid_ = true; // This is still a valid processing class 
*/
    return;
}

bool OneBitNormalizationParameters::isValid(void) const
{
    return pImpl_->isValid();
    //return isValid_;
}

void OneBitNormalizationParameters::setProcessingMode(
    const RTSeis::ProcessingMode mode)
{
    pImpl_->setProcessingMode(mode);
/*
    processingMode_ = mode;
    validate_();
*/
    return;
} 

RTSeis::ProcessingMode
OneBitNormalizationParameters::getProcessingMode(void) const
{
    return pImpl_->getProcessingMode();
    //return processingMode_;
}

RTSeis::Precision OneBitNormalizationParameters::getPrecision(void) const
{
    return pImpl_->getPrecision();
    //return precision_;
}

/*
void OneBitNormalizationParameters::validate_(void)
{
    isValid_ = false;
    if (getPrecision() != RTSeis::Precision::DOUBLE &&
        getPrecision() != RTSeis::Precision::FLOAT){return;}
    isValid_ = true;
    return;
}
*/
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
