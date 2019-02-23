#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/postProcessing/singleChannel/demean.hpp"

using namespace RTSeis::PostProcessing::SingleChannel;

class DemeanParameters::DemeanParms
{
    public:
        void clear(void)
        {
            precision_ = RTSeis::Precision::DOUBLE;
            mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
            linit_ = true;
            return;
        }
        RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
        RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        bool linit_ = true; // This module is always ready to roll
};

class Demean::DemeanImpl
{
    public:
        /// Default constructor
        DemeanImpl(void)
        {
            return;
        }
        /// Copy constructor
        DemeanImpl(const DemeanImpl &demean)
        {
            *this = demean;
            return;
        }
        /// Deep copy operator
        DemeanImpl& operator=(const DemeanImpl &demean)
        {
            if (&demean == this){return *this;}
            parms_ = demean.parms_;
            mean_ = demean.mean_;
            linit_ = demean.linit_;
            return *this;
        }
        /// Destructor
        ~DemeanImpl(void)
        {
            clear();
            return;
        }
        /// Resets the module
        void clear(void)
        {
            parms_.clear();
            mean_ = 0;
            linit_ = true; // Module always ready to roll 
            return;
        }
        /// Sets the parameters
        void setParameters(const DemeanParameters &parameters)
        {
            parms_ = parameters;
            return;
        }
        /// Removes mean from data
        int apply(const int nx, const double x[], double y[])
        {
            mean_ = 0;
            if (nx <= 0){return 0;}
            ippsMean_64f(x, nx, &mean_); // Compute mean of input 
            ippsSubC_64f(x, mean_, y, nx); // y - mean(x)
            return 0;
        }
        /// Removes mean from data
        int apply(const int nx, const float x[], float y[])
        {
            mean_ = 0;
            if (nx <= 0){return 0;} 
            float mean32;
            ippsMean_32f(x, nx, &mean32, ippAlgHintAccurate);
            mean_ = static_cast<double> (mean32);
            ippsSubC_32f(x, mean32, y, nx); // y - mean(x)
            return 0;
        }
        /// Utility to set the mean
        void setMean(const double mean)
        {
            mean_ = mean;
        }
    private: 
        DemeanParameters parms_;
        double mean_ = 0;
        bool linit_ = true; // This module is always ready to roll
};

DemeanParameters::DemeanParameters(const RTSeis::Precision precision) :
    pDemeanParmsImpl_(new DemeanParms())
{
    pDemeanParmsImpl_->precision_ = precision;
    return;
}

DemeanParameters::DemeanParameters(const DemeanParameters &parameters)
{
    *this = parameters;
    return;
}

DemeanParameters&
    DemeanParameters::operator=(const DemeanParameters &parameters)
{
    if (&parameters == this){return *this;}
    if (pDemeanParmsImpl_){pDemeanParmsImpl_->clear();}
    pDemeanParmsImpl_ = std::unique_ptr<DemeanParms>
                        (new DemeanParms(*parameters.pDemeanParmsImpl_));
    return *this;
}

DemeanParameters::~DemeanParameters(void)
{
    clear();
    return;
}

void DemeanParameters::clear(void)
{
    pDemeanParmsImpl_->clear();
    return;
}

RTSeis::Precision DemeanParameters::getPrecision(void) const
{
    return pDemeanParmsImpl_->precision_;
}

RTSeis::ProcessingMode DemeanParameters::getProcessingMode(void) const
{
    return pDemeanParmsImpl_->mode_;
}

bool DemeanParameters::isInitialized(void) const
{
    return pDemeanParmsImpl_->linit_;
}
//============================================================================//

Demean::Demean(void) :
    pDemean_(new DemeanImpl())
{
    return;
}

Demean::Demean(const Demean &demean)
{
    *this = demean;
    return;
}

Demean::Demean(const DemeanParameters &parameters) :
    pDemean_(new DemeanImpl())
{
    setParameters(parameters);
    return;
}

Demean& Demean::operator=(const Demean &demean)
{
    if (&demean == this){return *this;}
    if (pDemean_){pDemean_->clear();}
    pDemean_ = std::unique_ptr<DemeanImpl> (new DemeanImpl(*demean.pDemean_));
    return *this;
}

Demean::~Demean(void)
{
    clear();
    return;
}

void Demean::clear(void)
{
    pDemean_->clear();
    return;
}

void Demean::setParameters(const DemeanParameters &parameters)
{
    clear();
    if (!parameters.isInitialized())
    {
        RTSEIS_THROW_IA("%s", "Parameters not correctly initialized");
    }
    pDemean_->setParameters(parameters);
    return;
}

void Demean::apply(const int nx, const double x[], double y[])
{
    pDemean_->setMean(0);
    if (nx <= 0){return;}
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is null");}
        if (y == nullptr){RTSEIS_THROW_IA("%s", "y is null");}
        throw std::invalid_argument("Invalid input");
    }
    pDemean_->apply(nx, x, y);
    return;
}

void Demean::apply(const int nx, const float x[], float y[])
{
    pDemean_->setMean(0);
    if (nx <= 0){return;} // Nothing to do
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is null");}
        if (y == nullptr){RTSEIS_THROW_IA("%s", "y is null");}
        throw std::invalid_argument("Invalid input");
    }
    pDemean_->apply(nx, x, y);
    return;
}

