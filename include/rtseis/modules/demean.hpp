#ifndef RTSEIS_MODULES_DEMEAN_HPP
#define RTSEIS_MODULES_DEMEAN_HPP 1
#include <memory>
#include "rtseis/config.h"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Modules
{

class DemeanParameters
{
    public:
        DemeanParameters(const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE);
        DemeanParameters& operator=(const DemeanParameters &parameters);
        DemeanParameters(const DemeanParameters &parameters);
        ~DemeanParameters(void);
        void clear(void);
        enum rtseisPrecision_enum getPrecision(void) const;
        bool getIsRealTime(void) const;
        bool isInitialized(void) const;
    private:
        const enum rtseisPrecision_enum defaultPrecision_ = RTSEIS_DOUBLE;
        enum rtseisPrecision_enum precision_ = defaultPrecision_;
        bool lrt_ = false;
        bool linit_ = true; // This module is always ready to roll
};

class Demean
{
    public:
        Demean(void);
        Demean(const Demean &demean);
        Demean(const DemeanParameters &parameters);
        Demean& operator=(const Demean &demean);
        ~Demean(void);
        int setParameters(const DemeanParameters &parameters);
        int demean(const int nx, const double x[], double y[]);
        int demean(const int nx, const float  x[], float  y[]);
        void clear(void);
    private:
        int removeMean_(const int nx, const double x[], double y[]);
        int removeMean_(const int nx, const float x[], float y[]);
        int computeMean_(const int nx, const double x[]);
        int computeMean_(const int nx, const float x[]);
        /*!< The parameters. */ 
        DemeanParameters parms_;
        /*!< The mean of the data. */
        double mean_ = 0;
        bool linit_ = true; // This module is always ready to roll
};


};
};

#endif
