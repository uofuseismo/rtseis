#ifndef RTSEIS_MODULES_DETREND_HPP
#define RTSEIS_MODULES_DETREND_HPP 1
#include "rtseis/config.h"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Modules
{
 
class DetrendParameters
{
    public:
        DetrendParameters(const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE);
        DetrendParameters& operator=(const DetrendParameters &parameters);
        DetrendParameters(const DetrendParameters &parameters);
        ~DetrendParameters(void);
        void clear(void);
        //void setDefaults(void);
        //int setPrecision(const enum rtseisPrecision_enum precision);
        enum rtseisPrecision_enum getPrecision(void) const;
        bool getIsRealTime(void) const;
        bool isInitialized(void) const;
    private:
        const enum rtseisPrecision_enum defaultPrecision_ = RTSEIS_DOUBLE;
        enum rtseisPrecision_enum precision_ = defaultPrecision_;
        bool lrt_ = false;
        bool linit_ = true; // This module is always ready to roll
};

class Detrend
{
    public:
        Detrend(void);
        Detrend(const Detrend &detrend);
        Detrend(const DetrendParameters &parameters);
        Detrend &operator=(const Detrend &detrend);
        ~Detrend(void);
        int setParameters(const DetrendParameters &parameters);
        int detrend(const int nx, const double x[], double y[]);
        int detrend(const int nx, const float  x[], float  y[]);
        void clear(void);
    private:
        int removeTrend_(const int nx, const double x[], double y[]);
        int removeTrend_(const int nx, const float x[],  float y[]);
        int computeLinearRegressionCoeffs_(const int length, const double x[]);
        int computeLinearRegressionCoeffs_(const int length, const float x[]); 
        /*!< The detrend pararameters. */
        DetrendParameters parms_;
        /*!< The y-intercept */
        double b0_ = 0;
        /*!< The slope.  */
        double b1_ = 0;
        bool linit_ = true; // This module is always ready to roll
};


};
};

#endif
