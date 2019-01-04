#ifndef RTSEIS_MODULES_DETREND_H
#define RTSEIS_MODULES_DETREND_H 1
#include "rtseis/config.h"
#include "rtseis/enums.h"

#ifdef __cplusplus
namespace RTSeis
{
namespace Modules
{

class Detrend
{
    public:
        Detrend(void);
        ~Detrend(void);
        int setParameters(const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE);
        int detrend(const int nx, const double x[], double y[]);
        int detrend(const int nx, const float  x[], float  y[]);
    private:
        int removeTrend_(const int nx, const double x[], double y[]);
        int removeTrend_(const int nx, const float x[],  float y[]);
        int computeLinearRegressionCoeffs_(const int length, const double x[]);
        int computeLinearRegressionCoeffs_(const int length, const float x[]); 
        enum rtseisPrecision_enum precision_ = RTSEIS_DOUBLE;
        double b0_ = 0;
        double b1_ = 0;
        const bool lrt_ = false;
};


};
};
#endif

#endif
