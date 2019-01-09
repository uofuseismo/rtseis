#ifndef RTSEIS_MODULES_DETREND_H
#define RTSEIS_MODULES_DETREND_H 1
#include "rtseis/config.h"
#include "rtseis/enums.h"

#ifdef __cplusplus
namespace RTSeis
{
namespace Modules
{

class DetrendParameters
{
    public:
        DetrendParameters(void);
        DetrendParameters& operator=(const DetrendParameters &parameters)
        {
            precision_ = parameters.precision_;
            lrt_ = parameters.lrt_;
            return *this;
        }
        DetrendParameters(const DetrendParameters &parameters);
        ~DetrendParameters(void);
        int setParameters(const DetrendParameters &parameters)
        {
            *this = parameters;
            return 0;
        }
        void setDefaults(void);
        int setPrecision(const enum rtseisPrecision_enum precision);
        enum rtseisPrecision_enum getPrecision(void) const;
        bool getIsRealTime(void) const;
    private:
        const enum rtseisPrecision_enum defaultPrecision_ = RTSEIS_DOUBLE;
        enum rtseisPrecision_enum precision_ = defaultPrecision_;
        bool lrt_ = false;
};

class Detrend : public DetrendParameters
{
/*
    public:
        class Parameters
        {
            public:
                Parameters(void);
                Parameters& operator=(const DetrendParameters &parameters)
                {
                    precision_ = parameters.precision_;
                    return *this;
                }
                Parameters(const Parameters &parameters);
                ~Parameters(void);
                void setDefaults(void);
                int setPrecision(const enum rtseisPrecision_enum precision);
                enum rtseisPrecision_enum getPrecision(void) const;
                bool getIsRealTime(void) const;
            private:
                const enum rtseisPrecision_enum defaultPrecision_ = RTSEIS_DOUBLE;
                enum rtseisPrecision_enum precision_ = defaultPrecision_;
                bool lrt_ = false;
        };
*/
    public:
        Detrend(void);
        Detrend(const DetrendParameters &parameters)
        {
            parms_ = parameters;
        }
        ~Detrend(void);
        int detrend(const int nx, const double x[], double y[]);
        int detrend(const int nx, const float  x[], float  y[]);
    private:
        DetrendParameters parms_;
        int removeTrend_(const int nx, const double x[], double y[]);
        int removeTrend_(const int nx, const float x[],  float y[]);
        int computeLinearRegressionCoeffs_(const int length, const double x[]);
        int computeLinearRegressionCoeffs_(const int length, const float x[]); 
        double b0_ = 0;
        double b1_ = 0;
};


};
};
#endif

#endif
