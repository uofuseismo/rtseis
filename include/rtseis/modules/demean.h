#ifndef RTSEIS_MODULES_DEMEAN_H
#define RTSEIS_MODULES_DEMEAN_H 1
#include "rtseis/config.h"
#include "rtseis/enums.h"

#ifdef __cplusplus
namespace RTSeis
{
namespace Modules
{

class Demean
{
    public:
        Demean(void);
        ~Demean(void);
        int setParameters(const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE);
        int demean(const int nx, const double x[], double y[]);
        int demean(const int nx, const float  x[], float  y[]);
    private:
        int removeMean_(const int nx, const double x[], double y[]);
        int removeMean_(const int nx, const float x[], float y[]);
        int computeMean_(const int nx, const double x[]);
        int computeMean_(const int nx, const float x[]);
        enum rtseisPrecision_enum precision_ = RTSEIS_DOUBLE;
        double mean_ = 0;
        const bool lrt_ = false;
};


};
};
#endif

#endif
