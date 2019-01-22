#ifndef RTSEIS_MODULES_ONEBIT_HPP
#define RTSEIS_MODULES_ONEBIT_HPP 1
#include "rtseis/config.h"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Modules
{

class OneBitNormalizationParameters
{
    public:
        OneBitNormalizationParameters(
            const bool lrt = false,
            const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE);
        OneBitNormalizationParameters& operator=(const OneBitNormalizationParameters &parameters);
        OneBitNormalizationParameters(const OneBitNormalizationParameters &parameters);
        ~OneBitNormalizationParameters(void);
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

class OneBitNormalization
{
    public:
        OneBitNormalization(void);
        OneBitNormalization(const OneBitNormalization &onebit);
        OneBitNormalization(const OneBitNormalizationParameters &parameters);
        OneBitNormalization& operator=(const OneBitNormalization &onebit);
        ~OneBitNormalization(void);
        int setParameters(const OneBitNormalization &parameters);
        int apply(const int nx, const double x[], double y[]);
        int apply(const int nx, const float  x[], float  y[]);
        void clear(void);
    private:
        /*!< The parameters. */ 
        OneBitNormalizationParameters parms_;
        /*!< Flag indicating this is initialized. */
        bool linit_ = false;
};


};
};

#endif
