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
            const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        OneBitNormalizationParameters& operator=(const OneBitNormalizationParameters &parameters);
        OneBitNormalizationParameters(const OneBitNormalizationParameters &parameters);
        ~OneBitNormalizationParameters(void);
        virtual void clear(void);
        virtual bool isInitialized(void) const{return isInitialized_;}
        bool isRealTime(void) const{return isRealTime_;} 
        RTSeis::Precision getPrecision(void) const{return precision_;}
     private:
        /*!< Default precision. */
        const RTSeis::Precision defaultPrecision_ = RTSeis::Precision::DOUBLE;
        /*!< The precision of the module. */
        RTSeis::Precision precision_ = defaultPrecision_; 
        /*!< Flag indicating this module is for real-time. */
        bool isRealTime_ = false;
        /*!< Flag indicating the module is initialized. */
        bool isInitialized_ = false;
};

class OneBitNormalization : OneBitNormalizationParameters
{
    public:
        OneBitNormalization(void);
        OneBitNormalization(const OneBitNormalization &onebit);
        OneBitNormalization(const OneBitNormalizationParameters &parameters);
        OneBitNormalization& operator=(const OneBitNormalization &onebit);
        ~OneBitNormalization(void);
        int setInitialConditions(void);
        int apply(const int nx, const double x[], double y[]);
        int apply(const int nx, const float  x[], float  y[]);
        int resetInitialConditions(void);
        void clear(void) override;
        bool isInitialized(void) const override;
    private:
        /*!< The parameters. */ 
        OneBitNormalizationParameters parms_;
        /*!< Flag indicating the module is intialized. */
        bool isInitialized_ = false;
};


};
};

#endif
