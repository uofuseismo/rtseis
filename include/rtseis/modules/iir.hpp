#ifndef RTSEIS_MODULES_IIR_HPP
#define RTSEIS_MODULES_IIR_HPP 1
#include "rtseis/config.h"
#include "rtseis/enums.h"
#include "rtseis/ippsHelper.h"

namespace RTSeis
{
namespace Utilities::FilterRepresentations
{
class BA;
class ZPK;
};
namespace Modules
{

class IIRParameters
{
    public:
        enum class Prototype
        {
            /*!
             * @brief Defines the filter prototype.
             * @ingroup rtseis_modules_iir
             */
            BUTTERWORTH = 0,  /*!< Butterworth filter. */
            BESSEL = 1,       /*!< Bessel filter. */
            CHEBYSHEV1 = 2,   /*!< Chebyshev I filter. */
            CHEBYSHEV2 = 3,   /*!< Chebyshev II filter. */
            CUSTOM = 100      /*!< Indicates a custom filter. */
        };
        enum class Bandtype
        {
            /*!
             * @brief Defines the filter band type.
             * @ingroup rtseis_modules_iir
             */
            LOWPASS = 0,      /*!< Lowpass filter. */
            HIGHPASS = 1,     /*!< Highpass filter. */
            BANDPASS = 2,     /*!< Bandpass filter. */
            BANDSTOP = 3,     /*!< Bandstop filter. */
            CUSTOM = 100      /*!< Indicates a custom filter. */
        };
    public:
        IIRParameters(void);
        IIRParameters& operator=(const IIRParameters &parameters)
        {
            ba_ = parameters.ba_;
            prototype_ = parameters.prototype_;
            bandtype_ = parameters.bandtype_;
            w0_[0] = parameters.w0_[0];
            w0_[1] = parameters.w0_[1];
            fs_ = parameters.fs_;
            order_ = parameters.order_;
            precision_ = parameters.precision_;
            lrt_ = parameters.lrt_;
            lhaveFilter_ = parameters.lhaveFilter_;
            return *this;
        }
        IIRParameters(const IIRParameters &parameters);
        ~IIRParameters(void);
        int designFromAnalogPrototype(const int order, const double fs,
                                      const double f0,
                                      const Bandtype btype,
                                      const Prototype ftype,
                                      const double ripple = 5);
        int designFromAnalogPrototype(const int order, const double fs,
                                      const double f0, const double f1,
                                      const Bandtype btype,
                                      const Prototype ftype,
                                      const double ripple = 5);
        int setCustomFilter(const FilterRepresentations::ZPK &zpk);
        int setCustomFilter(const FilterRepresentations::BA &ba);
        //BA getFilter(void) const;
        void clear(void); 
        void resetFilterDesign(void);
        // Functions you probably shouldn't use
        int setPrecision(const enum rtseisPrecision_enum precision);
        enum rtseisPrecision_enum getPrecision(void) const;
        void setIsRealTime(const bool lrt);
        bool getIsRealTime(void) const;
        int setCornerFrequency(const double f0, const double fs);
        int setCornerFrequencies(const double f0, const double f1,
                                 const double fs);
    private:
        /*!< Default analog filter prototype. */
        const Prototype defaultPrototype_ = Prototype::BUTTERWORTH;
        /*!< Default filter bandtype. */
        const Bandtype defaultBandtype_ = Bandtype::LOWPASS;
        /*!< By default the module operates in double precision. */
        const enum rtseisPrecision_enum defaultPrecision_ = RTSEIS_DOUBLE;
        /*!< The default sampling frequency. */
        const double fsDefault_ = 1;
        /*!< The default filter order. */
        const int defaultOrder_ = 2;
        /*!< By default this is for post-processing. */
        const bool lrtDefault_ = false;
        /*!< The filter design must result in a digital filter. */
        const bool lanalog_ = false;
        /*!< The transfer function defining the filter. */
        Utilties::FilterRepresentations::BA ba_;
        /*!< The filter prototype. */
        Prototype prototype_ = defaultPrototype_;
        /*!< The filter bandtype. */
        Bandtype bandtype_ = defaultBandtype_;
        /*! The critical frequencies (Hz). */
        double w0_[2]; 
        /*!< The sampling frequench (Hz). */
        double fs_ = fsDefault_;
        /*!< The filter order. */
        int order_ = defaultOrder_;
        /*!< Preciison of module. */
        enum rtseisPrecision_enum precision_ = defaultPrecision_;
        /*!< Flag indicating this is for real-time. */
        bool lrt_ = lrtDefault_;
        /*!< Flag indicating a filter is initialized. */
        bool lhaveFilter_ = false;
};

class IIR
{
    public:
        IIR(void);
        IIR(const IIRParameters &parameters);
        int setParameters(const IIRParameters &parameters);
        ~IIR(void);
        int compute(const int nx, const double x[], double y[]);
        int compute(const int nx, const float  x[], float  y[]);
    private:
        IIRParameters parms_;
        struct ippsIIRFilter_struct ippsIIR_;
};


}; /* End Modules */
}; /* End RTSeis */

#endif
