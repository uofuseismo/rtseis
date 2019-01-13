#ifndef RTSEIS_MODULES_FIR_HPP
#define RTSEIS_MODULES_FIR_HPP 1
#include "rtseis/config.h"
#include "rtseis/enums.h"
#include "rtseis/ippsHelper.h"
#include "rtseis/utils/design.hpp"

namespace RTSeis
{
namespace Modules
{

class FIRParameters
{
    public:
        enum class Window
        {
            /*!
             * @brief Defines the window in the FIR window-based
                      design.
             * @ingroup rtseis_modules_fir
             */
            HAMMING  = 0,     /*!< Hamming window. */
            BARTLETT = 1,     /*!< Bartlett window. */
            HANN = 2,         /*!< Hann window. */
            BLACKMAN_OPT = 3, /*!< Optimal Blackman window. */
            CUSTOM = 100      /*!< Indicates a custom filter. */
        };
        enum class Bandtype
        {
            /*!
             * @brief Defines the filter band type.
             * @ingroup rtseis_modules_fir
             */
            LOWPASS = 0,      /*!< Lowpass filter. */
            HIGHPASS = 1,     /*!< Highpass filter. */
            BANDPASS = 2,     /*!< Bandpass filter. */
            BANDSTOP = 3,     /*!< Bandstop filter. */
            CUSTOM = 100      /*!< Indicates a custom filter. */
        };
    public:
        FIRParameters(void);
        FIRParameters& operator=(const FIRParameters &parameters)
        {
            ba_ = parameters.ba_;
            window_ = parameters.window_;
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
        int designFromWindow(const int order, const double fs,
                             const double f0,
                             const Bandtype btype,
                             const Window window = Window::HAMMING);
        int designFromWindow(const int order, const double fs,
                             const double f0, const double f1,
                             const Bandtype btype,
                             const Window window = Window::HAMMING);
        FIRParameters(const FIRParameters &parameters);
        ~FIRParameters(void);
        void clear(void);
        void resetFilterDesign(void);
    private:
        /*!< Default analog filter prototype. */
        const Window defaultWindow_ = Window::HAMMING;
        /*!< Default filter bandtype. */
        const Bandtype defaultBandtype_ = Bandtype::LOWPASS;
        /*!< By default the module operates in double precision. */
        const enum rtseisPrecision_enum defaultPrecision_ = RTSEIS_DOUBLE;
        /*!< The default sampling frequency (Hz). */
        const double fsDefault_ = 1;
        /*!< The default filter order. */
        const int defaultOrder_ = 30;
        /*!< By default this is for post-processing. */
        const bool lrtDefault_ = false;
        /*!< The filter taps (coefficients) defining the FIR filter. */
        BA ba_;
        /*!< The filter prototype. */
        Window window_ = defaultWindow_;
        /*!< The filter bandtype. */
        Bandtype bandtype_ = defaultBandtype_;
        /*! The critical frequencies. */
        double w0_[2];
        /*!< The sampling frequency (Hz). */
        double fs_ = fsDefault_;
        /*!< The filter order. */
        int order_ = defaultOrder_;
        /*!< Preciison of module. */
        enum rtseisPrecision_enum precision_ = defaultPrecision_;
        /*!< Flag indicating this is for real-time. */
        bool lrt_ = lrtDefault_;
        /*!< Flag indicating a filter is set. */
        bool lhaveFilter_ = false;
};

}; /* End Modules */
}; /* End RTSeis */

#endif
