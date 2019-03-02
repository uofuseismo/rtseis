#ifndef RTSEIS_POSTPROCESSING_SC_TAPER
#define RTSEIS_POSTPROCESSING_SC_TAPER 1
#include <memory>
#include <exception>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace PostProcessing
{
namespace SingleChannel
{

/*!
 * @defgroup rtseis_postprocessing_sc_taper Taper
 * @brief Utilities for tapering a signal.  This can be useful for mitigating
 *        wraparound artificats prior to application of Fourier methods.
 * @ingroup rtseis_postprocessing_sc
 */

/*!
 * @class TaperParameters taper.hpp "include/rtseis/processing/singleChannel/taper.hpp"
 * @brief Defines the parameters for tapering module.
 * @ingroup rtseis_postprocessing_sc_taper
 * @copyright Ben Baker distributed under the MIT license.
 */
class TaperParameters : public std::exception
{
    public:
        /*!
         * @brief Defines the supported tapers types that can be applied
         *        to a signal.
         */
        enum Type
        {
            HAMMING,   /*!< Apply Hamming window to ends of signal. */
            BARTLETT,  /*!< Apply Bartlett window to ends of signal. */
            HANN,      /*!< Apply Hann window to ends of signal. */
            BLACKMAN,  /*!< Apply Blackman window to ends of signal. */
            SINE       /*!< Apply a sine window to the ends of the signal. */
        };

    public:
        /*!
         * @brief Default constructor.
         * @param[in] pct   The percentage of the signal to which the taper
         *                  will be applied.  For example, if 5 percent then
         *                  the initial 2.5 percent and final 2.5 percent of 
         *                  the signal will be tapered with the given window.
         *                  This must be in the range [0,100].
         * @param[in] type  The window type.
         * @param[in] precision  The precision which the taper will be applied.
         * @throws std::invalid_argument if parameters are incorrect.
         */ 
        TaperParameters(const double pct = 5,
                        const Type type = Type::HAMMING,
                        const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        /*!
         * @brief Copy constructor.
         * @param[in] parms  Taper parameters class from which to initialize
         *                   this class.
         */
        TaperParameters(const TaperParameters &parms);
        /*!
         * @brief Move constructor.
         * @param[in,out] parms  Taper parameters class to move to this class.
         *                       On exit this class will no longer be usable.
         */
        TaperParameters(TaperParameters &&parms);
        /*!
         * @brief Assignment operator.
         * @param[in] parms  Taper class to copy.
         * @result A deep copy of the taper parameters class.
         */
        TaperParameters& operator=(const TaperParameters &parms);
        /*!
         * @brief Move operator.
         * @param[in,out] parms  Taper parameters to move.  On exit,
         *                       this class is no longer usable.
         */
        TaperParameters& operator=(TaperParameters &&parms);
        /*!
         * @brief Default destructor.
         */
        ~TaperParameters(void);
        /*!
         * @brief Sets the taper type.
         * @param[in] type   The taper type.
         */
        void setTaperType(const Type type);
        /*!
         * @brief Gets the taper type.
         * @result The taper type to apply to the signal.
         */
        Type getTaperType(void) const;
        /*!
         * @brief Sets the precision which the taper should be applied.
         * @param[in] precision  The precision which the taper should
         *                       be applied.
         */
        void setPrecision(RTSeis::Precision precision);
        /*!
         * @brief Gets the precision which the taper should be applied.
         * @result The precision which the taper should be applied.
         */
        RTSeis::Precision getPrecision(void) const;
        /*!
         * @brief Sets the percentage of the signal to taper.
         * @param[in] pct   The percentage of the signal to which the taper
         *                  will be applied.  For example, if 5 percent then
         *                  the initial 2.5 percent and final 2.5 percent of 
         *                  the signal will be tapered with the given window.
         *                  This must be in the range [0,100].
         * @throws std::invalid_argument if pct is out of range. 
         */
        void setPercentage(const double pct);
        /*!
         * @brief Gets the percentage of the signal to taper.
         * @result The percentage of the signal to which the taper will
         *         be applied.
         */
        double getPercentage(void) const;
        /*!
         * @brief Determines if the taper parameters are valid.
         * @retval True indicates that this is a correctly initialized
         *         parameter class.
         * @retval False indiates that this is an incorrectly initialized
         *         parameter class.
         */
        bool isValid(void) const;
        /*!
         * @brief Clears variables in class and restores defaults.
         */
        void clear(void);
    private:
        class TaperParametersImpl;
        std::unique_ptr<TaperParametersImpl> pImpl;
};

/*!
 * @class Taper taper.hpp "include/rtseis/processing/singleChannel/taper.hpp"
 * @brief Tapers a waveform.
 * @ingroup rtseis_postprocessing_sc_taper
 * @copyright Ben Baker distributed under the MIT license.
 */
class Taper : public std::exception
{
    public:
        /*!
         * @brief Default constructor.
         */
        Taper(void);
        /*!
         * @brief Copy constructor.
         * @param[in] taper   Taper class from which to initialize.
         */
        explicit Taper(const Taper &taper);
        /*!
         * @brief Constructs a taper command from the parameters.
         * @param[in] parameters  The taper parameters.
         * @throw std::invalid_argument if the parameters are invalid.
         */
        explicit Taper(const TaperParameters &parameters);
        /*!
         * @brief Copy assignment operator.
         * @param[in] demean   Demean class to copy.
         * @result A deep copy of the demean class.
         */
        Taper& operator=(const Taper &Taper);
        /*!
         * @brief Default constructor.
         */
        ~Taper(void);
        /*!
         * @brief Sets the taper parameters.
         * @param[in] parameters  A correctly initialized taper parameters
         *                        class.
         * @throw std::invalid_argument if the parameters are invalid.
         */
        void setParameters(const TaperParameters &parameters);
        /*!
         * @brief Applies the taper to the data.
         * @param[in] nx   Number of points in the signal.
         * @param[in] x    The signal to taper.  This has dimension [nx].
         * @param[out] y   The tapered signal.  This has dimension [nx].
         * @throw std::invalid_argument if the parameters are invalid. 
         */
        void apply(const int nx, const double x[], double y[]);
         
        /*!
         * @brief Clears the memory and restores the defaults.
         */
        void clear(void);
    private:
         class TaperImpl;
         std::unique_ptr<TaperImpl> pImpl;
};


}; // End Single channel
}; // End Post-processing
}; // End RTSeis
#endif
