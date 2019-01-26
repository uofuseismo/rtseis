#ifndef RTSEIS_MODULES_ONEBIT_HPP
#define RTSEIS_MODULES_ONEBIT_HPP 1
#include "rtseis/config.h"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Modules
{

/*!
 * @defgroup rtseis_modules_onebit_parameters Parameters
 * @brief Defines the parameters for the one-bit normalization module.
 * @ingroup rtseis_modules_onebit
 * @copyright Ben Baker distributed under the MIT license.
 */
class OneBitNormalizationParameters
{
    public:
        /*!
         * @brief Default constructor.
         * @param[in] lrt  Flag indicating whether or not this is for real-time.
         *                 By default this is for post-processing.
         * @param[in] precision  Defines the precision.  By default this
         *                       is double precision.
         */
        OneBitNormalizationParameters(
            const bool lrt = false,
            const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        /*!
         * @brief Copy assignement operator.
         * @param[in] parameters  Parameters class to copy.
         * @result A deep copy of the input parameters.
         */
        OneBitNormalizationParameters& operator=(const OneBitNormalizationParameters &parameters);
        /*! 
         * @brief Copy constructor.
         * @param[in] parameters  Parameters class to initialize from.
         */
        OneBitNormalizationParameters(const OneBitNormalizationParameters &parameters);
        /*!
         * @brief Default destructor.
         */
        ~OneBitNormalizationParameters(void);
        /*!
         * @brief Clears variables in class and restores defaults.
         *        This class will have to be re-initialized to use again.
         */
        virtual void clear(void);
        /*!
         * @brief Determines if the class is initialized.
         * @retval If true then the class is inititalized.
         * @retval If false then the class is not initialized.
         */
        virtual bool isInitialized(void) const;
        /*!
         * @brief Determines if the class is for real-time application.
         * @retval If true then the class is for real-time application.
         * @retval If false then the class is not for real-time application.
         */
        bool isRealTime(void) const;
        /*!
         * @brief Determines the precision of the class.
         * @result The precision with which the underlying copysign
         *         operation will be performed.
         */
        RTSeis::Precision getPrecision(void) const;
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
/*!
 * @defgroup rtseis_modules_onebit One-Bit Normalization
 * @brief Computes the one-bit normalization of the data using the copysign
 *        \f[
 *          \mathop{copysign}(x)
 *        = \left \{
 *          \begin{array}{cc}
 *            +1 & x \ge +0 \\
 *            -1 & x \le -0
 *          \end{array}
 *          \right .
 *        \f].
 * @ingroup rtseis_modules
 * @copyright Ben Baker distributed under the MIT license.
 */
class OneBitNormalization : OneBitNormalizationParameters
{
    public:
        /*!
         * @brief Default constructor.  This module will not yet be usable
         *        until the parameters are set.
         */
        OneBitNormalization(void);
        /*!
         * @brief Copy constructor.
         * @param[in] onebit   One-bit normalization class from which to
         *                     initialize.
         */
        OneBitNormalization(const OneBitNormalization &onebit);
        /*!
         * @brief Initializes from the onebit normalization parameters.
         * @param[in] parameters  One-bit normalization parameters from which
         *                        to initialize this class.
         */
        OneBitNormalization(const OneBitNormalizationParameters &parameters);
        /*!
         * @brief Copy assignment operator.
         * @param[in] onebit   One-bit normalization class to copy.
         * @result A deep copy of the one-bit normalization class.
         */
        OneBitNormalization& operator=(const OneBitNormalization &onebit);
        /*!
         * @brief Default destructor.
         */
        ~OneBitNormalization(void);
        /*!
         * @brief Initializes the class from the parameters.
         * @param[in] parameters  One-bit normalization parameters from which
         *                        to initialize this class.
         * @result 0 indicates success.
         */
        int initialize(const OneBitNormalizationParameters &parameters);
        /*!
         * @brief Sets the initial conditions for the one-bit normalization.
         *        Note, the copysign operation is applied to each sample
         *        so there are no initial conditions to set.
         * @result 0 indicates success.
         */
        int setInitialConditions(void);
        /*!
         * @brief Applies the one-bit normalization.
         * @param[in] nx   Number of points in signal.
         * @param[in] x    Signal to one-bit normalized.  This is an array of
         *                 dimension [nx].
         * @param[out] y   The one-bit normalized signal.  This has
         *                 dimension [nx].
         * @result 0 indicates success.
         */
        int apply(const int nx, const double x[], double y[]);
        /*!
         * @brief Applies the one-bit normalization.
         * @param[in] nx   Number of points in signal.
         * @param[in] x    Signal to one-bit normalized.  This is an array of
         *                 dimension [nx].
         * @param[out] y   The one-bit normalized signal.  This has
         *                 dimension [nx].
         * @result 0 indicates success.
         */
        int apply(const int nx, const float  x[], float  y[]);
        /*!
         * @brief Resets the filter to the initial conditions specified
         *        in setInitialConditions().  Note, the copysign operation
         *        is applied each sample so there are no initial
         *        parameters to restore.
         * @result 0 indicates success.
         */
        int resetInitialConditions(void);
        /*!
         * @brief Clears variables in class and restores defaults.
         *        This class will have to be re-initialized to use again.
         */
        void clear(void) override;
        /*!
         * @brief Determines if the class is for real-time application.
         * @retval If true then the class is for real-time application.
         * @retval If false then the class is not for real-time application.
         */
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
