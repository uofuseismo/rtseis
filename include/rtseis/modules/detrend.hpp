#ifndef RTSEIS_MODULES_DETREND_HPP
#define RTSEIS_MODULES_DETREND_HPP 1
#include <memory>
#include "rtseis/config.h"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Modules
{
 
/*!
 * @defgroup rtseis_modules_detrend_parameters Parameters
 * @brief Defines the parameters for the detrend module.
 * @ingroup rtseis_modules_detrend
 * @copyright Ben Baker distributed under the MIT license.
 */
class DetrendParameters
{
    public:
        /*!
         * @brief Default constructor.
         * @param[in] precision  Defines the precision.  By default this
         *                       is double.
         */
        DetrendParameters(const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        /*!
         * @brief Initializes parameters from parameters class.
         * @param[in] parameters   Parameters class to initialize from.
         */
        DetrendParameters(const DetrendParameters &parameters);
        /*!
         * @brief Copy assignement operator.
         * @param[in] parameters  Parameters from copy.
         * @result A deep copy of the input parameters.
         */
        DetrendParameters& operator=(const DetrendParameters &parameters);
        /*!
         * @brief Destructor.
         */
        ~DetrendParameters(void);
        /*!
         * @brief Clears variables in class and restores defaults.
         */
        void clear(void);
        /*!
         * @brief Gets the precision of the module.
         * @result The precision of the underlying detrend calculation.
         */
        RTSeis::Precision getPrecision(void) const;
        /*!
         * @brief Gets the precision of the module.
         * @result The precision of the module.
         */
        RTSeis::ProcessingMode getProcessingMode(void) const;
        /*!
         * @brief Returns whether or not the parameters class is ready to be
         *        passed onto the Detrend class for use.
         * @retval True indicates this is a correctly initialized parameter
         *         class.
          */
        bool isInitialized(void) const;
    private:
        RTSeis::Precision defaultPrecision_ = RTSeis::Precision::DOUBLE;
        RTSeis::Precision precision_ = defaultPrecision_;
        RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        bool linit_ = true; // This module is always ready to roll
};

/*!
 * @defgroup rtseis_modules_detrend Detrend
 * @brief Removes the trend from the data.
 * @ingroup rtseis_modules
 * @copyright Ben Baker distributed under the MIT license.
 */
class Detrend
{
    public:
        /*!
         * @brief Default constructor.
         */
        Detrend(void);
        /*!
         * @brief Copy constructor.
         * @param[in] detrend  Detrend class from which to initialize.
         */
        Detrend(const Detrend &detrend);
        /*!
         * @brief Initializes from the detrend parameters.
         * @param[in] parameters  Parameters from which to initialize.
         */
        Detrend(const DetrendParameters &parameters);
        /*!
         * @brief Copy assignment operator.
         * @param[in] detrend  Detrend class to copy to this class.
         * @result A deep copy of the input detrend class.
         */
        Detrend &operator=(const Detrend &detrend);
        /*!
         * @brief Default destructor.
         */
        ~Detrend(void);
        /*!
         * @brief Sets the parameters for the detrend class.
         * @param[in] parameters  Parameters to set.
         * @result 0 indicates success.
         */
        int setParameters(const DetrendParameters &parameters);
        /*!
         * @brief Removes the trend from the data by fitting a best-fitting
         *        line.
         * @param[in] nx   Number of points in x.
         * @param[in] x    Signal from which to remove trend.  This is an
         *                 array of dimension [nx].
         * @param[out] y   The detrended version of x.  This is an array of
         *                 dimension [nx].
         * @result 0 indicates success.
         */
        int detrend(const int nx, const double x[], double y[]);
        int detrend(const int nx, const float  x[], float  y[]);
        /*!
         * @brief Releases the memory and restores the defaults.
         */
        void clear(void);
    private:
        /// Forward declaration of implementation
        class DetrendImpl;
        /// PIMPL'd implementation
        std::unique_ptr<DetrendImpl> pDetrend_;
};


};
};

#endif
