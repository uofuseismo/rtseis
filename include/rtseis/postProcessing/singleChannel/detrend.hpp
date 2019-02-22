#ifndef RTSEIS_POSTPROCESSING_SC_DETREND
#define RTSEIS_POSTPROCESSING_SC_DETREND 1
#include <memory>
#include <exception>
#include "rtseis/config.h"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace PostProcessing
{
namespace SingleChannel
{
 
/*!
 * @defgroup rtseis_postprocessing_sc_detrend Detrend
 * @brief Utilities for removing a best-fitting trend line from the data.
 * @ingroup rtseis_postprocessing_sc
 */
/*!
 * @class DetrendParameters detrend.hpp "include/rtseis/processing/singleChannel/detrend.hpp"
 * @brief Defines the parameters for the detrend module.
 * @ingroup rtseis_postprocessing_sc_detrend
 * @copyright Ben Baker distributed under the MIT license.
 */
class DetrendParameters : public std::exception
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
        class DetrendParms;
        std::unique_ptr<DetrendParms> pDetrendParmsImpl_; 
};

/*!
 * @class Detrend detrend.hpp "include/rtseis/processing/singleChannel/detrend.hpp"
 * @brief Removes the trend from the data.
 * @ingroup rtseis_postprocessing_sc_detrend
 * @copyright Ben Baker distributed under the MIT license.
 */
class Detrend : public std::exception
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
         * @throws std::invalid_error if parameters are invalid.
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
         * @throws std::invalid_error if parameters is invalid.
         */
        void setParameters(const DetrendParameters &parameters);
        /*!
         * @brief Removes the trend from the data by fitting a best-fitting
         *        line.
         * @param[in] nx   Number of points in x.  This must be at least 2.
         * @param[in] x    Signal from which to remove trend.  This is an
         *                 array of dimension [nx].
         * @param[out] y   The detrended version of x.  This is an array of
         *                 dimension [nx].
         * @throws std::invalid_error if the arguments are invalid.
         */
        void apply(const int nx, const double x[], double y[]);
        /*! @copydoc apply */
        void apply(const int nx, const float  x[], float  y[]);
        /*! @} */
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


}; // End sc
}; // End pp
}; // End rtseis

#endif
