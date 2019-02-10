#ifndef RTSEIS_MODULES_DEMEAN_HPP
#define RTSEIS_MODULES_DEMEAN_HPP 1
#include <memory>
#include "rtseis/config.h"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Modules
{

/*!
 * @defgroup rtseis_modules_demean_parameters Parameters
 * @brief Defines the parameters for the demean module.
 * @ingroup rtseis_modules_demean
 * @copyright Ben Baker distributed under the MIT license.
 */
class DemeanParameters
{
    public:
        /*!
         * @brief Default construtor.
         * @param[in] precision  Defines the precision.  By default this is
         *                       double.
         */
        DemeanParameters(const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        /*!
         * @brief Copy assignement operator.
         * @param[in] parameters  Parameters from copy.
         * @result A deep copy of the input parameters.
         */
        DemeanParameters& operator=(const DemeanParameters &parameters);
        /*!
         * @brief Copy constructor.
         * @param[in] parameters  Parmaters class from which to initialize.
         */
        DemeanParameters(const DemeanParameters &parameters);
        /*!
         * @brief Default constructor.
         */
        ~DemeanParameters(void);
        /*!
         * @brief Resets variables and releases memory in class.
         */
        void clear(void);
        /*!
         * @brief Gets the precision of the module.
         * @result The precision of the module.
         */
        RTSeis::Precision getPrecision(void) const;
        /*!
         * @brief Gets the processing mode of the module.
         * @result The processing mode.
         */
        RTSeis::ProcessingMode getProcessingMode(void) const;
        /*!
         * @brief Determines if the class is initialized.
         * @result True indicates that this is a correctly initialized
         *         parameter class.
         */
        bool isInitialized(void) const;
    private:
        RTSeis::Precision defaultPrecision_ = RTSeis::Precision::DOUBLE;
        RTSeis::Precision precision_ = defaultPrecision_;
        RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        bool linit_ = true; // This module is always ready to roll
};

/*!
 * @defgroup rtseis_modules_demean Demean
 * @brief Removes the mean from the data.
 * @ingroup rtseis_modules
 * @copyright Ben Baker distributed under the MIT license.
 */
class Demean
{
    public:
        /*!
         * @brief Default constructor.
         */
        Demean(void);
        /*!
         * @brief Copy constructor.
         * @param[in] demean   Demean class from which to initialize.
         */
        Demean(const Demean &demean);
        /*!
         * @brief Copy assignment operator.
         * @param[in] demean   Demean class to copy.
         * @result A deep copy of the demean class.
         */
        Demean& operator=(const Demean &demean);
        /*!
         * @brief Default destructor.
         */
        ~Demean(void);
        /*!
         * @brief Sets the parameters for the demean function.
         * @param[in] parameters  An initialized parameters class from which
         *                        to initialize the demean class.
         * @result 0 indicates success.
         */
        int setParameters(const DemeanParameters &parameters);
        /*!
         * @brief Removes the mean from the data.
         * @param[in] nx   Number of points in x.
         * @param[in] x    Signal from which to remove mean.  This is an array
         *                 of dimension [nx].
         * @param[out] y   The demeaned version of x.  This is an array of
         *                 dimension [nx].
         * @result 0 indicates success.
         */
        int demean(const int nx, const double x[], double y[]);
        int demean(const int nx, const float  x[], float  y[]);
        /*!
         * @brief Clears the memory and restores the defaults.
         */
        void clear(void);
    private:
        /// Forward declaration of implemenation.
        class DemeanImpl;
        /// PIMPL'd implementation
        std::unique_ptr<DemeanImpl> pDemean_;
};


};
};

#endif
