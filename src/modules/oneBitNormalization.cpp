#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/oneBitNormalization.hpp"

using namespace RTSeis::Modules;

/*!
 * @defgroup rtseis_modules_onebit One-Bit Normalization
 * @brief Computes the one-bit normalization of the data using the copysign
 *        \f$ 
 *          \mathop{copysign}(x)
 *        = \begin{cases}
 *            +1 & x \ge +0 \\
 *            -1 & x \le -0
 *          \end{cases}
 *        \f$.
 * @ingroup rtseis_modules
 * @copyright Ben Baker distributed under the MIT license.
 */
/*!
 * @defgroup rtseis_modules_onebit_parameters Parameters
 * @brief Defines the parameters for the one-bit normalization module.
 * @ingroup rtseis_modules_onebit
 * @copyright Ben Baker distributed under the MIT license.
 */

/*!
 * @brief Default construtor.
 * @param[in] lrt        Flag indicating whether or not this is for real-time.
 * @param[in] precision  Defines the precision.  By default this is double.
 * @ingroup rtseis_modules_onebit_parameters 
 */
OneBitNormalizationParameters::OneBitNormalizationParameters(
    const bool lrt,
    const RTSeis::Precision prec) :
    precision_(prec),
    isRealTime_(lrt),
    isInitialized_(true)
{
    return;
}
/*!
 * @brief Initializes parameters from parameters class.
 * @param[in] parameters   Parameters class to initialize from.
 * @ingroup rtseis_modules_onebit_parameters
 */
OneBitNormalizationParameters::OneBitNormalizationParameters(
    const OneBitNormalizationParameters &parameters)
{
    *this = parameters;
    return;
}
/*!
 * @brief Copy assignement operator.
 * @param[in] parameters  Parameters from copy.
 * @result A copy of the input parameters.
 * @ingroup rtseis_modules_onebit_parameters
 */
OneBitNormalizationParameters&
OneBitNormalizationParameters::operator=(
    const OneBitNormalizationParameters &parameters)
{
    if (&parameters == this){return *this;}
    precision_ = parameters.precision_;
    isRealTime_ = parameters.isRealTime_;
    isInitialized_ = parameters.isInitialized_;
    return *this;
}
/*!
 * @brief Destructor.
 * @ingroup rtseis_modules_onebit_parameters
 */
OneBitNormalizationParameters::~OneBitNormalizationParameters(void)
{
    clear();
    return;
}
/*!
 * @brief Clears variables in class and restores defaults.
 * @ingroup rtseis_modules_onebit_parameters
 */
void OneBitNormalizationParameters::clear(void)
{
    precision_ = defaultPrecision_;
    isRealTime_ = false;
    isInitialized_ = false;
    return;
}
//============================================================================//
//                                   End Parameters                           //
//============================================================================//
/*!
 * @brief Default constructor.
 * @ingroup rtseis_modules_onebit
 */
OneBitNormalization::OneBitNormalization(void) :
    isInitialized_(false)
{
    clear();
    return;
}
/*!
 * @brief Copy constructor.
 * @param[in] onebit   One-bit normalization class from which to initialize.
 * @ingroup rtseis_modules_onebit
 */
OneBitNormalization::OneBitNormalization(const OneBitNormalization &onebit)
{
    *this = onebit;
    return;
}
/*!
 * @brief Initializes from the onebit normalization parameters.
 * @param[in] parameters  One-bit normalization parameters from which
 *                        to initialized.
 * @ingroup rtseis_modules_onebit
 */
OneBitNormalization::OneBitNormalization(
     const OneBitNormalizationParameters &parameters)
{
    clear();
    if (!parameters.isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Parameters are not yet set");
        return;
    }
    parms_ = parameters;
    isInitialized_ = true;
    return;
}
/*!
 * @brief Copy assignment operator.
 * @param[in] onebit   One-bit normalization class to copy.
 * @result A deep copy of the one-bit normalization class.
 * @ingroup rtseis_modules_onebit
 */
OneBitNormalization&
OneBitNormalization::operator=(const OneBitNormalization &onebit)
{
    if (&onebit == this){return *this;}
    clear();
    parms_ = onebit.parms_;
    isInitialized_ = onebit.isInitialized_;
    return *this;
}
/*!
 * @brief Resets the initial conditions.  Note, this filter does not require
 *        initial conditions so this will not actually do anything.
 * @ingroup rtseis_modules_onebit
 */
int OneBitNormalization::setInitialConditions(void)
{
    return 0;
}
/*!
 * @brief Resets the filter to the initial conditions.  Note, this filter does
 *        not use initial conditions so this function will not do anything.
 * @ingroup rtseis_modules_onebit
 */
int OneBitNormalization::resetInitialConditions(void)
{
    return 0;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_modules_onebit
 */
OneBitNormalization::~OneBitNormalization(void)
{
    clear();
    return;
}
/*!
 * @brief Clears the memory and restores the defaults.
 * @ingroup rtseis_modules_onebit
 */
void OneBitNormalization::clear(void)
{
    parms_.clear();
    isInitialized_ = false;
    return;
}
/*!
 * @brief Applies the one-bit normalization.
 * @param[in] nx   Number of points in signal.
 * @param[in] x    Signal to one-bit normalized.  This is an array of
 *                 dimension [nx].
 * @param[out] y   The one-bit normalized signal.
 * @ingroup rtseis_modules_onebit
 */
int OneBitNormalization::apply(const int nx, const double x[], double y[])
{
    if (nx <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class is not initialied");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is null");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is null");}
        return -1;
    }
    #pragma omp simd
    for (int i=0; i<nx; i++)
    {
        y[i] = std::copysign(1.0, x[i]);
    }
    return 0;
}
/*!
 * @brief Applies the one-bit normalization.
 * @param[in] nx   Number of points in signal.
 * @param[in] x    Signal to one-bit normalized.  This is an array of
 *                 dimension [nx].
 * @param[out] y   The one-bit normalized signal.
 * @ingroup rtseis_modules_onebit
 */
int OneBitNormalization::apply(const int nx, const float x[], float y[])
{
    if (nx <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class is not initialied");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is null");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is null");}
        return -1; 
    }
    #pragma omp simd
    for (int i=0; i<nx; i++)
    {
        y[i] = std::copysign(1.0f, x[i]);
    }
    return 0;
}
/*!
 * @brief Returns a flag indicating that this class is initialized.
 * @ingroup rtseis_modules_onebit
 */
bool OneBitNormalization::isInitialized(void) const
{
    return isInitialized_;
}
