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
 * @brief Computes the copysign of the data:
 *        \f$ 
 *          \mathop{copysign}(x)
 *        = \begin{cases}
 *            +1 & x \ge +0 \\
 *            -1 & x \le -0
 *          \end{cases}
 *        \f$.
 * @ingroup rtseis_modules
 */
/*!
 * @defgroup rtseis_modules_onebit_parameters Parameters
 * @brief Defines the parameters for the one-bit normalization module.
 * @ingroup rtseis_modules_onebit
 */

/*!
 * @brief Default construtor.
 * @param[in] lrt        Flag indicating whether or not this is for real-time.
 * @param[in] precision  Defines the precision.  By default this is double.
 * @ingroup rtseis_modules_onebit_parameters 
 */
OneBitNormalizationParameters::OneBitNormalizationParameters(
    const bool lrt,
    const enum rtseisPrecision_enum precision) :
    precision_(precision),
    lrt_(lrt),
    linit_(true)
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
    lrt_ = parameters.lrt_;
    linit_ = parameters.linit_;
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
    lrt_ = false;
    linit_ = true; // One-bit normalization is always ready to roll
    return;
}
/*!
 * @brief Gets the precision of the module.
 * @result The precision of the module.
 * @ingroup rtseis_modules_onebit_parameters
 */
enum rtseisPrecision_enum
    OneBitNormalizationParameters::getPrecision(void) const
{
    return precision_;
}
/*!
 * @brief Returns whether or not the module is for real-time application.
 * @retval A flag indicating whether or not this is for real-time use.
 * @ingroup rtseis_modules_onebit_parameters
 */
bool OneBitNormalizationParameters::getIsRealTime(void) const
{
    return lrt_;
}
/*!
 * @brief Returns whether or not the parameters class is ready to be
 *        passed onto the OneBitNormalization class for use.
 * @retval True indicates this is a correctly initialized parameter class.
 */
bool OneBitNormalizationParameters::isInitialized(void) const
{
    return linit_;
}
/*!
 * @brief Default constructor.
 * @ingroup rtseis_modules_onebit
 */
OneBitNormalization::OneBitNormalization(void)
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
    linit_ = onebit.linit_;
    return *this;
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
    linit_ = false;
    parms_.clear();
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
