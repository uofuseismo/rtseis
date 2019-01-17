#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/demean.hpp"

using namespace RTSeis::Modules;

/*!
 * @defgroup rtseis_modules_demean Demean
 * @brief Removes the mean from the data.
 * @ingroup rtseis_modules
 */
/*!
 * @defgroup rtseis_modules_demean_parameters Parameters
 * @brief Defines the parameters for the demean module.
 * @ingroup rtseis_modules_demean
 */

/*!
 * @brief Default construtor.
 * @param[in] precision  Defines the precision.  By default this is double.
 * @ingroup rtseis_modules_detrend_parameters 
 */
DemeanParameters::DemeanParameters(const enum rtseisPrecision_enum precision) :
    precision_(precision),
    lrt_(false),
    linit_(true)
{
    return;
}
/*!
 * @brief Initializes parameters from parameters class.
 * @param[in] parameters   Parameters class to initialize from.
 * @ingroup rtseis_modules_detrend_parameters
 */
DemeanParameters::DemeanParameters(const DemeanParameters &parameters)
{
    *this = parameters;
    return;
}
/*!
 * @brief Copy assignement operator.
 * @param[in] parameters  Parameters from copy.
 * @result A copy of the input parameters.
 * @ingroup rtseis_modules_detrend_parameters
 */
DemeanParameters&
    DemeanParameters::operator=(const DemeanParameters &parameters)
{
    if (&parameters == this){return *this;}
    precision_ = parameters.precision_;
    lrt_ = parameters.lrt_;
    linit_ = parameters.linit_;
    return *this;
}
/*!
 * @brief Destructor.
 * @ingroup rtseis_modules_detrend_parameters
 */
DemeanParameters::~DemeanParameters(void)
{
    clear();
    return;
}
/*!
 * @brief Clears variables in class and restores defaults.
 * @ingroup rtseis_modules_detrend_parameters
 */
void DemeanParameters::clear(void)
{
    precision_ = defaultPrecision_;
    lrt_ = false;
    linit_ = true; // Demeaning is always ready to roll
    return;
}
/*!
 * @brief Gets the precision of the module.
 * @result The precision of the module.
 * @ingroup rtseis_modules_detrend_parameters
 */
enum rtseisPrecision_enum DemeanParameters::getPrecision(void) const
{
    return precision_;
}
/*!
 * @brief Returns whether or not the module is for real-time application.
 * @retval A flag indicating whether or not this is for real-time use.
 * @ingroup rtseis_modules_detrend_parameters
 */
bool DemeanParameters::getIsRealTime(void) const
{
    return lrt_;
}
/*!
 * @brief Returns whether or not the parameters class is ready to be
 *        passed onto the Detrend class for use.
 * @retval True indicates this is a correctly initialized parameter class.
 */
bool DemeanParameters::isInitialized(void) const
{
    return linit_;
}
/*!
 * @brief Default constructor.
 * @ingroup rtseis_modules_demean
 */
Demean::Demean(void)
{
    clear();
    return;
}
/*!
 * @brief Copy constructor.
 * @param[in] demean   Demean class from which to initialize.
 * @ingroup rtseis_modules_demean
 */
Demean::Demean(const Demean &demean)
{
    *this = demean;
    return;
}
/*!
 * @brief Copy assignment operator.
 * @param[in] demean   Demean class to copy.
 * @result A deep copy of the demean class.
 * @ingroup rtseis_modules_demean
 */
Demean& Demean::operator=(const Demean &demean)
{
    if (&demean == this){return *this;}
    clear();
    mean_ = demean.mean_;
    linit_ = demean.linit_;
    return *this;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_modules_demean
 */
Demean::~Demean(void)
{
    clear();
    return;
}
/*!
 * @brief Clears the memory and restores the defaults.
 */
void Demean::clear(void)
{
    mean_ = 0;
    linit_ = true; // MOdule is always ready to roll
    parms_.clear();
    return;
}
/*!
 * @brief Sets the parameters for the demean function.
 * @param[in] precision  RTSEIS_DOUBLE indicates that double precision
 *                       calculations will be used.
 * @param[in] precision  RTSEIS_FLOAT indicates that float precision
 *                       calculations will be used.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_demean
 */
int Demean::setParameters(const DemeanParameters &parameters)
{
    clear();
    int ierr = setParameters(parameters);
    if (ierr != 0){clear();}
    return 0;
}
/*!
 * @brief Removes the mean from the data.
 * @param[in] nx   Number of points in x.
 * @param[in] x    Signal from which to remove mean.  This is an array of
 *                 dimension [nx].
 * @param[out] y   The demeaned version of x.  This is an array of
 *                 dimension [nx].
 * @result 0 indicates success.
 * @ingroup rtseis_modules_demean
 */
int Demean::demean(const int nx, const double x[], double y[])
{
    mean_ = 0;
    if (nx < 1 || x == nullptr || y == nullptr)
    {
        if (nx < 1){RTSEIS_ERRMSG("%s", "No points");}
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is null");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is null");}
    }
    if (parms_.getPrecision() == RTSEIS_DOUBLE)
    {
        computeMean_(nx, x);
        removeMean_(nx, x, y);
    }
    else
    {
        float *x32 = ippsMalloc_32f(nx);
        float *y32 = ippsMalloc_32f(nx);
        ippsConvert_64f32f(x, x32, nx);
        int ierr = demean(nx, x32, y32);
        if (ierr == 0){ippsConvert_32f64f(y32, y, nx);} 
        ippsFree(x32);
        ippsFree(y32);
        return ierr;
    }
    return 0;
}
/*!
 * @brief Removes the mean from the data.
 * @param[in] nx   Number of points in x.
 * @param[in] x    Signal from which to remove mean.  This is an array of
 *                 dimension [nx].
 * @param[out] y   The demeaned version of x.  This is an array of
 *                 dimension [nx].
 * @result 0 indicates success.
 * @ingroup rtseis_modules_demean
 */
int Demean::demean(const int nx, const float x[], float y[])
{
    mean_ = 0; 
    if (nx < 1 || x == nullptr || y == nullptr)
    {
        if (nx < 1){RTSEIS_ERRMSG("%s", "No points");}
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is null");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is null");}
    }
    if (parms_.getPrecision() == RTSEIS_FLOAT)
    {
        computeMean_(nx, x);
        removeMean_(nx, x, y);
    }
    else
    {
        double *x64 = ippsMalloc_64f(nx);
        double *y64 = ippsMalloc_64f(nx);
        ippsConvert_32f64f(x, x64, nx);
        int ierr = demean(nx, x64, y64);
        if (ierr == 0){ippsConvert_64f32f(y64, y, nx);}
        ippsFree(x64);
        ippsFree(y64);
        return ierr;
    }
    return 0;
}
/*!
 * @brief Computes the mean of the data.
 * @param[in] nx  Number of points in x.
 * @param[in] x   Signal of which to compute mean.  This has dimension [nx].
 * @result 0 indicates success.
 * @ingroup rtseis_modules_demean
 */
int Demean::computeMean_(const int nx, const double x[])
{
    mean_ = 0;
    ippsMean_64f(x, nx, &mean_); // Compute mean of input 
    return 0;
}
/*!
 * @copydoc Demean::computeMean_
 * @ingroup rtseis_modules_demean
 */
int Demean::computeMean_(const int nx, const float x[])
{
    float mean32;
    mean_ = 0;
    ippsMean_32f(x, nx, &mean32, ippAlgHintAccurate); // Compute mean of x
    mean_ = static_cast<double> (mean32);
    return 0;
}
/*!
 * @brief Removes the mean from the data.
 * @param[in] nx   Number of points in x.
 * @param[in] x    Signal from which to remove mean.  This is an array of
 *                 dimension [nx].
 * @param[out] y   The demeaned version of x.  This is an array of
 *                 dimension [nx].
 * @result 0 indicates success.
 * @ingroup rtseis_modules_demean
 */ 
int Demean::removeMean_(const int nx, const double x[], double y[])
{
    ippsSubC_64f(x, mean_, y, nx); // y - mean(x)
    return 0;
}
/*!
 * @copydoc removeMean_
 * @ingroup rtseis_modules_demean
 */
int Demean::removeMean_(const int nx, const float x[], float y[])
{
    float mean32 = static_cast<float> (mean_);
    ippsSubC_32f(x, mean32, y, nx); // y - mean(x)
    return 0;
}

