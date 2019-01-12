#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/detrend.hpp"

using namespace RTSeis::Modules;

/*!
 * @defgroup rtseis_modules_detrend Detrend
 * @brief Removes the trend from the data.
 * @ingroup rtseis_modules
 */
/*!
 * @defgroup rtseis_modules_detrend_parameters Parameters
 * @brief Defines the parameters for the detrend module.
 * @ingroup rtseis_modules_detrend
 */

/*!
 * @brief Default construtor.
 * @ingroup rtseis_modules_detrend_parameters 
 */
DetrendParameters::DetrendParameters(void)
{
    setDefaults();
    return;
}
/*!
 * @brief Initializes parameters from parameters class.
 * @param[in] parameters   Parameters class to initialize from.
 * @ingroup rtseis_modules_detrend_parameters
 */
DetrendParameters::DetrendParameters(const DetrendParameters &parameters)
{
    *this = parameters;
    return;
}
/*!
 * @brief Destructor.
 * @ingroup rtseis_modules_detrend_parameters
 */
DetrendParameters::~DetrendParameters(void)
{
    return;
}
/*!
 * @brief Sets the default parameters.
 * @ingroup rtseis_modules_detrend_parameters
 */
void DetrendParameters::setDefaults(void)
{
    precision_ = defaultPrecision_;
    return;
}
/*!
 * @brief Sets the precision.
 * @param[in] precision  Precision of calculations in module.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_detrend_parameters
 */
int DetrendParameters::setPrecision(
    const enum rtseisPrecision_enum precision)
{
    precision_ = defaultPrecision_;
    if (precision != RTSEIS_DOUBLE && precision != RTSEIS_FLOAT)
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        return -1;
    }
    precision_ = precision;
    return 0;
}
/*!
 * @brief Gets the precision of the module.
 * @result The precision of the module.
 * @ingroup rtseis_modules_detrend_parameters
 */
enum rtseisPrecision_enum DetrendParameters::getPrecision(void) const
{
    return precision_;
}
/*!
 * @brief Returns whether or not the module is for real-time application.
 * @retval A flag indicating whether or not this is for real-time use.
 */
bool DetrendParameters::getIsRealTime(void) const
{
    return lrt_;
}
//============================================================================//
/*!
 * @brief Default constructor.
 * @ingroup rtseis_modules_detrend
 */
Detrend::Detrend(void)
{
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_modules_detrend
 */
Detrend::~Detrend(void)
{
    b0_ = 0;
    b1_ = 0;
    return;
}
/*!
 * @brief Sets the parameters for the detrend class.
 * @param[in] precision  RTSEIS_DOUBLE indicates that double precision
 *                       calculations will be used.
 * @param[in] precision  RTSEIS_FLOAT indicates that float precision
 *                       calculations will be used.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_detrend
 */
/*
int Detrend::setParameters(const enum rtseisPrecision_enum precision)
{
    precision_ = RTSEIS_DOUBLE;
    if (precision != RTSEIS_DOUBLE && precision != RTSEIS_FLOAT)
    {
        RTSEIS_ERRMSG("Invalid precision %d", (int) precision);
        return -1;
    }
    precision_ = precision;
    return 0;
}
*/
/*!
 * @brief Removes the trend from the data by fitting a best-fitting line.
 * @param[in] nx   Number of points in x.
 * @param[in] x    Signal from which to remove trend.  This is an array of
 *                 dimension [nx].
 * @param[out] y   The detrended version of x.  This is an array of
 *                 dimension [nx].
 * @result 0 indicates success.
 * @ingroup rtseis_modules_detrend
 */
int Detrend::detrend(const int nx, const double x[], double y[])
{
    b0_ = 0;
    b1_ = 0;
    if (nx < 1 || x == nullptr || y == nullptr)
    {
        if (nx < 1){RTSEIS_ERRMSG("%s", "No points");}
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is null");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is null");}
    }
    if (parms_.getPrecision() == RTSEIS_DOUBLE)
    {
        computeLinearRegressionCoeffs_(nx, x);
        removeTrend_(nx, x, y);
    }
    else
    {
        float *x32 = ippsMalloc_32f(nx);
        float *y32 = ippsMalloc_32f(nx);
        ippsConvert_64f32f(x, x32, nx);
        int ierr = detrend(nx, x32, y32);
        if (ierr == 0){ippsConvert_32f64f(y32, y, nx);}
        ippsFree(x32);
        ippsFree(y32);
        return ierr;
    }
    return 0;
}
/*!
 * @copydoc Detrend::compute
 * @ingroup rtseis_modules_detrend
 */
int Detrend::detrend(const int nx, const float x[], float y[])
{
    b0_ = 0;
    b1_ = 0;
    if (nx < 1 || x == nullptr || y == nullptr)
    {
        if (nx < 1){RTSEIS_ERRMSG("%s", "No points");}
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is null");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is null");}
    }
    if (parms_.getPrecision() == RTSEIS_FLOAT)
    {
        computeLinearRegressionCoeffs_(nx, x);
        removeTrend_(nx, x, y);
    }
    else
    {
        double *x64 = ippsMalloc_64f(nx);
        double *y64 = ippsMalloc_64f(nx);
        ippsConvert_32f64f(x, x64, nx);
        int ierr = detrend(nx, x64, y64);
        if (ierr == 0){ippsConvert_64f32f(y64, y, nx);}
        ippsFree(x64);
        ippsFree(y64);
        return ierr;
    }
    return 0;
}
/*!
 * @brief Computes the coefficients for a linear regression
 *          \f$ \hat{y}_i = b_0 + b_1 x_i \f$
 *        using IPP functions where
 *          \f$ b_1
 *           = \frac{\mathrm{Cov}(x,y)}{\mathrm{Var}(x,y)} \f$
 *        and 
 *          \f$ b_0 = \bar{y} - b_0 \bar{x} \f$.
 *        The code is modified from:
 *        https://software.intel.com/en-us/forums/intel-integrated-performance-primitives/topic/299457 
 *        where I have assumed the data is generated by linearly spaced samples.
 * @param[in] length  Length of array y.
 * @param[in] pSrcY   Points at evaluated at indices.
 * @ingroup rtseis_modules_detrend
 */
int Detrend::computeLinearRegressionCoeffs_(
    const int length, const double pSrcY[])
{
    double cov_xy;
    double mean_x;
    double mean_y;
    double var_x;
    // Mean of x - analytic formula for evenly spaced samples starting at indx 0
    // This is computed by simplifying Gauss's formula.
    uint64_t len64 = static_cast<uint64_t> (length);
    mean_x = 0.5*static_cast<double> (len64 - 1);
    // Note, the numerator is the sum of consecutive squared numbers.
    // In addition we simplify.
    //var_x = static_cast<double> ( ((len64 - 1)*len64)*(2*(len64 - 1) + 1) )
    //       /(static_cast<double> ( 6*len64 ) ) - mean_x*mean_x;
    var_x = static_cast<double> ( ((len64 - 1))*(2*(len64 - 1) + 1) )/6.0
          - mean_x*mean_x;
    ippsMean_64f(pSrcY, length, &mean_y);
    cov_xy = 0;
    #pragma omp simd reduction(+:cov_xy)
    for (int i=0; i<length; i++)
    {
        cov_xy = cov_xy + static_cast<double> (i)*pSrcY[i];
    }
    // This is computed by expanding (x_i - bar(x))*(y_i - bar(y)), using 
    // the definition of the mean, and simplifying
    cov_xy = (cov_xy/static_cast<double> (length)) - mean_x*mean_y;
    b1_ = cov_xy/var_x;
    b0_  = mean_y - b1_*mean_x;
    return 0;
}
/*!
 * @copydoc Detrend::computeLinearRegressionCoeffs_
 * @ingroup rtseis_modules_detrend
 */
int Detrend::computeLinearRegressionCoeffs_(
    const int length, const float pSrcY[])
{
    double cov_xy;
    double mean_x;
    double mean_y;
    double var_x;
    float mean_y32;
    // Mean of x - analytic formula for evenly spaced samples starting at indx 0
    // This is computed by simplifying Gauss's formula.
    uint64_t len64 = static_cast<uint64_t> (length);
    mean_x = 0.5*static_cast<double> (len64 - 1); 
    // Note, the numerator is the sum of consecutive squared numbers.
    // In addition we simplify.
    //var_x = static_cast<double> ( ((len64 - 1)*len64)*(2*(len64 - 1) + 1) )
    //       /static_cast<double> ( 6*len64 ) - mean_x*mean_x;
    var_x = static_cast<double> ( ((len64 - 1))*(2*(len64 - 1) + 1) )/6.0
          - mean_x*mean_x;
    ippsMean_32f(pSrcY, length, &mean_y32, ippAlgHintAccurate);
    mean_y = static_cast<double> (mean_y32);
    cov_xy = 0;
    #pragma omp simd reduction(+:cov_xy)
    for (int i=0; i<length; i++)
    {
        cov_xy = cov_xy
              + static_cast<double> (i)*static_cast<double> (pSrcY[i]);
    }
    // This is computed by expanding (x_i - bar(x))*(y_i - bar(y)), using 
    // the definition of the mean, and simplifying
    cov_xy = (cov_xy/static_cast<double> (length)) - mean_x*mean_y;
    b1_ = cov_xy/var_x;
    b0_  = mean_y - b1_*mean_x;
    return 0;
}
/*!
 * @brief Removes the trend from the data.
 * @param[in] length    Length of the time series.
 * @param[in] x         Time series from which to remove trend.  This has
 *                      dimension [length].
 * @param[out] y        The time series with the trend removed.  This has
 *                      dimension [length].
 * @result 0 indicates success.
 * @ingroup rtseis_modules_detrend
 */
int Detrend::removeTrend_(const int length, const double x[], double y[])
{
    #pragma omp simd
    for (int i=0; i<length; i++)
    {
        y[i] = x[i] - (b0_ + b1_*static_cast<double> (i));
    }
    return 0;
}
/*!
 * @copydoc Detrend::removeTrend_
 * @ingroup rtseis_modules_detrend
 */
int Detrend::removeTrend_(const int length, const float x[], float y[])
{
    float b0 = static_cast<float> (b0_);
    float b1 = static_cast<float> (b1_);
    #pragma omp simd
    for (int i=0; i<length; i++)
    {
        y[i] = x[i] - (b0 + b1*static_cast<float> (i));
    }
    return 0;
}
//============================================================================//
//                             ctypes interface                               //
//============================================================================//
/*
void *rtseis_modules_detrend_parameters_initialize(void)
{
    DetrendParameters::DetrendParameters *parms = new DetrendParameters::DetrendParameters; 
    parms->setDefaults();
    return parms;
}

void rtseis_modules_detrend_parameters_setDefault(void *parameters)
{
    DetrendParameters::DetrendParameters *parms
        = static_cast<DetrendParameters::DetrendParameters *> (parameters);
    parms->setDefaults();
    return;
}

void rtseis_modules_detrend_parameters_finalize(void *parameters)
{
    DetrendParameters::DetrendParameters *parms
        = static_cast<DetrendParameters::DetrendParameters *> (parameters);
    if (parms != nullptr){delete parms;}
    parameters = nullptr;
    return;
}

void *rtseis_modules_detrend_initialize(void *parameters)
{
    Detrend::Detrend *detrend = new Detrend::Detrend;
    DetrendParameters::DetrendParameters *parms
        = static_cast<DetrendParameters::DetrendParameters *> (parameters);
    detrend->setParameters(*parms);
    return static_cast<void *> (detrend);
}

int rtseis_modules_detrend_detrend(const int nx, const double x[],
                                   double y[], 
                                   void *detrend)
{
    Detrend::Detrend *module = static_cast<Detrend::Detrend *> (detrend);
    int ierr = module->detrend(nx, x, y);
    return ierr;
} 

void rtseis_modules_detrend_finalize(void *detrend)
{
    Detrend::Detrend *module = static_cast<Detrend::Detrend *> (detrend);
    delete module;
    detrend = nullptr;
    return;
}
*/
