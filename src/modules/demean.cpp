#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/demean.h"

using namespace RTSeis::Modules;

/*!
 * @defgroup rtseis_modules_demean Demean
 * @brief Removes the mean from the data.
 * @ingroup rtseis_modules
 */

/*!
 * @brief Default constructor.
 * @ingroup rtseis_modules_demean
 */
Demean::Demean(void)
{
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_modules_demean
 */
Demean::~Demean(void)
{
    mean_ = 0;
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
int Demean::setParameters(const enum rtseisPrecision_enum precision)
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
    if (precision_ == RTSEIS_DOUBLE)
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
    if (precision_ == RTSEIS_FLOAT)
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

