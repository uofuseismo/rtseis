#include <stdio.h>
#include <stdlib.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#define IPPS_CORE_SRC 1
#include "rtseis/utils/filters.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utils::Filters;

/*!
 * @defgroup rtseis_utils_filters_downsample Downsample
 * @brief This is the core implementation for downsampling a signal.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
/*!
 * @brief Default constructor.
 * @ingroup rtseis_utils_filters_downsample
 */
Downsample::Downsample(void)
{
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_utils_filters_downsample
 */
Downsample::~Downsample(void)
{
    clear();
    return;
}
/*!
 * @brief Copy constructor.
 * @param[in] downsample  Class from which to initialize.
 * @ingroup rtseis_utils_filters_downsample
 */
Downsample::Downsample(const Downsample &downsample)
{
    *this = downsample;
    return;
}
/*!
 * @brief Copy operator.
 * @param[in] downsample  Downsample class to copy.
 * @result A deep copy of the input class.
 * @ingroup rtseis_utils_filters_downsample
 */
Downsample& Downsample::operator=(const Downsample &downsample)
{
    if (&downsample == this){return *this;}
    clear();
    phase0_     = downsample.phase0_;
    downFactor_ = downsample.downFactor_;
    phase_      = downsample.phase_;
    linit_      = downsample.linit_;
    return *this;
}
/*!
 * @brief Initializes the downsampler. 
 * @param[in] downFactor   The downsampling factor.  This must be positive.
 * @param[in] lisRealTime  If true then the downsampler is for real-time
 *                         applications.
 * @param[in] lisRealTime  Otherwise, it is for post-processing.
 * @param[in] precision    Precision of the downsampler.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_downsample
 */
int Downsample::initialize(const int downFactor,
                               const bool lisRealTime,
                               const enum rtseisPrecision_enum precision)
{
    clear();
    if (downFactor < 1)
    {
        RTSEIS_ERRMSG("Downsampling factor=%d must be positive", downFactor_);
        return -1;
    }
    setPrecision(precision);
    toggleRealTime(lisRealTime);
    downFactor_ = downFactor;
    phase0_ = 0;
    phase_ = 0;
    linit_ = true; 
    return 0;
}
/*!
 * @brief Releases memory on module and sets defaults.
 * @ingroup rtseis_utils_filters_downsample
 */
void Downsample::clear(void)
{
    setPrecision(RTSEIS_DOUBLE);
    toggleRealTime(false);
    phase0_ = 0;
    downFactor_ = 0;
    phase_ = 0;
    linit_ = false;
    return;
}
/*!
 * @brief Sets the initial conditions of the downsampler which is the phase.
 * @param[in] phase  Phase of downsampler.  This must be in the range
 *                  [0, getDownsampleFactor()].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_downsample
 */
int Downsample::setInitialConditions(const int phase)
{
    if (!linit_)
    {   
        RTSEIS_ERRMSG("%s", "Downsampler not initialized");
        return -1;
    }
    resetInitialConditions();
    if (phase < 0 || phase > downFactor_ - 1)
    {
        RTSEIS_ERRMSG("phase=%d must be in range[0,%d]", phase, downFactor_-1);
        return -1; 
    }
    phase0_ = phase;
    phase_  = phase;
    return 0;
}
/*!
 * @brief Resets the initial conditions to the phase set in 
 *        setInitialConditions.  If setInitialConditions was not called
 *        then this will set the phase to 0.
 * @ingroup rtseis_utils_filters_downsample 
 */
int Downsample::resetInitialConditions(void)
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Downsampler not initialized");
        return -1;
    }
    phase_ = phase0_;
    return 0;
}
/*!
 * @brief Applies the downsampler to the data.
 * @param[in] nx       The number data points in x.
 * @param[in] x        The signal to downsample.
 * @param[in] ny       The maximum number of samples in y.  One can estimate ny
 *                     by using estimateSpace(). 
 * @param[out] nyDown  The number of defined downsampled points in y. 
 * @param[out] y       The downsampled signal.  This has dimension [ny] however
 *                     only the first [nyDown] points are defined.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_downsample
 */
int Downsample::apply(const int nx, const double x[], const int ny,
                      int *nyDown, double y[])
{
    *nyDown = 0;
    if (nx <= 0){return 0;} // Nothing to do
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "ippsDS structure not intitialized");
        return -1;
    }
    int pDstLen = estimateSpace(nx); 
    if (ny < pDstLen || pDstLen < 0)
    {
        if (pDstLen < 1)
        {
            RTSEIS_ERRMSG("%s", "Space estimate error");
            return -1;
        }
        RTSEIS_ERRMSG("ny=%d must be at least length=%d", ny, pDstLen);
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    // Handle the float case
    if (isFloatPrecision())
    {
        float *x32 = ippsMalloc_32f(nx);
        float *y32 = ippsMalloc_32f(ny);
        ippsConvert_64f32f(x, x32, nx);
        int ierr = apply(nx, x32, ny, nyDown, y32);
        ippsFree(x32);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Switched precision failed");
            ippsFree(y32);
            return -1;
        }
        ippsConvert_32f64f(y32, y, *nyDown);
        ippsFree(y32);
        return 0;
    }
    // There's really nothing to do
    if (downFactor_ == 1)
    {
        ippsCopy_64f(x, y, nx);
        *nyDown = nx;
        return 0;
    }
    // Apply downsampler
    int phase = 0;
    int pEst = pDstLen;
    if (isRealTime()){phase = phase_;}
    IppStatus status = ippsSampleDown_64f(x, nx, y, &pEst,
                                          downFactor_, &phase);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to downsample signal");
        return -1;
    }
    if (pEst < 0 || pEst >= nx){pEst = pDstLen;} // Weird unhandled exception
    *nyDown = pEst; //std::min(nx, std::max(0, pEst));
    if (isRealTime()){phase_ = phase;}
    return 0;
}
/*!
 * @brief Applies the downsampler to the data.
 * @param[in] nx       The number data points in x.
 * @param[in] x        The signal to downsample.
 * @param[in] ny       The maximum number of samples in y.  One can estimate ny
 *                     by using estimateSpace(). 
 * @param[out] nyDown  The number of defined downsampled points in y. 
 * @param[out] y       The downsampled signal.  This has dimension [ny] however
 *                     only the first [nyDown] points are defined.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_downsample
 */
int Downsample::apply(const int nx, const float x[], const int ny, 
                          int *nyDown, float y[])
{
    *nyDown = 0;
    if (nx <= 0){return 0;} // Nothing to do
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "ippsDS structure not intitialized");
        return -1;
    }
    int pDstLen = estimateSpace(nx); 
    if (ny < pDstLen || pDstLen < 0)
    {
        if (pDstLen < 1)
        {
            RTSEIS_ERRMSG("%s", "Space estimate error");
            return -1; 
        }
        RTSEIS_ERRMSG("ny=%d must be at least length=%d", ny, pDstLen);
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1; 
    }
    // Handle the double case
    if (isDoublePrecision())
    {
        double *x64 = ippsMalloc_64f(nx);
        double *y64 = ippsMalloc_64f(ny);
        ippsConvert_32f64f(x, x64, nx);
        int ierr = apply(nx, x64, ny, nyDown, y64);
        ippsFree(x64);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Switched precision failed");
            ippsFree(y64);
            return -1;
        }
        ippsConvert_64f32f(y64, y, *nyDown);
        ippsFree(y64);
        return 0;
    }
    // There's really nothing to do
    if (downFactor_ == 1)
    {   
        ippsCopy_32f(x, y, nx);
        *nyDown = nx; 
        return 0;
    }
    // Apply downsampler
    int phase = 0;
    int pEst = pDstLen;
    if (isRealTime()){phase = phase_;}
    IppStatus status = ippsSampleDown_32f(x, nx, y, &pEst,
                                          downFactor_, &phase);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to downsample signal");
        return -1; 
    }
    if (pEst < 0 || pEst >= nx){pEst = pDstLen;} // Weird unhandled exception
    *nyDown = pEst; //std::min(nx, std::max(0, pEst));
    if (isRealTime()){phase_ = phase;}
    return 0;
}
/*!
 * @brief Estimates the length of the downsampled signal.
 * @param[in] n   Number of points in the signal to downsample.
 * @result The number of points required to store the output signal.
 *         If negative then there was a failure.
 * @ingroup rtseis_utils_filters_downsample
 */
int Downsample::estimateSpace(const int n) const
{
    if (!linit_ || n <= 0)
    {
        if (!linit_){RTSEIS_ERRMSG("%s", "Class not initialized");}
        if (n < 0){RTSEIS_ERRMSG("n=%d cannot be negative", n);}
        return -1;
    }
    int phase = 0;
    if (isRealTime()){phase = phase_;}
    int pDstLen = (n + downFactor_ - 1 - phase)/downFactor_;
    return pDstLen;
}
/*!
 * @brief Service routine to get the downsampling factor.
 * @result The downsampling factor.
 * @ingroup rtseis_utils_filters_downsample
 */
int Downsample::getDownsampleFactor(void) const
{
    return downFactor_;
}
