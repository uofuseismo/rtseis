#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/iir.hpp"
#include "rtseis/ippsHelper.h"

using namespace RTSeis::Modules;

/*!
 * @defgroup rtseis_modules_iir IIR
 * @brief Infinite impulse response filtering.
 * @ingroup rtseis_modules
 * @author Ben Baker
 * @copyright Ben Baker distributed under the MIT license.
 */
/*!
 * @defgroup rtseis_modules_iir_parameters Parameters
 * @brief Defines the parameters for the IIR filter.
 * @ingroup rtseis_modules_iir
 * @author Ben Baker
 * @copyright Ben Baker distributed under the MIT license.
 */

/*!
 * @brief Default constructor.
 * @ingroup rtseis_modules_iir_parameters
 */
IIRParameters::IIRParameters(void)
{
    clear();
    return;
}
/*!
 * @brief Constructor which copies the input parameter class.
 * @ingroup rtseis_modules_iir_parameters
 */
IIRParameters::IIRParameters(const IIRParameters &parameters)
{
    *this = parameters;
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_modules_iir_parameters
 */
IIRParameters::~IIRParameters(void)
{
    clear();
    return;
}
/*!
 * @brief Resets the current filter design and deletes saved filter designs.
 * @ingroup rtseis_modules_iir_parameters
 */
void IIRParameters::clear(void)
{
    resetFilterDesign();
    return;
}
/*!
 * @brief Resets the filter design parameters.  The saved filter designs
 *        will be preserved. 
 * @ingroup rtseis_modules_iir_parameters
 */
void IIRParameters::resetFilterDesign(void)
{
    ba_.clear();
    prototype_ = defaultPrototype_;
    bandtype_ = defaultBandtype_;
    w0_[0] = 0;
    w0_[1] = 0;
    fs_ = fsDefault_; 
    order_ = defaultOrder_;
    precision_ = defaultPrecision_;
    lrt_ = lrtDefault_;
    return;
}
/*!
 * @brief Creates a filter from an analog prototype.
 * @ingroup rtseis_modules_iir_parameters
 */
/*
int IIRParameters::design(
    const int order, )
{
    resetFilterDesign(); 
}
*/
/*!
 * @brief Sets a custom filter in zero, pole, gain format.
 * @param[in] zpk   The custom filter to set.
 * @result 0 indicates success.
 * @ingroup 0 indicates success.
 * @ingroup rtseis_modules_iir_parameters
 */
int IIRParameters::setCustomFilter(const ZPK zpk)
{
    resetFilterDesign();
    int ierr = RTSeis::Utils::FilterDesign::IIR::zpk2tf(zpk, ba_);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to convert filter to transfer fn");
        ba_.clear();
        return -1; 
    }
    prototype_ = Prototype::CUSTOM;
    bandtype_  = Bandtype::CUSTOM;
    return 0;
}
/*!
 * @brief Sets a custom filter in numerator/denominator format.
 * @param[in] ba   The custom filter to set.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_iir_parameters
 */
int IIRParameters::setCustomFilter(const BA ba)
{
    resetFilterDesign();
    if (ba.getNumberOfNumeratorCoefficients() < 1)
    {
        RTSEIS_ERRMSG("%s", "No numerator coefficients");
        return -1;
    }
    if (ba.getNumberOfDenominatorCoefficients() < 2)
    {
        if (ba.getNumberOfDenominatorCoefficients() < 1)
        {
            RTSEIS_ERRMSG("%s", "No denominator coefficients");
            return -1;
        }
        else
        {
            RTSEIS_WARNMSG("%s", "Consider using an FIR filter");
        }
    }
    if (ba.isZeroDenominator())
    {
        RTSEIS_ERRMSG("%s", "All denominator coefficients are zero");
        return -1;
    }
    ba_ = ba;
    prototype_ = Prototype::CUSTOM;
    bandtype_  = Bandtype::CUSTOM;
    return 0;
}
/*!
 * @brief Sets the precision.
 * @param[in] precision  Precision of calculations in module.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_iir_parameters
 */
int IIRParameters::setPrecision(
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
 * @brief Sets the corner frequency for lowpass or highpass filter
 *        design.
 * @param[in] f0    The corner frequency in Hz.  This must be positive
 *                  but less than the Nyquist frequency.
 * @param[in] fs    The sampling rate in Hz.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_iir_parameters
 */
int IIRParameters::setCornerFrequency(const double f0, const double fs)
{
    w0_[0] = 0;
    w0_[1] = 0;
    if (f0 <= 0 || fs <= 0)
    {
        if (f0 <= 0){RTSEIS_ERRMSG("f0=%lf must be positive", f0);}
        if (fs <= 0){RTSEIS_ERRMSG("fs=%lf must be postiive", fs);}
        return -1;
    }
    double fnyq = fs/2; // Nyquist
    if (f0 >= fnyq)
    {
        RTSEIS_ERRMSG("f0=%lf cannot exceed nyquist=%lf", f0, fnyq);
        return -1;
    }
    w0_[0] = f0/fnyq;
    return 0; 
}
/*!
 * @brief Sets the corner frequencies for bandpass or bandstop filter
 *        design.
 * @param[in] f0    The low-corner frequency in Hz.  This must be positive
 *                  but less f1.
 * @param[in] f1    The high-corner frequency in Hz.  This must be greater
 *                  than f0 but less than the Nyquist frequency.
 * @param[in] fs    The sampling rate in Hz.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_iir_parameters
 */
int IIRParameters::setCornerFrequencies(const double f0, const double f1,
                                        const double fs) 
{
    w0_[0] = 0;
    w0_[1] = 0;
    if (f0 <= 0 || fs <= 0 || f1 <= f0)
    {   
        if (f0 <= 0){RTSEIS_ERRMSG("f0=%lf must be positive", f0);}
        if (fs <= 0){RTSEIS_ERRMSG("fs=%lf must be postiive", fs);}
        if (f1 <= f0)
        {
            RTSEIS_ERRMSG("f1=%lf must be greater than f0=%lf", f0, f1);
        }
        return -1; 
    }   
    double fnyq = fs/2; // Nyquist
    if (f1 >= fnyq)
    {   
        RTSEIS_ERRMSG("f1=%lf cannot exceed nyquist=%lf", f1, fnyq);
        return -1; 
    }
    w0_[0] = f0/fnyq;
    w0_[1] = f1/fnyq;
    return 0;  
}
/*!
 * @brief Gets the precision of the module.
 * @result The precision of the module.
 * @ingroup rtseis_modules_iir_parameters
 */
enum rtseisPrecision_enum IIRParameters::getPrecision(void) const
{
    return precision_;
}
/*!
 * @brief Toggles the filter for real-time application.
 * @param[in] lrt   Flag indicating whether or not the filter is for
 *                  real-time processing.
 * @ingroup rtseis_modules_iir_parameters
 */
void IIRParameters::setIsRealTime(const bool lrt)
{
    lrt_ = lrt;
    return;
}
/*!
 * @brief Determines if the filter is for real-time processing or not.
 * @retval If true then the filter is for real-time processing.
 * @retval If false then the filter is for post-processing.
 * @ingroup rtseis_modules_iir_parameters 
 */
bool IIRParameters::getIsRealTime(void) const
{
    return lrt_;
}
/*!
 * @brief Gets the filter.
 * @result The IIR filter.
 * @ingroup rtseis_modules_iir_parameters
 */
BA IIRParameters::getFilter(void) const
{
   return ba_;
}
//============================================================================//
//                           The filter implementation                        //
//============================================================================//
/*!
 * @brief Default constructor.
 * @ingroup rtseis_modules_iir
 */
IIR::IIR(void)
{
    memset(&ippsIIR_, 0, sizeof(struct ippsIIRFilter_struct));
    return;
}
/*!
 * @brief Constructs the class from the IIR filtering parameters.
 * @ingroup rtseis_modules_iir
 */
IIR::IIR(const IIRParameters &parameters)
{
    memset(&ippsIIR_, 0, sizeof(struct ippsIIRFilter_struct)); 
    int ierr = setParameters(parameters);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize IIR filter");
    }
    return; 
}
/*!
 * @brief Class destructor.
 * @ingroup rtseis_modules_iir
 */
IIR::~IIR(void)
{
    rtseis_ippsIIRFilter_finalize(&ippsIIR_);
    return;
}
/*!
 * @brief Sets the IIR filter parameters.
 * @param[in] parameters  The filtering parameters from which to initialize
 *                        the IIR filter.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_iir
 */
int IIR::setParameters(const IIRParameters &parameters)
{
    rtseis_ippsIIRFilter_finalize(&ippsIIR_);
    BA ba = parameters.getFilter();
    if (ba.getNumberOfNumeratorCoefficients() == 0)
    {
        RTSEIS_ERRMSG("%s", "No numerator coefficients");
        return -1;
    }
    if (ba.getNumberOfDenominatorCoefficients() == 0)
    {
        RTSEIS_ERRMSG("%s", "No denominator coefficients");
        return -1;
    }
    if (ba.isZeroDenominator())
    {
        RTSEIS_ERRMSG("%s", "All denominator coefficients are zero");
        return -1;
    }
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients(); 
    int nb = static_cast<int> (ba.getNumberOfNumeratorCoefficients());
    int na = static_cast<int> (ba.getNumberOfDenominatorCoefficients());
    const double *bPtr = b.data();
    const double *aPtr = a.data();
    bool lrt = parms_.getIsRealTime();
    enum rtseisPrecision_enum precision = parms_.getPrecision();
    int ierr = rtseis_ippsIIRFilter_initialize(nb, bPtr, na, aPtr, lrt,
                                               precision, &ippsIIR_);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize IIR filter");
        rtseis_ippsIIRFilter_finalize(&ippsIIR_); 
        return -1;
    }
    return 0;
}
