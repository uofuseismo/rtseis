#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/fir.hpp"

using namespace RTSeis::Utilities::FilterRepresentations;
using namespace RTSeis::Modules;

/*!
 * @defgroup rtseis_modules_fir FIR
 * @brief Finite impulse response filtering.
 * @ingroup rtseis_modules
 * @author Ben Baker
 * @copyright Ben Baker distributed under the MIT license.
 */
/*!
 * @defgroup rtseis_modules_fir_parameters Parameters
 * @brief Defines the parameters for the FIR filter.
 * @ingroup rtseis_modules_fir
 * @author Ben Baker
 * @copyright Ben Baker distributed under the MIT license.
 */

/*!
 * @brief Default constructor.
 * @ingroup rtseis_modules_fir_parameters
 */
FIRParameters::FIRParameters(void)
{
    clear();
    return;
}
/*!
 * @brief Constructor which copies the input parameter class.
 * @ingroup rtseis_modules_fir_parameters
 */
FIRParameters::FIRParameters(const FIRParameters &parameters)
{
    *this = parameters;
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_modules_fir_parameters
 */
FIRParameters::~FIRParameters(void)
{
    clear();
    return;
}
/*!
 * @brief Resets the current filter design and deletes saved filter designs.
 * @ingroup rtseis_modules_fir_parameters 
 */
void FIRParameters::clear(void)
{
    resetFilterDesign();
    return;
}
/*!
 * @brief Resets the filter design parameters.  The saved filter designs
 *        will be preserved. 
 * @ingroup rtseis_modules_fir_parameters
 */
void FIRParameters::resetFilterDesign(void)
{
    ba_.clear();
    window_   = defaultWindow_;
    bandtype_ = defaultBandtype_;
    w0_[0] = 0;
    w0_[1] = 0;
    fs_ = fsDefault_; 
    order_ = defaultOrder_;
    precision_ = defaultPrecision_;
    lrt_ = lrtDefault_;
    lhaveFilter_ = false;
    return;
}
/*!
 * @brief Designs a lowpass or highpass FIR filter using an FIR
 *        window-based design.
 * @param[in] order   The order of the filter.  This must be at least 4.
 * @param[in] fs      The sampling frequency in Hz.
 * @param[in] f0      The corner frequency of the filter in Hz.  This
 *                    must be greater than 0 and less than the Nyquist.
 * @param[in] btype   The bandtype, i.e., lowpass or highpass.
 * @param[in] window  The window on which to base the design.  The default
 *                    is a Hanning window.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_fir_parameters
 */
int FIRParameters::designFromWindow(const int order, const double fs, 
                                    const double f0, 
                                    const Bandtype btype,
                                    const Window window)
{
    // Check inputs
    resetFilterDesign();
    if (f0 <= 0 || fs <= 0 || order < 4 || window == Window::CUSTOM)
    {
        if (f0 <= 0){RTSEIS_ERRMSG("f0=%lf must be positive", f0);}
        if (fs <= 0){RTSEIS_ERRMSG("fs=%lf must be positive", fs);}
        if (order < 4){RTSEIS_ERRMSG("order=%d must exceed 4", order);}
        if (window == Window::CUSTOM){RTSEIS_ERRMSG("%s", "Invalid window");}
        return -1;
    }
    double fnyq = fs/2;
    if (f0 >= fnyq)
    {
        RTSEIS_ERRMSG("f0=%lf cannot exceed Nyquist=%lf", fnyq, f0);
        return -1;
    }
    // Create normalized frequency
    double r = f0/fnyq;
    // Get the window
    Utilities::FilterDesign::FIR::Window windowPass;
    if (window == Window::HAMMING)
    {
        windowPass = Utilities::FilterDesign::FIR::Window::HAMMING;
    }
    else if (window == Window::BARTLETT)
    {
        windowPass = Utilities::FilterDesign::FIR::Window::BARTLETT;
    }
    else if (window == Window::HANN)
    {
        windowPass = Utilities::FilterDesign::FIR::Window::HANN;
    }
    else if (window == Window::BLACKMAN_OPT)
    {
        windowPass = Utilities::FilterDesign::FIR::Window::BLACKMAN_OPT;
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Cannot design a custom window");
        return -1;
    }
    // Design the appropriate filter
    std::vector<double> pTaps;
    int ierr = 0;
    if (btype == Bandtype::LOWPASS)
    {
        ierr = Utilities::FilterDesign::FIR::FIR1Lowpass(order, r, pTaps,
                                                         windowPass);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Lowpass filter design failed");
            return -1;
        }
    }
    else if (btype == Bandtype::HIGHPASS)
    {
        ierr = Utilities::FilterDesign::FIR::FIR1Highpass(order, r, pTaps,
                                                          windowPass);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Highpass filter design failed");
            return -1;
        }
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Two corner freqs needed for bandstop/pass design");
        return -1;
    }
    // Save the filter design
    ba_ = Utilities::FilterRepresentations::BA(pTaps);
    window_ = window;
    bandtype_ = btype;
    w0_[0] = f0;
    w0_[1] = 0;
    fs_ = fs;
    order_ = order;
    lhaveFilter_ = true;
    return ierr;
}
/*!
 * @brief Designs a bandpass or bandstop FIR filter using an FIR
 *        window-based design.
 * @param[in] order   The order of the filter.  This must be at least 4.
 * @param[in] fs      The sampling frequency in Hz.
 * @param[in] f0      The low-corner frequency of the filter in Hz.  This
 *                    must be greater than 0 and f1.
 * @param[in] f1      The high-corner frequency of the filter in Hz.  This
 *                    must be greater than f0 and less than the Nyquist.
 * @param[in] btype   The bandtype, i.e., bandpass or bandstop.
 * @param[in] window  The window on which to base the design.  The default
 *                    is a Hanning window.
 * @result 0 indicates success.
 * @ingroup rtseis_modules_fir_parameters
 */
int FIRParameters::designFromWindow(const int order, const double fs, 
                                    const double f0, const double f1,
                                    const Bandtype btype,
                                    const Window window)
{
    // Check inputs
    resetFilterDesign();
    if (f0 <= 0 || fs <= 0 || order < 4 || f1 <= f0 || window == Window::CUSTOM)
    {
        if (f0 <= 0){RTSEIS_ERRMSG("f0=%lf must be positive", f0);}
        if (fs <= 0){RTSEIS_ERRMSG("fs=%lf must be positive", fs);}
        if (order < 4){RTSEIS_ERRMSG("order=%d must exceed 4", order);}
        if (f1 <= f0)
        {
            RTSEIS_ERRMSG("f1=%lf must be greatre than f0=%lf", f1, f0);
        }
        if (window == Window::CUSTOM){RTSEIS_ERRMSG("%s", "Invalid window");}
        return -1;
    }
    double fnyq = fs/2;
    if (f1 >= fnyq)
    {
        RTSEIS_ERRMSG("f1=%lf cannot exceed Nyquist=%lf", fnyq, f1);
        return -1;
    }
    // Create normalized frequency
    double r0 = f0/fnyq;
    double r1 = f0/fnyq;
    // Get the window
    Utilities::FilterDesign::FIR::Window windowPass;
    if (window == Window::HAMMING)
    {
        windowPass = Utilities::FilterDesign::FIR::Window::HAMMING;
    }
    else if (window == Window::BARTLETT)
    {
        windowPass = Utilities::FilterDesign::FIR::Window::BARTLETT;
    }
    else if (window == Window::HANN)
    {
        windowPass = Utilities::FilterDesign::FIR::Window::HANN;
    }
    else if (window == Window::BLACKMAN_OPT)
    {
        windowPass = Utilities::FilterDesign::FIR::Window::BLACKMAN_OPT;
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Cannot design a custom window");
        return -1;
    }
    // Design the appropriate filter
    std::vector<double> pTaps;
    int ierr = 0;
    if (btype == Bandtype::LOWPASS)
    {
        ierr = Utilities::FilterDesign::FIR::FIR1Lowpass(order, r0, pTaps,
                                                         windowPass);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Lowpass filter design failed");
            return -1;
        }
        w0_[0] = f0;
        w0_[1] = 0;
    }
    else if (btype == Bandtype::HIGHPASS)
    {
        ierr = Utilities::FilterDesign::FIR::FIR1Highpass(order, r0, pTaps,
                                                          windowPass);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Highpass filter design failed");
            return -1;
        }
        w0_[0] = f0;
        w0_[1] = 0;
    }
    else if (btype == Bandtype::BANDPASS)
    {
        std::pair<double,double> r(r0, r1);
        ierr = Utilities::FilterDesign::FIR::FIR1Bandpass(order, r, pTaps,
                                                          windowPass);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Bandpass filter design failed");
            return -1;
        }
        w0_[0] = f0;
        w0_[1] = f1;
    }
    else if (btype == Bandtype::BANDSTOP)
    {
        std::pair<double,double> r(r0, r1);
        ierr = Utilities::FilterDesign::FIR::FIR1Bandstop(order, r, pTaps,
                                                          windowPass);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Bandstop filter design failed");
            return -1;
        }
        w0_[0] = f0;
        w0_[1] = f1;
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid window");
        return -1;
    }
    // Save the filter design
    ba_ = Utilities::FilterRepresentations::BA(pTaps);
    window_ = window;
    bandtype_ = btype;
    fs_ = fs;
    order_ = order;
    lhaveFilter_ = true;
    return ierr;
}
//============================================================================//

FIRFilter::FIRFilter(const FIRParameters parms)
{
    Utilities::Filters::FIRFilter *fir = firFilter_.get();
    fir->clear();
    return;
}
