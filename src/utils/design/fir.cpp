#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include <vector>
#include <cmath>
#include "rtseis/utils/design.hpp"
#include "rtseis/log.h"
#include <ipps.h>

using namespace RTSeis::Utilities::FilterDesign;

static IppWinType classifyWindow(const FIR::Window window);
/*!
 * @defgroup rtseis_utils_design_fir FIR Design
 * @brief Utility functions for FIR filter design.  This code is originally
 *        from ISTI's ISPL and has been modified to conform with C++.
 *        Function names have also been changed to conform with rtseis's
 *        naming conventions.
 * @copyright ISTI distributed under the Apache 2 license.
 * @ingroup rtseis_utils_design
 */


/*!
 * @brief Designs an FIR lowpass filter using the window method.
 * @param[in] order   Order of filter.  The number of taps is order + 1.
 *                    This must be at least 4.
 * @param[in] r       The normalized cutoff frequency where 1 is the Nyquist.
 * @param[out] taps   The filter coefficients the FIR lowpass filter.
 *                    has dimension [order+1].
 * @param[in] window  FIR window design.  The default is a Hamming window.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design_fir
 */
int FIR::FIR1Lowpass(const int order, const double r,
                     std::vector<double> &taps,
                     const Window window)
{
    taps.clear();
    // Check inputs
    if (order < 4 || r <= 0.0 || r >= 1.0)
    {
        if (order < 4){RTSEIS_ERRMSG("order=%d must be at least 4", order);}
        if (r < 0.0 || r >= 1.0)
        {
            RTSEIS_ERRMSG("r=%f must be in range (0,1)", r);
        }
        return -1;
    }
    IppWinType winType = classifyWindow(window);
    IppBool doNormal = ippTrue; // Normalize filter coefficients
    double rFreq = r/2.0;       // IPP uses 0.5 as Nyquist
    int tapsLen = order + 1;    // Length of filter is order + 1
    int bufSize;
    IppStatus status = ippsFIRGenGetBufferSize(tapsLen, &bufSize);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Error getting buffer size");
        return -1;
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);  
    status = ippsFIRGenLowpass_64f(rFreq, pTaps, tapsLen, winType,
                                   doNormal, pBuffer);
    if (status == ippStsNoErr)
    {
        taps.assign(pTaps, pTaps+static_cast<size_t> (tapsLen));
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Error generating lowpass filter coefficients");
    }
    ippsFree(pBuffer);
    ippsFree(pTaps);
    return 0;
}
/*!
 * @brief Designs an FIR highpass filter using the window method.
 * @param[in] order   Order of filter.  The number of taps is order + 1.
 *                    This must be at least 4.
 * @param[in] r       The normalized cutoff frequency where 1 is the Nyquist.
 * @param[out] taps   The filter coefficients the FIR highpass filter.
 *                    has dimension [order+1].
 * @param[in] window  FIR window design.  The default is a Hamming window.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design_fir
 */
int FIR::FIR1Highpass(const int order, const double r,
                      std::vector<double> &taps,
                      const Window window)
{
    taps.clear();
    // Check inputs
    if (order < 4 || r <= 0.0 || r >= 1.0)
    {   
        if (order < 4){RTSEIS_ERRMSG("order=%d must be at least 4", order);}
        if (r < 0.0 || r >= 1.0)
        {
            RTSEIS_ERRMSG("r=%f must be in range (0,1)", r);
        }
        return -1;
    }
    IppWinType winType = classifyWindow(window);
    IppBool doNormal = ippTrue; // Normalize filter coefficients
    double rFreq = r/2.0;       // IPP uses 0.5 as Nyquist
    int tapsLen = order + 1;    // Length of filter is order + 1
    int bufSize;
    IppStatus status = ippsFIRGenGetBufferSize(tapsLen, &bufSize);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Error getting buffer size");
        return -1;
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);
    status = ippsFIRGenHighpass_64f(rFreq, pTaps, tapsLen, winType,
                                    doNormal, pBuffer);
    if (status == ippStsNoErr)
    {
        taps.assign(pTaps, pTaps+static_cast<size_t> (tapsLen));
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Error generating highpass filter coefficients");
    }
    ippsFree(pBuffer);
    ippsFree(pTaps);
    return 0;
}
/*!
 * @brief Designs an FIR bandpass filter using the window method.
 * @param[in] order   Order of filter.  The number of taps is order + 1.
 *                    This must be at least 4.
 * @param[in] r       Normalized low cutoff frequency and high cutoff
 *                    frequencies where 1 is the Nyquist.  This is an array
 *                    with dimension [2].
 * @param[out] taps   The filter coefficients the FIR bandpass filter.
 *                    This has dimension [order+1].
 * @param[in] window  FIR window design.  The default is a Hamming window.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design_fir
 */
int FIR::FIR1Bandpass(const int order, const double r[2],
                      std::vector<double> &taps,
                      const Window window)
{
    taps.clear();
    // Check inputs
    if (order < 4 || r[0] <= 0.0 || r[0] >= 1.0 || r[1] < r[0])
    {
        if (order < 4){RTSEIS_ERRMSG("order=%d must be at least 4", order);}
        if (r[0] < 0.0 || r[0] >= 1.0)
        {
            RTSEIS_ERRMSG("r[0]=%f must be in range (0,1)", r[0]);
        }
        if (r[1] < r[0])
        {
            RTSEIS_ERRMSG("r[1]=%lf < r[0]=%lf", r[0], r[1]);
        }
        if ((int) window < 0 || (int) window > 3)
        {
            RTSEIS_ERRMSG("Invalid window=%d", (int) window);
        }
    }
    IppWinType winType = classifyWindow(window);
    IppBool doNormal = ippTrue; // Normalize filter coefficients
    double rLo = r[0]/2.0;      // IPP uses 0.5 as Nyquist
    double rHi = r[1]/2.0;      // Ipp uses 0.5 as Nyquist
    int tapsLen = order + 1;    // Length of filter is order + 1
    int bufSize;
    IppStatus status = ippsFIRGenGetBufferSize(tapsLen, &bufSize);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Error getting buffer size");
        return -1;
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);
    status = ippsFIRGenBandpass_64f(rLo, rHi, pTaps, tapsLen, winType,
                                    doNormal, pBuffer);
    if (status == ippStsNoErr)
    {
        taps.assign(pTaps, pTaps+static_cast<size_t> (tapsLen));
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Error generating bandpass filter coefficients");
    }
    ippsFree(pBuffer);
    ippsFree(pTaps);
    return 0;
}
/*!
 * @brief Designs an FIR bandstop (notch) filter using the window method.
 * @param[in] order   Order of filter.  The number of taps is order + 1.
 *                    This must be at least 4.
 * @param[in] r       Normalized low cutoff frequency and high cutoff
 *                    frequencies where 1 is the Nyquist.  This is an array
 *                    with dimension [2].
 * @param[out] taps   The filter coefficients the FIR bandpass filter.
 *                    has dimension [order+1].
 * @param[in] window  FIR window design.  The default is a Hamming window.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design_fir
 */
int FIR::FIR1Bandstop(const int order, const double r[2],
                      std::vector<double> &taps,
                      const Window window)
{
    taps.clear();
    // Check inputs
    if (order < 4 || r[0] <= 0.0 || r[0] >= 1.0 || r[1] < r[0])
    {
        if (order < 4){RTSEIS_ERRMSG("order=%d must be at least 4", order);}
        if (r[0] < 0.0 || r[0] >= 1.0)
        {
            RTSEIS_ERRMSG("r[0]=%f must be in range (0,1)", r[0]);
        }
        if (r[1] < r[0])
        {
            RTSEIS_ERRMSG("r[1]=%lf < r[0]=%lf", r[0], r[1]);
        }
        if ((int) window < 0 || (int) window > 3)
        {
            RTSEIS_ERRMSG("Invalid window=%d", (int) window);
        }
    }
    IppWinType winType = classifyWindow(window);
    IppBool doNormal = ippTrue; // Normalize filter coefficients
    double rLo = r[0]/2.0;      // IPP uses 0.5 as Nyquist
    double rHi = r[1]/2.0;      // Ipp uses 0.5 as Nyquist
    int tapsLen = order + 1;    // Length of filter is order + 1
    int bufSize;
    IppStatus status = ippsFIRGenGetBufferSize(tapsLen, &bufSize);
    if (status != ippStsNoErr)
    {   
        RTSEIS_ERRMSG("%s", "Error getting buffer size");
        return -1; 
    }   
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);
    status = ippsFIRGenBandstop_64f(rLo, rHi, pTaps, tapsLen, winType,
                                    doNormal, pBuffer);
    if (status == ippStsNoErr)
    {
        taps.assign(pTaps, pTaps+static_cast<size_t> (tapsLen));
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Error generating bandstop filter coefficients");
    }
    ippsFree(pBuffer);
    ippsFree(pTaps);
    return 0;
}
/*!
 * @brief Utility function to return the IPP window from the given
 *        design window.
 * @para[in] window  The desired window for the FIR window-based design.
 * @retval The corresponding IPP window for FIR design.
 * @ingroup rtseis_utils_design_fir
 */
IppWinType classifyWindow(const FIR::Window window)
{
    IppWinType winType = ippWinHamming; // Matlab default
    if (window == FIR::Window::HAMMING)
    {
        winType = ippWinHamming;
    }
    else if (window == FIR::Window::BARTLETT)
    {
        winType = ippWinBartlett;
    }
    else if (window == FIR::Window::HANN)
    {
        winType = ippWinHann;
    }
    else if (window == FIR::Window::BLACKMAN_OPT)
    {
        winType = ippWinBlackman;
    }
    return winType;
}
