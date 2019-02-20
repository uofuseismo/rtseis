#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/design/fir.hpp"
#include "rtseis/log.h"
#include <ipps.h>

using namespace RTSeis::Utilities::FilterDesign;

static IppWinType classifyWindow(const FIR::Window window);

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

int FIR::FIR1Bandpass(const int order, const std::pair<double, double> &r,
                      std::vector<double> &taps,
                      const Window window)
{
    taps.clear();
    // Check inputs
    const double r0 = r.first;
    const double r1 = r.second;
    if (order < 4 || r0 <= 0.0 || r0 >= 1.0 || r1 < r0)
    {
        if (order < 4){RTSEIS_ERRMSG("order=%d must be at least 4", order);}
        if (r0 < 0.0 || r0 >= 1.0)
        {
            RTSEIS_ERRMSG("r.first=%f must be in range (0,1)", r0);
        }
        if (r1 < r0)
        {
            RTSEIS_ERRMSG("r.second=%lf < r.first=%lf", r0, r1);
        }
        if ((int) window < 0 || (int) window > 3)
        {
            RTSEIS_ERRMSG("Invalid window=%d", (int) window);
        }
    }
    IppWinType winType = classifyWindow(window);
    IppBool doNormal = ippTrue; // Normalize filter coefficients
    double rLo = r0/2.0;        // IPP uses 0.5 as Nyquist
    double rHi = r1/2.0;        // Ipp uses 0.5 as Nyquist
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

int FIR::FIR1Bandstop(const int order, const std::pair<double, double> &r,
                      std::vector<double> &taps,
                      const Window window)
{
    taps.clear();
    // Check inputs
    const double r0 = r.first;
    const double r1 = r.second;
    if (order < 4 || r0 <= 0.0 || r0 >= 1.0 || r1 < r0)
    {
        if (order < 4){RTSEIS_ERRMSG("order=%d must be at least 4", order);}
        if (r0 < 0.0 || r0 >= 1.0)
        {
            RTSEIS_ERRMSG("r.first=%f must be in range (0,1)", r0);
        }
        if (r1 < r0)
        {
            RTSEIS_ERRMSG("r.second=%lf < r.first=%lf", r0, r1);
        }
        if ((int) window < 0 || (int) window > 3)
        {
            RTSEIS_ERRMSG("Invalid window=%d", (int) window);
        }
    }
    IppWinType winType = classifyWindow(window);
    IppBool doNormal = ippTrue; // Normalize filter coefficients
    double rLo = r0/2.0;        // IPP uses 0.5 as Nyquist
    double rHi = r1/2.0;        // Ipp uses 0.5 as Nyquist
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
