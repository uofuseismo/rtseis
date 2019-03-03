#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <stdexcept>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/design/fir.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/log.h"
#include <ipps.h>

using namespace RTSeis::Utilities;
using namespace RTSeis::Utilities::FilterDesign;

static IppWinType classifyWindow(const FIR::Window window);

void FIR::FIR1Lowpass(const int order, const double r,
                      FilterRepresentations::FIR &fir,
                      const Window window)
{
    fir.clear();
    // Check inputs
    if (order < 4 || r <= 0.0 || r >= 1.0)
    {
        if (order < 4)
        {
            RTSEIS_ERRMSG("order=%d must be at least 4", order);
            throw std::invalid_argument("order must be at least 4");
        }
        if (r < 0.0 || r >= 1.0)
        {
            RTSEIS_ERRMSG("r=%f must be in range (0,1)", r);
            throw std::invalid_argument("r must be in range (0,1)");
        }
        throw std::invalid_argument("invalid inputs");
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
        throw std::invalid_argument("Algorithmic failure in buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);  
    status = ippsFIRGenLowpass_64f(rFreq, pTaps, tapsLen, winType,
                                   doNormal, pBuffer);
    ippsFree(pBuffer);
    if (status == ippStsNoErr)
    {
        fir.setFilterTaps(static_cast<size_t> (tapsLen), pTaps);
    }
    else
    {
        ippsFree(pTaps);
        RTSEIS_ERRMSG("%s", "Error generating lowpass filter coefficients");
        throw std::invalid_argument("Internal error calling lowpass");
    }
    ippsFree(pTaps);
    return;
}

void FIR::FIR1Highpass(const int order, const double r,
                       FilterRepresentations::FIR &fir,
                       const Window window)
{
    fir.clear();
    // Check inputs
    if (order < 4 || r <= 0.0 || r >= 1.0)
    {   
        if (order < 4)
        {
            RTSEIS_ERRMSG("order=%d must be at least 4", order);
            throw std::invalid_argument("r must be at least 4");
        }
        if (r < 0.0 || r >= 1.0)
        {
            RTSEIS_ERRMSG("r=%f must be in range (0,1)", r);
            throw std::invalid_argument("r must be in range (0,1)");
        }
        throw std::invalid_argument("Invalid arguments");
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
        throw std::invalid_argument("Algorithmic failure in buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);
    status = ippsFIRGenHighpass_64f(rFreq, pTaps, tapsLen, winType,
                                    doNormal, pBuffer);
    ippsFree(pBuffer);
    if (status == ippStsNoErr)
    {
        fir.setFilterTaps(static_cast<size_t> (tapsLen), pTaps);
    }
    else
    {
        ippsFree(pTaps);
        RTSEIS_ERRMSG("%s", "Error generating highpass filter coefficients");
        throw std::invalid_argument("Internal error calling highass");
    }
    ippsFree(pTaps);
    return;
}

void FIR::FIR1Bandpass(const int order, const std::pair<double, double> &r,
                       FilterRepresentations::FIR &fir,
                       const Window window)
{
    fir.clear();
    // Check inputs
    const double r0 = r.first;
    const double r1 = r.second;
    if (order < 4 || r0 <= 0.0 || r0 >= 1.0 || r1 < r0)
    {
        if (order < 4)
        {
            RTSEIS_ERRMSG("order=%d must be at least 4", order);
            throw std::invalid_argument("Order must be at least 4");
        }
        if (r0 < 0.0 || r0 >= 1.0)
        {
            RTSEIS_ERRMSG("r.first=%f must be in range (0,1)", r0);
            throw std::invalid_argument("r0 must be in range (0,1)");
        }
        if (r1 < r0)
        {
            RTSEIS_ERRMSG("r.second=%lf < r.first=%lf", r0, r1);
            throw std::invalid_argument("r1 must be in range (r0,1)");
        }
        throw std::invalid_argument("Invalid arguments");
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
        throw std::invalid_argument("Algorithmic failure in buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);
    status = ippsFIRGenBandpass_64f(rLo, rHi, pTaps, tapsLen, winType,
                                    doNormal, pBuffer);
    ippsFree(pBuffer);
    if (status == ippStsNoErr)
    {
        fir.setFilterTaps(static_cast<size_t> (tapsLen), pTaps);
    }
    else
    {
        ippsFree(pTaps);
        RTSEIS_ERRMSG("%s", "Error generating bandpass filter coefficients");
        throw std::invalid_argument("Error calling bandpass");
    }
    ippsFree(pTaps);
    return;
}

void FIR::FIR1Bandstop(const int order, const std::pair<double, double> &r,
                       FilterRepresentations::FIR &fir,
                       const Window window)
{
    fir.clear();
    // Check inputs
    const double r0 = r.first;
    const double r1 = r.second;
    if (order < 4 || r0 <= 0.0 || r0 >= 1.0 || r1 < r0)
    {
        if (order < 4)
        {
            RTSEIS_ERRMSG("order=%d must be at least 4", order);
            throw std::invalid_argument("Order must be at least 4");
        }
        if (r0 < 0.0 || r0 >= 1.0)
        {
            RTSEIS_ERRMSG("r.first=%f must be in range (0,1)", r0);
            throw std::invalid_argument("r0 must be in range (0,1)");
        }
        if (r1 < r0)
        {
            RTSEIS_ERRMSG("r.second=%lf < r.first=%lf", r0, r1);
            throw std::invalid_argument("r1 must be in range (r0,1)");
        }
        throw std::invalid_argument("Invalid arguments");
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
        throw std::invalid_argument("Algorithmic failure in buffer size");
    }   
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    Ipp64f *pTaps = ippsMalloc_64f(tapsLen);
    status = ippsFIRGenBandstop_64f(rLo, rHi, pTaps, tapsLen, winType,
                                    doNormal, pBuffer);
    ippsFree(pBuffer);
    if (status == ippStsNoErr)
    {
        fir.setFilterTaps(static_cast<size_t> (tapsLen), pTaps);
    }
    else
    {
        ippsFree(pTaps);
        RTSEIS_ERRMSG("%s", "Error generating bandstop filter coefficients");
        throw std::invalid_argument("Error calling bandstop");
    }
    ippsFree(pTaps);
    return;
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
