#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#ifndef NDEBUG
#include <cassert>
#endif
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/filterDesign/fir.hpp"
#include "rtseis/filterDesign/enums.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include "rtseis/utilities/windowFunctions.hpp"

using namespace RTSeis::Utilities;
using namespace RTSeis::FilterDesign;

static IppWinType classifyWindow(const FIRWindow window);
static void sinc(const int n, const double x[], double sinc[]);

RTSeis::FilterRepresentations::FIR
FIR::FIR1Lowpass(const int order, const double r,
                 const FIRWindow window)
{
    RTSeis::FilterRepresentations::FIR fir;
    // Check inputs
    if (order < 4 || r <= 0.0 || r >= 1.0)
    {
        if (order < 4)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                     + " must be at least 4");
        }
        throw std::invalid_argument("r = " + std::to_string(r)
                                 + " must be in range (0,1)");
    }
    IppWinType winType = classifyWindow(window);
    IppBool doNormal = ippTrue; // Normalize filter coefficients
    double rFreq = r/2.0;       // IPP uses 0.5 as Nyquist
    int tapsLen = order + 1;    // Length of filter is order + 1
    int bufSize;
    IppStatus status = ippsFIRGenGetBufferSize(tapsLen, &bufSize);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("Error getting buffer size");
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
        throw std::runtime_error(
            "Error generating lowpass filter coefficients");
    }
    ippsFree(pTaps);
    return fir;
}

RTSeis::FilterRepresentations::FIR
FIR::FIR1Highpass(const int order, const double r,
                  const FIRWindow window)
{
    RTSeis::FilterRepresentations::FIR fir;
    // Check inputs
    if (order < 4 || r <= 0.0 || r >= 1.0)
    {   
        if (order < 4)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                     + " must be at least 4");
        }
        throw std::invalid_argument("r = " + std::to_string(r)
                                 + " must be in range (0,1)");
    }
    IppWinType winType = classifyWindow(window);
    IppBool doNormal = ippTrue; // Normalize filter coefficients
    double rFreq = r/2.0;       // IPP uses 0.5 as Nyquist
    int tapsLen = order + 1;    // Length of filter is order + 1
    int bufSize;
    IppStatus status = ippsFIRGenGetBufferSize(tapsLen, &bufSize);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("Error getting buffer size");
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
        throw std::runtime_error(
            "Error generating highpass filter coefficients");
    }
    ippsFree(pTaps);
    return fir;
}

RTSeis::FilterRepresentations::FIR
FIR::FIR1Bandpass(const int order, const std::pair<double, double> &r,
                  const FIRWindow window)
{
    RTSeis::FilterRepresentations::FIR fir;
    // Check inputs
    const double r0 = r.first;
    const double r1 = r.second;
    if (order < 4 || r0 <= 0.0 || r0 >= 1.0 || r1 < r0)
    {
        if (order < 4)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                     + " must be at least 4");
        }
        if (r0 < 0.0 || r0 >= 1.0)
        {
            throw std::invalid_argument("r.first = " + std::to_string(r0)
                                      + " must be in range (0,1)");
        }
        throw std::invalid_argument("r.second = " + std::to_string(r0)
                                  + " < r.first = " + std::to_string(r1));
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
        throw std::runtime_error("Error getting buffer size");
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
        throw std::runtime_error(
            "Error generating bandpass filter coefficients");
    }
    ippsFree(pTaps);
    return fir;
}

RTSeis::FilterRepresentations::FIR
FIR::FIR1Bandstop(const int order, const std::pair<double, double> &r,
                  const FIRWindow window)
{
    RTSeis::FilterRepresentations::FIR fir;
    // Check inputs
    const double r0 = r.first;
    const double r1 = r.second;
    if (order < 4 || r0 <= 0.0 || r0 >= 1.0 || r1 < r0)
    {
        if (order < 4)
        {
            throw std::invalid_argument("order = " + std::to_string(order)
                                     + " must be at least 4");
        }
        if (r0 < 0.0 || r0 >= 1.0)
        {
            throw std::invalid_argument("r.first = " + std::to_string(r0)
                                      + " must be in range (0,1)");
        }
        throw std::invalid_argument("r.second = " + std::to_string(r0)
                                  + " < r.first = " + std::to_string(r1));
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
        throw std::runtime_error("Error getting buffer size");
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
        throw std::runtime_error(
          "Error generating bandstop filter coefficients");
    }
    ippsFree(pTaps);
    return fir;
}

/// Hilbert transformer
std::pair<RTSeis::FilterRepresentations::FIR,RTSeis::FilterRepresentations::FIR>
FIR::HilbertTransformer(const int order, const double beta)
{
    if (order < 0){throw std::invalid_argument("order cannot be negative");}
    int n = order + 1;
    // Special case
    if (n == 1)
    {
        std::vector<double> hfiltR(1, 1);
        std::vector<double> hfiltI(1, 0);
        RTSeis::FilterRepresentations::FIR realFIR(hfiltR);
        RTSeis::FilterRepresentations::FIR imagFIR(hfiltI); 
        return std::pair(realFIR, imagFIR);
    }
    // Create a kaiser window
    std::vector<double> kaiser(n);
    double *kaiserData = kaiser.data();
    WindowFunctions::kaiser(n, &kaiserData, beta);
    // Compute the sinc function where fc = 1 and fc/2 = 0.5
    // Part 1: t = fc/2*((1-n)/2:(n-1)/2)
    std::vector<double> t(n);
    auto term1  = 0.25*static_cast<double> (1 - n);
    #pragma omp simd
    for (int i=0; i<n; i++)
    {
        auto di = static_cast<double> (i);
        t[i] = term1 + 0.5*di;
    } 
    std::vector<double> sinct(n);
    sinc(n, t.data(), sinct.data());
    // Create the complex valued FIR coefficients of 12.66 of Oppenheim and 
    // Schafer.  Note, that hfilt = sinc(t)*exp(i*pi*t) is what Matlab uses.
    // However Oppenheim and Schafer use sin*sinc in 12.67.  This uses the
    // Matlab implementation.
    std::vector<double> hfiltR(n, 0);
    std::vector<double> hfiltI(n, 0);
    double gain = 0;
    if (n%2 == 1)
    {
        // Type III has many zeros
        gain = 1; 
        hfiltR[n/2] = 1;
        // Two cases
        int istart = 0;
        if ((n/2)%2 == 0){istart = 1;}
        for (int i=istart; i<n; i=i+2)
        {
            double ks = kaiser[i]*sinct[i];
            hfiltI[i] = ks*(sin(M_PI*t[i]));
        }
    }
    else
    {
        // Type IV is, in general, non-zero
        #pragma omp simd reduction(+:gain)
        for (int i=0; i<n; ++i)
        {
            double ks = kaiser[i]*sinct[i];
            hfiltR[i] = ks*cos(M_PI*t[i]);
            gain = gain + hfiltR[i];
            hfiltI[i] = ks*sin(M_PI*t[i]);
        }
#ifndef NDEBUG
        assert(gain != 0);
#endif
        // Normalize
        ippsDivC_64f_I(gain, hfiltR.data(), n); 
        ippsDivC_64f_I(gain, hfiltI.data(), n);
    }
    // Set the filter
    RTSeis::FilterRepresentations::FIR realFIR;
    RTSeis::FilterRepresentations::FIR imagFIR; 
    realFIR.setFilterTaps(hfiltR);
    imagFIR.setFilterTaps(hfiltI);
    return std::pair(realFIR, imagFIR);
}

void sinc(const int n, const double x[], double sincx[])
{
    double pix;
    #pragma omp simd
    for (int i=0; i<n; i++)
    {
        sincx[i] = 1;
        if (x[i] != 0)
        {
            pix = M_PI*x[i];
            sincx[i] = std::sin(pix)/pix;
        }
    }
    return;
}
/*!
 * @brief Utility function to return the IPP window from the given
 *        design window.
 * @param[in] window  The desired window for the FIR window-based design.
 * @retval The corresponding IPP window for FIR design.
 * @ingroup rtseis_utils_design_fir
 */
IppWinType classifyWindow(const FIRWindow window)
{
    IppWinType winType = ippWinHamming; // Matlab default
    if (window == FIRWindow::HAMMING)
    {
        winType = ippWinHamming;
    }
    else if (window == FIRWindow::BARTLETT)
    {
        winType = ippWinBartlett;
    }
    else if (window == FIRWindow::HANN)
    {
        winType = ippWinHann;
    }
    else if (window == FIRWindow::BLACKMAN_OPT)
    {
        winType = ippWinBlackman;
    }
    return winType;
}

/*!
 * @brief Converts a linear-phase FIR filter to minimum phase.
 * @param[in] fir      The FIR phase filter to convert to minimum phase.
 *                     This filter should be symmetric.
 * @param[out] minfir  The corresponding minimum phase FIR filter.  This
 *                     will (\c fir.getNumberOfFilterTaps() + 1)/2 filter
 *                     taps.
 * @param[in] method   Defines the minimum phase conversion method. 
 * @param[in] nfft     The number of points to use for the FFT.  This should
 *                     be at least a few times larger than 
 *                     \c fir.getNumberOfFilterTaps(). 
 *                     If this is negative then it will be selected.
 * @throws std::invalid_argument if any of the arguments are invalid.
 */
/*
void FIR::minimumPhase(
    const RTSeis::FilterRepresentations::FIR &fir, 
   RTSeis::FilterRepresentations::FIR &minfir,
    const MinimumPhaseMethod method= MinimumPhaseMethod::HOMOMORPHIC,
     const int nfft =-1)
{
    const double atol = 1.e-8;
    const double rtol = 1.e-5;
    std::vector h = getFilterTaps(); 
    int lenh = h.size();
    int nhalf = lenh/2;
    int lwarn = 0;
    #pragma omp simd reduction(max: lwarn)  
    for (int i=0; i<nhalf; i++)
    {
        if (std::abs(h[i] - h[lenh-1-i]) >=  atol + rtol*std::abs(h[i]))
        {
            lwarn = 1;
        }
    }
    if (lwarn == 1)
    {
        std::cerr << "Filter not symmetric - this may fail" << std::endl;
    }
    nfftUse = nfft;
    if (nfft < 1)
    {
        double dexp = std::ceil(2*static_cast<double> (lenh - 1)/0.01);
        nfftUse = static_cast<int> (std::pow(2, dexp));
    }
    if (nfft < h.size())
    {
        throw std::invalid_argument("nfft = " + std::to_string(nfft)
                                   + " must be at least "
                                   + std::to_string(h.size()));
    }
    if (method == MinimumPhaseMethod::HILBERT)
    {
        std::vector<std::complex<double>> hft;
        double xscal = (2*M_PI)/static_cast<double> (nfftUse)
                      *static_cast<double> (nhalf);
        #pragma omp simd
        for (int i=0; i<nfftUse; i++)
        {
            double arg = static_cast<double> (i)*xscal;
            hftr[i] = std::real(hft[i]*std::exp(std::complex<double>(0, arg)));
        }
        double dp = max(hftr) - 1;
        double ds =-min(hftr);
        double xden = std::sqrt(1 + dp + ds) + std::sqrt(1 - dp + ds);
        double s = 4/(xden*xden);
        ippsAddC_64f_I(ds, hftr.data(), nfft);
        ippsMulC_64f_I(s,  hftr.data(), nfft);
        ippsSqrt_64f_I(hftr.data(), nfft);
        ippsMinThresh_64f_I(1.e-10, hftr.data(), nfft); 
    }
    else
    {
    }
    return;
}
*/
