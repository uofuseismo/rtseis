#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#define RTSEIS_LOGGING 1
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "rtseis/utils/design.hpp"
#include "rtseis/utils/polynomial.hpp"
#include "rtseis/utils/ipps.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utils::FilterDesign;
using namespace RTSeis::Utils::Math;

static int getNextWorstRealPoleIndex(const std::vector<bool> pmask,
                                     const std::vector<bool> isreal_pole,
                                     const std::vector<std::complex<double>> p);
static int getWorstPoleIndex(const std::vector<bool> pmask,
                             const std::vector<std::complex<double>> p);
static int isreal_p_sum(const std::vector<bool> pmask,
                        const std::vector<bool> isreal_pole);
static int nearest_real_complex_idx(const std::vector<bool> zmask,
                                    const std::vector<bool> isreal_zero,
                                    const std::vector<std::complex<double>> z,
                                    const std::complex<double> p1,
                                    const bool lreal);
static int nearest_realOrComplex_idx(const std::vector<bool> zmask,
                                     const std::vector<std::complex<double>> z,
                                     const std::complex<double> p1);
static int cmplxreal(const std::vector<std::complex<double>> z,
                     std::vector<std::complex<double>> &p,
                     std::vector<bool> &isreal,
                     const double *tolIn = nullptr);
/*!
 * @defgroup rtseis_utils_design_iir IIR Design
 * @brief Utility functions for IIR filter design.  This code is originally
 *        from ISTI's ISCL and has been modified to conform with C++.
 *        Function names have also been changed to conform with rtseis's
 *        naming conventions.
 * @copyright ISTI distributed under the Apache 2 license.
 * @ingroup rtseis_utils_design
 */

/*!
 * @param[in] Convenience function to design a digital or analog filter from
 *            an analog prototype.
 * @param[in] n        Order of filter to design.
 * @param[in] W        A scalar or length-2 array defining the critical
 *                     frequencies.  For digital filters, these are normalized
 *                     frequencies in the range [0,1] where 1 is the Nyquist
 *                     frequency in pi radians/sample (thereby making Wn
 *                     half-cycles per sample.)   For analog filters, Wn is the
 *                     angular frequency in rad/s.
 * @param[in] rp       For Chebyshev I filters this specifies the maximum ripple
 *                     in the passband in decibels.
 * @param[in] rs       For Chebyshev II filters this specifies the minimum
 *                     attenutation in the stop band in decibels.
 * @param[in] btype    The type of filter, e.g., lowpass, highpass, bandpass,
 *                     or bandstop.
 * @param[in] ftype    The IIR filter protoytpe, e.g., Butterworth, Bessel, 
 *                     Chebyshev1, or Chebyshev2.
 * @param[out] sos     The corresponding IIR filter as a cascaded series of
 *                     second order sections.
 * @param[in] lanlaog  If true then this designs an analog filter.  The
 *                     default is a digital filter.
 * @param[in] pairing  Defines the pairing strategy.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::iirfilter(const int n, const double *W,
                   const double rp, const double rs,
                   const Bandtype btype,
                   const Prototype ftype,
                   SOS &sos,
                   const bool lanalog,
                   const Pairing pairing)
{
     sos.clear();
     ZPK zpk;
     int ierr = IIR::iirfilter(n, W, rp, rs, btype, ftype, zpk, lanalog);
     if (ierr != 0)
     {
         RTSEIS_ERRMSG("%s", "Failed to create iirfilter");
         return -1;
     }
     ierr = IIR::zpk2sos(zpk, sos, pairing);
     if (ierr != 0)
     {
         RTSEIS_ERRMSG("%s", "Failed to convert zpk to sos");
         sos.clear();
         return -1;
     }
     return 0;
}
/*!
 * @param[in] Convenience function to design a digital or analog filter from
 *            an analog prototype.
 * @param[in] n        Order of filter to design.
 * @param[in] W        A scalar or length-2 array defining the critical
 *                     frequencies.  For digital filters, these are normalized
 *                     frequencies in the range [0,1] where 1 is the Nyquist
 *                     frequency in pi radians/sample (thereby making Wn
 *                     half-cycles per sample.)   For analog filters, Wn is the
 *                     angular frequency in rad/s.
 * @param[in] rp       For Chebyshev I filters this specifies the maximum ripple
 *                     in the passband in decibels.
 * @param[in] rs       For Chebyshev II filters this specifies the minimum
 *                     attenutation in the stop band in decibels.
 * @param[in] btype    The type of filter, e.g., lowpass, highpass, bandpass,
 *                     or bandstop.
 * @param[in] ftype    The IIR filter protoytpe, e.g., Butterworth, Bessel, 
 *                     Chebyshev1, or Chebyshev2.
 * @param[out] ba      The corresponding filter specified in terms of 
 *                     feedforward and feedback coefficients.
 * @param[in] lanlaog  If true then this designs an analog filter.  The
 *                     default is a digital filter.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::iirfilter(const int n, const double *W, 
                   const double rp, const double rs, 
                   const Bandtype btype,
                   const Prototype ftype,
                   BA &ba,
                   const bool lanalog)
{
     ba.clear();
     ZPK zpk;
     int ierr = IIR::iirfilter(n, W, rp, rs, btype, ftype, zpk, lanalog);
     if (ierr != 0)
     {
         RTSEIS_ERRMSG("%s", "Failed to create iirfilter");
         return -1;
     }
     ierr = IIR::zpk2tf(zpk, ba);
     if (ierr != 0)
     {
         RTSEIS_ERRMSG("%s", "Failed to convert zpk to ba");
         ba.clear();
         return -1;
     }
     return 0;
}
/*!
 * @param[in] Convenience function to design a digital or analog filter from
 *            an analog prototype.
 * @param[in] n        Order of filter to design.
 * @param[in] W        A scalar or length-2 array defining the critical
 *                     frequencies.  For digital filters, these are normalized
 *                     frequencies in the range [0,1] where 1 is the Nyquist
 *                     frequency in pi radians/sample (thereby making Wn
 *                     half-cycles per sample.)   For analog filters, Wn is the
 *                     angular frequency in rad/s.
 * @param[in] rp       For Chebyshev I filters this specifies the maximum ripple
 *                     in the passband in decibels.
 * @param[in] rs       For Chebyshev II filters this specifies the minimum
 *                     attenutation in the stop band in decibels.
 * @param[in] btype    The type of filter, e.g., lowpass, highpass, bandpass,
 *                     or bandstop.
 * @param[in] ftype    The IIR filter protoytpe, e.g., Butterworth, Bessel, 
 *                     Chebyshev1, or Chebyshev2.
 * @param[out] zpk     The corresonding filter stored in zero, pole, and gain
 *                     format.
 * @param[in] lanlaog  If true then this designs an analog filter.  The
 *                     default is a digital filter.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::iirfilter(const int n, const double *W,
                   const double rp, const double rs,
                   const Bandtype btype,
                   const Prototype ftype,
                   ZPK &zpk,
                   const bool lanalog)
{
    zpk.clear();
    if (n < 1)
    {
        RTSEIS_ERRMSG("%d must be positive", n);
        return -1;
    }
    if (W == nullptr)
    {
        RTSEIS_ERRMSG("%s", "W cannot be null");
        return -1;
    }
    // Check ripples
    if (ftype == Prototype::CHEBYSHEV1 && rp <= 0)
    {
        RTSEIS_ERRMSG("rp %lf must be positive", rp);
        return -1;
    }
    if (ftype == Prototype::CHEBYSHEV2 && rs <= 0)
    {
        RTSEIS_ERRMSG("rs %lf must be positive", rs);
        return -1;
    }
    if (lanalog){RTSEIS_WARNMSG("%s", "Analog filter design is untested");}
    // Check the low-corner frequency
    double wn[2];
    wn[0] = W[0];
    wn[1] = 0;
    if (wn[0] < 0)
    {
        RTSEIS_ERRMSG("W[0]=%lf cannot be negative", wn[0]);
        return -1;
    } 
    if (!lanalog && wn[0] > 1)
    {
        RTSEIS_ERRMSG("W[0]=%lf cannot exceed 1", wn[0]);
        return -1;
    }
    // Check the high-corner frequency
    if (btype == Bandtype::BANDPASS || btype == Bandtype::BANDSTOP)
    {
        wn[1] = W[1];
        if (wn[1] < wn[0])
        {
            RTSEIS_ERRMSG("W[1]=%lf can't be less than W[0]=%lf", wn[0], wn[1]);
            return -1;
        }
        if (wn[1] > 1)
        {
            RTSEIS_ERRMSG("W[1]=%lf cannot exced 1", wn[1]);
            return -1;
        } 
    }
    // Create the analog prototype
    int ierr;
    ZPK zpkAp;
    if (ftype == Prototype::BUTTERWORTH)
    {
        ierr = AnalogPrototype::butter(n, zpkAp);
    }
    else if (ftype == Prototype::BESSEL)
    {
        ierr = AnalogPrototype::bessel(n, zpkAp);
    }
    else if (ftype == Prototype::CHEBYSHEV1)
    {
        ierr = AnalogPrototype::cheb1ap(n, rp, zpkAp); 
    }
    else
    {
        ierr = AnalogPrototype::cheb2ap(n, rs, zpkAp);
    }
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Analog prototype failed");
        return -1;
    }
    // Pre-warp the frequencies
    double warped[2];
    double fs;
    if (lanalog != 1)
    {
        fs = 2.0;
        warped[0] = 2.0*fs*tan(M_PI*wn[0]/fs);
        warped[1] = 2.0*fs*tan(M_PI*wn[1]/fs);
    }
    else
    {
        fs = 1.0;
        warped[0] = wn[0];
        warped[1] = wn[1];
        if (wn[1] > 1.0)
        {
            RTSEIS_ERRMSG("%s", "wn[1] > 0; result may be unstable");
        }
    }
    // Transform to lowpass, bandpass, highpass, or bandstop 
    ZPK zpktf;
    if (btype == Bandtype::LOWPASS)
    {
        ierr = IIR::zpklp2lp(zpkAp, warped[0], zpktf);
    }
    else if (btype == Bandtype::HIGHPASS)
    {
        ierr = IIR::zpklp2hp(zpkAp, warped[0], zpktf);
    }
    else if (btype == Bandtype::BANDPASS)
    {
        double bw = warped[1] - warped[0];
        double w0 = std::sqrt(warped[0]*warped[1]);
        ierr = IIR::zpklp2bp(zpkAp, w0, bw, zpktf);
    }
    else if (btype == Bandtype::BANDSTOP)
    {
        double bw = warped[1] - warped[0];
        double w0 = std::sqrt(warped[0]*warped[1]);
        ierr = IIR::zpklp2bs(zpkAp, w0, bw, zpktf);
    }
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to convert filter");
        return -1;
    }
    // Find the discrete equivalent
    if (!lanalog)
    {
        ierr = zpkbilinear(zpktf, fs, zpk);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed in bilinear transform");
            return -1;
        }
    }
    else
    {
        zpk = zpktf;
    }
    return 0;
}
/*!
 * @brief Convert a filter specified as zeros, poles, and a gain to 
 *        a filter consisting of cascaded second order sections.
 * @param[in] zpk      ZPK filter to convert to second order sections.
 * @param[out] sos     The corresponding filter stored as cascaded
 *                     section order sections.
 * @param[in] pairint  The pairing strategy.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::zpk2sos(const ZPK zpk, SOS &sos, const Pairing pairing)
{
    sos.clear();
    ZPK zpkWork = zpk;
    double k = zpkWork.getGain();
    int npoles = zpkWork.getNumberOfPoles();
    int nzeros = zpkWork.getNumberOfZeros();
    if (npoles <= 0 || nzeros <= 0)
    {
        RTSEIS_WARNMSG("sos will simply scale data by %lf", k); 
        std::vector<double> bs({k, 0, 0});
        std::vector<double> as({1, 0, 0});
        sos = SOS(1, bs, as);
        return 0;
    }
    // Get sizes and handles to zeros and poles
    int np = npoles + std::max(nzeros - npoles, 0);
    int nz = nzeros + std::max(npoles - nzeros, 0);
    int nSections = std::max(np, nz + 1)/2;
    std::vector<std::complex<double>> zin = zpkWork.getZeros();
    std::vector<std::complex<double>> pin = zpkWork.getPoles();
    // Base case
    if (nSections == 1)
    {
        // (z - z_1)*(z - 0) = z^2 - z*z_1 + 0
        std::vector<double> bs({1, 0, 0});
        if (nzeros == 1)
        {
            RTSEIS_WARNMSG("%s", "nzeros = 1 not yet tested");
            if (std::abs(std::imag(zin[0])) > 1.e-15)
            {
                RTSEIS_WARNMSG("%s", "Taking real part of zin");
            }
            bs[1] =-std::real(zin[0]);
        }
        else if (nzeros == 2)
        {
            bs[1] =-2.0*std::real(zin[0]);
            bs[2] = std::pow(std::abs(zin[0]), 2);
        }
        // Introduce gain
        bs[0] = bs[0]*k;
        bs[1] = bs[1]*k;
        bs[2] = bs[2]*k;
        // (p - p_1)*(p - 0) = p^2 - p*p_1 + 0
        std::vector<double> as({1, 0, 0});
        if (npoles == 1)
        {
            RTSEIS_WARNMSG("%s", "npoles = 1 not yet tested");
            if (std::abs(std::imag(pin[0])) > 1.e-15)
            {
                RTSEIS_WARNMSG("%s", "Taking real part of pin");
            }
            as[1] =-std::real(pin[0]);
        }
        // Conjugate pairs: (p - p_1)*(p - conj(p_1)) 
        else if (npoles == 2)
        {
            as[1] =-2.0*std::real(pin[0]);
            as[2] = std::pow(std::abs(pin[0]), 2);
        }
        sos = SOS(1, bs, as);
        RTSEIS_WARNMSG("%s", "Summary of design");
        sos.print(stdout);
        return 0;
    }
    if (np%2 == 1 && pairing == Pairing::NEAREST)
    {
        np = np + 1;
        nz = nz + 1;
    }
    if (np != nz || np < 1 || nz < 1)
    {
        if (np != nz){RTSEIS_ERRMSG("Error size inconsistent %d,%d!", np, nz);}
        if (np < 1){RTSEIS_ERRMSG("Error no zeros in use %d", nz);}
        if (nz < 1){RTSEIS_ERRMSG("Error no poles in use %d", np);}
        return -1;
    }
    std::vector<std::complex<double>> z;
    std::vector<bool> isreal_zero;
    if (cmplxreal(zin, z, isreal_zero, nullptr) != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to sort zeros");
        return -1;
    }
    std::vector<std::complex<double>> p;
    std::vector<bool> isreal_pole;
    if (cmplxreal(pin, p, isreal_pole, nullptr) != 0) 
    {
        RTSEIS_ERRMSG("%s", "Failed to sort poles");
        return -1;
    }
    std::vector<bool> zmask(z.size(), false);
    std::vector<bool> pmask(p.size(), false);
    std::vector<BA> basAll;
    for (int is=0; is<nSections; is++)
    {
        std::complex<double> z1 = std::complex<double> (0, 0);
        // Select the `worst' pole
        int p1_idx = getWorstPoleIndex(pmask, p);
        if (p1_idx ==-1)
        {
            RTSEIS_ERRMSG("%s", "Pole index error");
            return -1;
        }
        std::complex<double> p1 = p[p1_idx];
        pmask[p1_idx] = true;
        int psum = isreal_p_sum(pmask, isreal_pole);
        // Initialize complex conjugates for real case
        std::complex<double> p2 = std::complex<double> (0, 0);
        std::complex<double> z2 = std::complex<double> (0, 0);
        // Pair that pole with a zero
        if (isreal_pole[p1_idx] && psum == 0)
        {
            // Special case to set a first order section
            int z1_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                  z, p1, true);
            if (z1_idx ==-1)
            {
                RTSEIS_ERRMSG("Zero index error in section %d", is+1);
                return -1;
            }
            z1 = z[z1_idx];
            zmask[z1_idx] = true;
        }
        else
        {
            int zsum = isreal_p_sum(zmask, isreal_zero);
            int z1_idx =-1;
            if (!isreal_pole[p1_idx] && zsum == 1)
            {
                // Special case to ensure choose a complex zero to pair
                // with so later (setting up a first-order section)
                z1_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                  z, p1, false);
                if (z1_idx ==-1)
                {
                    RTSEIS_ERRMSG("%s", "Zero index error 2");
                    return -1;
                }
                if (isreal_zero[z1_idx])
                {
                    RTSEIS_ERRMSG("%s", "Error this should be complex");
                    return -1;
                }
            }
            else
            {
                // Pair that pole with the closest zero (real or complex)
                z1_idx = nearest_realOrComplex_idx(zmask, z, p1);
                if (z1_idx ==-1)
                {
                    RTSEIS_ERRMSG("%s", "Zero index error 3");
                    return -1;
                } 
            }
            z1 = z[z1_idx];
            zmask[z1_idx] = true;
            // Now that we have p1 and z1, figure out what p2 and z2 need to be
            if (!isreal_pole[p1_idx])
            {
                // Complex pole, complex zero -> complex conjugates
                if (!isreal_zero[z1_idx])
                {
                    p2 = std::conj(p1);
                    z2 = std::conj(z1);
                }
                // Complex pole and real zero
                else
                {
                    p2 = std::conj(p1);
                    int z2_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                          z, p1, true);
                    if (z2_idx ==-1)
                    {
                        RTSEIS_ERRMSG("%s", "Zero2 index error 1");
                        return -1;
                    }
                    if (!isreal_zero[z2_idx])
                    {
                        RTSEIS_ERRMSG("%s", "Should be real zero!");
                        return -1;
                    }
                    z2 = z[z2_idx];
                    zmask[z2_idx] = true;
                }
            }
            // Real pole
            else
            {
                // Real pole, complex zero
                if (!isreal_zero[z1_idx])
                {
                    z2 = conj(z1);
                    int p2_idx = nearest_real_complex_idx(pmask, isreal_pole,
                                                          p, z1, true);
                    if (p2_idx ==-1)
                    {
                        RTSEIS_ERRMSG("%s", "Pole2 index error 1");
                        return -1;
                    }
                    if (!isreal_pole[p2_idx])
                    {
                        RTSEIS_ERRMSG("%s", "Should be real pole");
                        return -1;
                    }
                    p2 = p[p2_idx];
                    pmask[p2_idx] = true;
                }
                // Real pole, real zero
                else
                {
                    // Pick the next `worst' pole to use
                    int p2_idx = getNextWorstRealPoleIndex(pmask,
                                                           isreal_pole, p);
                    if (p2_idx ==-1)
                    {
                        RTSEIS_ERRMSG("%s", "Error couldn't find next real pole");
                        return -1;
                    }
                    if (!isreal_pole[p2_idx])
                    {
                        RTSEIS_ERRMSG("%s", "Error next pole must be real");
                        return -1;
                    }
                    p2 = p[p2_idx];
                    // Find a real zero to match the added pole
                    int z2_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                          z, p2, true);
                    if (z2_idx ==-1)
                    {
                        RTSEIS_ERRMSG("%s", "Zero3 index error 1");
                        return -1;
                    }
                    if (!isreal_zero[z2_idx])
                    {
                        RTSEIS_ERRMSG("%s", "Should be real zero!");
                        return -1;
                    }
                    z2 = z[z2_idx];
                    zmask[z2_idx] = true; 
                    pmask[p2_idx] = true;
                }
            }
        }
        std::vector<std::complex<double>> ptemp({p1, p2});
        std::vector<std::complex<double>> ztemp({z1, z2});
        double ktemp = 1;
        if (is == nSections - 1){ktemp = k;}
        ZPK zpkTemp = ZPK(ztemp, ptemp, ktemp);
        BA baTemp;
        int ierr = zpk2tf(zpkTemp, baTemp);  
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to create transfer function");
            return -1;
        }
        // Save the transfer function of the secon dorder section
        basAll.push_back(baTemp);
    } // Loop on sections
    // Reality check
    for (size_t ip=0; ip<pmask.size(); ip++)
    {
        if (!pmask[ip])
        {
            RTSEIS_ERRMSG("Failed to find pole %ld; nsections=%d",
                          ip, nSections);
            return -1;
        }
    }
    for (size_t iz=0; iz<zmask.size(); iz++)
    {
        if (!zmask[iz])
        {
            RTSEIS_ERRMSG("Failed to find zero %ld", iz);
            return -1;
        }
    }
    // Construct the cascading transfer function with the `worst' pole last
    int js = nSections;
    std::vector<double> bs;
    std::vector<double> as;
    for (int is=0; is<nSections; is++)
    {
        js = js - 1;
        std::vector<double> bsTemp = basAll[js].getNumeratorCoefficients(); 
        std::vector<double> asTemp = basAll[js].getDenominatorCoefficients();
        for (size_t k=0; k<3; k++)
        {
            bs.push_back(bsTemp[k]);
            as.push_back(asTemp[k]);
        }
    }
    // Finally pack the second order sections
    sos = SOS(nSections, bs, as);
    if (sos.getNumberOfSections() == 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to pack sos");
        return -1;
    }
    return 0;
}
/*!
 * @brief Computes the polynomial transfer function from a pole-zero
 *        representation.
 * @param[in] zpk  Input transfer function specified in terms of zeros,
 *                 poles, and gain.
 * @param[out] ba  Corresponding transfer function specified in terms of
 *                 numerator and denominator coefficients.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::zpk2tf(const ZPK zpk, BA &ba)
{
    ba.clear();
    double k = zpk.getGain();
    if (k == 0){RTSEIS_WARNMSG("%s", "System gain is zero");}
    // Compute polynomial representation of zeros by expanding:
    //  (z - z_1)*(z - z_2)*...*(z - z_nzeros)
    std::vector<std::complex<double>> z = zpk.getZeros();
    std::vector<std::complex<double>> bz;
    int ierr = Polynomial::poly(z, bz); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to expand numerator polynomial");
        return -1;
    }
    // Compute polynomial representation of zeros by expanding:
    //  (z - z_1)*(z - z_2)*...*(z - z_nzeros)
    std::vector<std::complex<double>> p = zpk.getPoles();
    std::vector<std::complex<double>> az;
    ierr = Polynomial::poly(p, az); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to expand denominator polynomial");
        return -1;
    }
    // Introduce gain into the numerator zeros
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<bz.size(); i++){bz[i] = k*bz[i];}
    // Take the real
    std::vector<double> b;
    b.resize(bz.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<bz.size(); i++){b[i] = std::real(bz[i]);}
    std::vector<double> a;
    a.resize(az.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<az.size(); i++){a[i] = std::real(az[i]);}
    // Pack into a transfer function struct
    ba = BA(b, a);
    return 0;
}
/*!
 * @brief Converts an analog filter to a digital filter using hte bilinear
 *        transform.  This works by transforming a set of poles and zeros from
 *        the s-plane to the digital z-plane using Tustin's method, which
 *        substitutes \f$ s = \frac{z-1}{z+1} \f$ threby maintaining the
 *        the shape of the freuqency response.
 *
 * @param[in] zpk     The analog zeros, poles, and gain to transform to the
 *                    z-plane.
 * @param[in] fs      The sampling rate in Hz.  Note, that pre-warping is
 *                    performed in this function.
 * @param[out] zpkbl  The bilinear-transformed variant of zpk.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::zpkbilinear(const ZPK zpk, const double fs, ZPK &zpkbl)
{
    zpkbl.clear();
    size_t nzeros = zpk.getNumberOfZeros();
    size_t npoles = zpk.getNumberOfPoles();
    if (nzeros > npoles)
    {
        RTSEIS_ERRMSG("%s", "Cannot have more zeros than poles");
        return -1;
    }
    std::vector<std::complex<double>> z = zpk.getZeros();
    std::vector<std::complex<double>> p = zpk.getPoles();
    // Bilinear transform
    std::complex<double> fs2 = std::complex<double> (2*fs, 0);
    std::vector<std::complex<double>> z_bl;
    z_bl.resize(npoles);
    for (size_t i=0; i<nzeros; i++)
    {
        std::complex<double> znum = fs2 + z[i];
        std::complex<double> zden = fs2 - z[i];
        z_bl[i] = znum/zden; //(fs2 + z)/(fs2 - z)
    }
    // Any zeros at infinity get moved to Nyquist frequency
    for (size_t i=nzeros; i<npoles; i++)
    {
        z_bl[i] = std::complex<double> (-1, 0);
    }
    std::vector<std::complex<double>> p_bl;
    p_bl.resize(npoles);
    for (size_t i=0; i<npoles; i++)
    {
        std::complex<double> znum = fs2 + p[i];
        std::complex<double> zden = fs2 - p[i];
        p_bl[i] = znum/zden; //(fs2 + p)/(fs2 - p) 
    }
    // Compensate for gain change
    std::complex<double> znum = std::complex<double> (1, 0);
    for (size_t i=0; i<nzeros; i++)
    {
        std::complex<double> zdif = fs2 - z[i];
        znum = znum*zdif;
    }
    std::complex<double> zden = std::complex<double> (1, 0);
    for (size_t i=0; i<npoles; i++)
    {
        std::complex<double> zdif = fs2 - p[i];
        zden = zden*zdif;
    }
    std::complex<double> zk = znum/zden;
    double k = zpk.getGain();
    double k_bl = k*std::real(zk);
    zpkbl = ZPK(z_bl, p_bl, k_bl);
    return 0;
}
/*!
 * @brief Transforms a lowpass filter prototype to a bandpass filter.
 *        The passband width is defined by [w0,w1] = [w0,w0+bw] in rad/s.
 * @param[in] zpkIn    Input lowpass filter prototype to convert.
 * @param[in] w0       Desired cutoff frequency in rad/s.
 * @param[in] bw       The desired passband width in rad/s.
 * @param[out] zpkOut  The corresponding bandpass filter.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::zpklp2bp(const ZPK zpkIn, const double w0,
                  const double bw, ZPK &zpkOut)
{
    zpkOut.clear();
    size_t nzeros = zpkIn.getNumberOfZeros();
    size_t npoles = zpkIn.getNumberOfPoles();
    if (w0 < 0 || bw < 0)
    {
        if (w0 < 0){RTSEIS_ERRMSG("w0=%lf must be non-negative", w0);}
        if (bw < 0){RTSEIS_ERRMSG("bw=%lf must be non-negative", bw);}
        return -1;
    }
    size_t nzeros_bp = nzeros + npoles;
    size_t npoles_bp = 2*npoles;
    // Define some constants
    std::complex<double> zw02 = std::complex<double> (w0*w0 , 0);
    std::complex<double> zbw2 = std::complex<double> (bw/2.0, 0);
    // Duplicate and zeros and shift from baseband to +w0 and -w0
    std::vector<std::complex<double>> z_bp;
    z_bp.resize(nzeros_bp);
    std::vector<std::complex<double>> z = zpkIn.getZeros();
    size_t j = nzeros;
    for (size_t i=0; i<nzeros; i++)
    {
        std::complex<double> zlp = zbw2*z[i];        //Scale zeros to desired bw
        std::complex<double> zlp2 = zlp*zlp;         //zlp^2
        std::complex<double> zdif = zlp2 - zw02;     //zlp^2 - w0^2
        std::complex<double> zsqrt = std::sqrt(zdif);//sqrt(zlp^2 - w0^2)
        z_bp[i] = zlp + zsqrt; //zlp+sqrt(zlp^2-w0^2)
        z_bp[j] = zlp - zsqrt; //zlp-sqrt(zlp^2-w0^2)
        j = j + 1;
    }
    // Move degree zeros to origin leaving degree zeros at infinity for BPF 
    std::complex<double> zero(0, 0);
    for (size_t i=2*nzeros; i<nzeros_bp; i++)
    {
        z_bp[i] = zero;
    }
    // Duplicate poles and shift from baseband to +w0 and -w0
    std::vector<std::complex<double>> p_bp;
    p_bp.resize(npoles_bp);
    std::vector<std::complex<double>> p = zpkIn.getPoles();
    j = npoles;
    for (size_t i=0; i<npoles; i++)
    {   
        std::complex<double> plp = zbw2*p[i];        //Scale poles to desired bw
        std::complex<double> plp2 = plp*plp;         //plp^2
        std::complex<double> pdif = plp2 - zw02;     //plp^2 - w0^2
        std::complex<double> psqrt = std::sqrt(pdif);//sqrt(plp^2 - w0^2)
        p_bp[i] = plp + psqrt; //plp+sqrt(plp^2-w0^2)
        p_bp[j] = plp - psqrt; //plp-sqrt(plp^2-w0^2)
        j = j + 1;
    }
    // Cancel out gain change from frequency scaling
    double k = zpkIn.getGain();
    double k_bp = k*std::pow(bw,npoles-nzeros); //k*bw*bw^degree
    zpkOut = ZPK(z_bp, p_bp, k_bp);
    return 0;
}
/*!
 * @brief Transforms a lowpass filter prototype to a bandstop filter.
 *        The stopband width is defined by [w0,w1] = [w0,w0+bw] in rad/s.
 * @param[in] zpkIn    Input lowpass filter prototype to convert.
 * @param[in] w0       Desired cutoff frequency in rad/s.
 * @param[in] bw       The desired stopband width in rad/s.
 * @param[out] zpkOut  The corresponding stopband filter.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::zpklp2bs(const ZPK zpkIn, const double w0,
                  const double bw, ZPK &zpkOut)
{
    zpkOut.clear();
    size_t nzeros = zpkIn.getNumberOfZeros();
    size_t npoles = zpkIn.getNumberOfPoles();
    if (w0 < 0 || bw < 0)
    {
        if (w0 < 0){RTSEIS_ERRMSG("w0=%lf must be non-negative", w0);}
        if (bw < 0){RTSEIS_ERRMSG("bw=%lf must be non-negative", bw);}
        return -1;
    }
    size_t nzeros_bs = 2*npoles;
    size_t npoles_bs = 2*npoles;
    // Define some constants
    std::complex<double> zw02 = std::complex<double> (w0*w0 , 0);
    std::complex<double> zbw2 = std::complex<double> (bw/2.0, 0);
    // Duplicate the zeros and shift from baseband to +w0 and -w0
    std::vector<std::complex<double>> z = zpkIn.getZeros();
    std::vector<std::complex<double>> z_bs;
    z_bs.resize(nzeros_bs);
    std::complex<double> zprod = std::complex<double> (1, 0); // Init product
    size_t j = nzeros;
    for (size_t i=0; i<nzeros; i++)
    {
        // Invert to a highpass filter with desired bandwidth
        std::complex<double> zhp = zbw2/z[i];         // bw2/2/z
        std::complex<double> zhp2 = zhp*zhp;          // zhp^2
        std::complex<double> zdif = zhp2 - zw02;      // zhp^2 - w0^2
        std::complex<double> zsqrt = std::sqrt(zdif); // sqrt(zhp^2 - w0^2)
        z_bs[i] = zhp + zsqrt; // zhp+sqrt(zhp^2-w0^2)
        z_bs[j] = zhp - zsqrt; // zhp-sqrt(zhp^2-w0^2)
        // Update product 
        zprod = (-z[i])*zprod;
        j = j + 1;
    }
    // Move any zeros at infinity to center of stopband
    j = 2*nzeros;
    for (size_t i=0; i < npoles - nzeros; i++)
    {
        z_bs[j] = std::complex<double> (0, w0);
        j = j + 1;
    }
    for (size_t i=0; i < npoles - nzeros; i++)
    {
        z_bs[j] = std::complex<double> (0,-w0);
        j = j + 1;
    }
    // Duplicate the poles and shift from baseband to +w0 and -w0
    std::vector<std::complex<double>> p = zpkIn.getPoles();
    std::complex<double> pprod = std::complex<double> (1, 0);
    std::vector<std::complex<double>> p_bs;
    p_bs.resize(npoles_bs);
    j = npoles;
    for (size_t i=0; i<npoles; i++)
    {
        // Invert to a highpass filter with desired bandwidth
        std::complex<double> php = zbw2/p[i];         // bw2/2/p
        std::complex<double> php2 = php*php;          // php^2
        std::complex<double> pdif = php2 - zw02;      // php^2 - w0^2
        std::complex<double> psqrt = std::sqrt(pdif); // sqrt(php^2 - w0^2)
        p_bs[i] = php + psqrt; // php+sqrt(php^2-w0^2)
        p_bs[j] = php - psqrt; // php-sqrt(php^2-w0^2)
        // Update product
        pprod = (-p[i])*pprod;
        j = j + 1;
    }
    // Cancel out gain change caused by inversion
    double k = zpkIn.getGain();
    std::complex<double> zdiv =  zprod/pprod;
    double k_bs = k*std::real(zdiv);
    zpkOut = ZPK(z_bs, p_bs, k_bs);
    return 0;
}
/*!
 * @brief Converts a lowpass filter prototype to a different cutoff frequency.
 * @param[in] zpkIn    Input lowpass filter prototype to convert.
 * @param[in] w0       Desired cutoff frequency (rad/s).
 * @param[out] zpkOut  The corresponding lowpass filter.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int IIR::zpklp2lp(const ZPK zpkIn, const double w0, ZPK &zpkOut)
{
    zpkOut.clear();
    size_t nzeros = zpkIn.getNumberOfZeros();
    size_t npoles = zpkIn.getNumberOfPoles();
    if (w0 < 0)
    {
        RTSEIS_ERRMSG("w0=%lf must be non-negative", w0);
        return -1;
    }
    if (npoles < nzeros)
    {
        RTSEIS_ERRMSG("%s", "npoles < nzeros not implemented");
        return -1;
    }
    std::vector<std::complex<double>> z = zpkIn.getZeros();
    std::vector<std::complex<double>> p = zpkIn.getPoles();
    // Scale all points radially from origin to shift cutoff frequency
    std::vector<std::complex<double>> z_lp(nzeros);
    for (size_t i=0; i<z.size(); i++)
    {
        z_lp[i] = w0*z[i];
    }
    std::vector<std::complex<double>> p_lp(npoles);
    for (size_t i=0; i<p.size(); i++)
    {
        p_lp[i] = w0*p[i];
    }
    // Each shifted pole decreases gain by w0, each shifted zero increases
    // it.  Cancel out the net change to keep overall gain the same
    double k = zpkIn.getGain(); 
    int ndeg = npoles - nzeros;
    double k_lp = k*std::pow(w0, ndeg);
    zpkOut = ZPK(z_lp, p_lp, k_lp);
    return 0;
}
/*!
 * @brief Converts a lowpass filter prototype to a highpss filter.
 * @param[in] zpkIn    Input lowpass filter prototype to convert.
 * @param[in] w0       Desired cutoff frequency (rad/s).
 * @param[out] zpkOut  The corresponding highpass filter.
 * @result 0 indicates success. 
 * @ingroup rtseis_utils_design
 */
int IIR::zpklp2hp(const ZPK zpkIn, const double w0, ZPK &zpkOut)
{
    zpkOut.clear();
    size_t nzeros = zpkIn.getNumberOfZeros();
    size_t npoles = zpkIn.getNumberOfPoles();
    if (w0 < 0)
    {
        RTSEIS_ERRMSG("w0=%lf must be non-negative", w0);
        return -1;
    }
    if (npoles < nzeros)
    {
        RTSEIS_ERRMSG("%s", "npoles < nzeros not implemented");
        return -1;
    }
    std::vector<std::complex<double>> z = zpkIn.getZeros();
    std::vector<std::complex<double>> p = zpkIn.getPoles();
    // Invert positions radially about unit circle to convert LPF TO HPF
    // Scale all points radially from origin to shift cutoff frequency
    std::vector<std::complex<double>> z_hp(std::max(nzeros, npoles));
    std::complex<double> zprod = std::complex<double> (1, 0);
    for (size_t i=0; i<z.size(); i++)
    {
        z_hp[i] = w0/z[i];
        zprod = (-z[i])*zprod;
    }
    // If lowpass had zeros at infinity, inverting moves them to origin
    for (size_t i=nzeros; i<npoles; i++)
    {
        z_hp[i] = std::complex<double> (0, 0);
    }
    // Invert positions radially about unit circle to convert LPF TO HPF
    // Scale all points radially from origin to shift cutoff frequency
    std::vector<std::complex<double>> p_hp(npoles);
    std::complex<double> pprod = std::complex<double> (1, 0);
    for (size_t i=0; i<p.size(); i++)
    {
        p_hp[i] = w0/p[i];
        pprod = (-p[i])*pprod;
    }
    // Compute scale factors and cancel out gain caused by inversion
    double k = zpkIn.getGain();
    std::complex<double> zt = zprod/pprod;
    double k_hp = k*std::real(zt);
    zpkOut = ZPK(z_hp, p_hp, k_hp);
    return 0;
}
//============================================================================//
//                            Utility functions for zpk2sos                   //
//============================================================================//
/*!
 * @brief Utility function for getting index of next closest pole to
 *        unit circle.
 */
static int getNextWorstRealPoleIndex(const std::vector<bool> pmask,
                                     const std::vector<bool> isreal_pole,
                                     const std::vector<std::complex<double>> p)
{
    const std::complex<double> zone = std::complex<double> (1, 0);
    double difMin = DBL_MAX;
    int indx =-1;
    for (size_t i=0; i<p.size(); i++)
    {
        if (!pmask[i] && isreal_pole[i])
        {
            double dif = std::abs(p[i] - zone);
            if (dif < difMin)
            {
                dif = difMin;
                indx = i;
            }
        }
    }
    return indx;
}
/*!
 * @brief Utility function for getting next pole closest to unit circle.
 */
static int getWorstPoleIndex(const std::vector<bool> pmask,
                             const std::vector<std::complex<double>> p)
{
    double difMin = DBL_MAX;
    int indx =-1;
    for (int i=0; i<pmask.size(); i++)
    {
        if (!pmask[i])
        {
            double dif = std::abs(1.0 - std::abs(p[i]));
            if (dif < difMin)
            {
                difMin = dif;
                indx = i;
            }
        }
    }
    return indx;
}
/*!
 * @brief Utility function for computing the sum of the reals of the poles
 *        array
 */
static int isreal_p_sum(const std::vector<bool> pmask,
                        const std::vector<bool> isreal_pole)
{
    int xsum = 0;
    for (size_t i=0; i<isreal_pole.size(); i++)
    {
        if (!pmask[i] && isreal_pole[i]){xsum = xsum + 1;}
    }
    return xsum;
}
/*!
 * @brief Utility function for finding the closest real or complex zero
 *        to the pole p1
 */
static int nearest_real_complex_idx(const std::vector<bool> zmask,
                                    const std::vector<bool> isreal_zero,
                                    const std::vector<std::complex<double>> z,
                                    const std::complex<double> p1,
                                    const bool lreal)
{
    double difMin = DBL_MAX;
    int indx =-1;
    // Look for the closest real zero to the pole p1
    if (lreal)
    {
        for (size_t i=0; i<z.size(); i++)
        {
            if (isreal_zero[i] && !zmask[i])
            {
                double dif = std::abs(z[i] - p1);
                if (dif < difMin)
                {
                    difMin = dif;
                    indx = i;
                }
            }
        }
    }
    // Look for the closest complex zero to the pole p1
    else
    {
        for (size_t i=0; i<z.size(); i++)
        {
            if (!isreal_zero[i] && !zmask[i])
            {
                double dif = std::abs(z[i] - p1);
                if (dif < difMin)
                {
                    difMin = dif;
                    indx = i;
                }
            }
        }
    }
    return indx;
}

/*!
 * @brief Find the index of the next nearest real or complex zero to the
 *        given pole.
 */
static int nearest_realOrComplex_idx(const std::vector<bool> zmask,
                                     const std::vector<std::complex<double>> z,
                                     const std::complex<double> p1)
{
    double difMin = DBL_MAX;
    int indx =-1;
    for (size_t i=0; i<z.size(); i++)
    {
        if (!zmask[i])
        {
            double dif = std::abs(p1 - z[i]);
            if (dif < difMin)
            {
                difMin = dif;
                indx = i;
            }
        }
    }
    return indx;
}
static int cmplxreal(const std::vector<std::complex<double>> z,
                     std::vector<std::complex<double>> &p,
                     std::vector<bool> &isreal,
                     const double *tolIn)
{
    p.resize(0);
    isreal.resize(0);
    if (z.size() < 1)
    {
        RTSEIS_ERRMSG("%s", "No elements in z");
        return -1;
    }
    // Set the tolerance
    double tol = 100.0*DBL_EPSILON;
    if (tolIn != nullptr){tol = *tolIn;}
    // Get the real poles and complex poles 
    std::vector<double> reals;
    std::vector<std::complex<double>> cmplxs;
    for (size_t i=0; i<z.size(); i++)
    {
        if (std::imag(z[i]) == 0)
        {
            reals.push_back(std::real(z[i]));
        }
        else
        {
            cmplxs.push_back(z[i]);
        }
    }
    size_t nreal  = reals.size();
    size_t ncmplx = cmplxs.size();
    if (nreal + ncmplx != z.size())
    {
        RTSEIS_ERRMSG("%s", "Algorithmic error");
        return -1;
    }
    // Sort the reals into ascending order and save them 
    nreal = 0;
    if (reals.size() > 0)
    {
        std::sort(reals.begin(), reals.end());
        for (size_t i=0; i<reals.size(); i++)
        {
            p.push_back(std::complex<double> (reals[i], 0));
            isreal.push_back(true);
            nreal = nreal + 1;
        }
    }
    // Pair the complex conjugate poles up and save them
    size_t nconj = 0;
    if (cmplxs.size() > 0)
    {
        // Sort complex numbers based on magnitudes
        std::sort(cmplxs.begin(), cmplxs.end(),
                  [](std::complex<double> a,  
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
        // Verify everyone has a buddy (conjugate pair)
        std::vector<int> lskip(cmplxs.size(), 0);
        for (size_t i=0; i<cmplxs.size(); i++)
        {
            if (lskip[i] == 1){continue;}
            // Find it's complex conjugate
            for (size_t j=i+1; j<cmplxs.size(); j++)
            {
                if (lskip[j] == 1){continue;}
                if (std::abs(cmplxs[i] - std::conj(cmplxs[j])) < tol)
                {
                    // Balance poles for numerical accuracy
                    std::complex<double> zhalf
                         = (cmplxs[i] + std::conj(cmplxs[j]))/2.0;
                    // Retain the `positive' pole for consistency 
                    if (std::imag(zhalf) >= 0)
                    {
                        p.push_back(zhalf); 
                        //p.push_back(std::conj(zhalf));
                    }
                    else
                    {
                        p.push_back(std::conj(zhalf)); 
                        //p.push_back(zhalf);
                    }
                    nconj = nconj + 1;
                    isreal.push_back(false);
                    //isreal.push_back(false);
                    lskip[i] = 1;
                    lskip[j] = 1;
                }
            }
        }
        // I'm not dealing with unpaired
        int npaired = std::accumulate(lskip.begin(), lskip.end(), 0);
        if (npaired != static_cast<int> (cmplxs.size()))
        {
            RTSEIS_ERRMSG("All poles need complex conjugate pairs; %d %ld",
                          npaired, cmplxs.size());
            p.resize(0);
            return -1;
        }
    }
    if (2*nconj + nreal != z.size())
    {
        RTSEIS_ERRMSG("%s", "Algorithmic failure");
        return -1;
    }
/*
printf("Input\n");
for (size_t i=0; i<z.size(); i++)
{
 printf("%+lf%+lfi\n", std::real(z[i]), std::imag(z[i]));
}
printf("Output\n");
for (size_t i=0; i<p.size(); i++)
{
 printf("%+lf%+lfi\n", std::real(p[i]), std::imag(p[i]));
}
*/
    return 0;
}
