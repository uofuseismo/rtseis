#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include <vector>
#include <cmath>
#include "rtseis/utils/design.hpp"
#include "rtseis/utils/polynomial.hpp"
#include "rtseis/utils/ipps.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utils::FilterDesign;
using namespace RTSeis::Utils::Math;

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
 * @param[in] lanlaog  If true then this designs an analog filter.  The
 *                     default is a digital filter.
 * @param[out] ba      The corresponding filter specified in terms of numerators
 *                     and denominators.
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
 * @param[in] lanlaog  If true then this designs an analog filter.  The
 *                     default is a digital filter.
 * @param[out] zpk     The corresonding filter stored in zero, pole, and gain
 *                     format.
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
