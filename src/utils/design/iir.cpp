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
 * @brief Computes the polynomial transfer function from a pole-zero
 *        representation.
 * @param[in] zpk  Input transfer function specified in terms of zeros,
 *                 poles, and gain.
 * @param[out] ba  Corresponding transfer function specified in terms of
 *                 numerator and denominator coefficients.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
/*
int Design::zpk2tf(const ZPK zpk, BA &ba)
{
    double k = zpk.getGain();
    if (k == 0){RTSEIS_WARNMSG("%s", "System gain is zero");}
    size_t nzeros = zpk.getNumberOfZeros();
    size_t npoles = zpk.getNumberOfPoles();
    size_t nb = nzeros + 1;
    size_t na = npoles + 1; 
    std::vector<std::complex<double>> z = zpk.getZeros();
    
    return 0;
}
*/
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
    for (size_t i=0; i<nzeros; i++)
    {
        z_lp[i] = w0*z[i];
    }
    std::vector<std::complex<double>> p_lp(npoles);
    for (size_t i=0; i<npoles; i++)
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
