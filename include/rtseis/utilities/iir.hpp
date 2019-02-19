#ifndef RTSEIS_UTILS_DESIGN_IIR_HPP
#define RTSEIS_UTILS_DESIGN_IIR_HPP 1
#include <complex>
#include <vector>
#include "rtseis/config.h"

namespace RTSeis
{
class SOS;
class BA;
class ZPK;
namespace Utilities
{
namespace FilterDesign
{

/*!
 * @defgroup rtseis_utils_design_iir IIR Design
 * @brief Utility functions for IIR filter design.  This code is originally
 *        from ISTI's ISCL and has been modified to conform with C++.
 *        Function names have also been changed to conform with rtseis's
 *        naming conventions.
 * @copyright ISTI distributed under the Apache 2 license.
 * @ingroup rtseis_utils_design
 */
namespace IIR
{
    /*!
     * @brief Defines the analog prototype for the IIR filter design.
     * @ingroup rtseis_utils_design_iir
     */
    enum Prototype
    {
        BUTTERWORTH, /*!< Butterworth filter design. */
        BESSEL,      /*!< Bessel filter design. */
        CHEBYSHEV1,  /*!< Chebyshev I filter design. */
        CHEBYSHEV2   /*!< Chebyshev II filter design. */
    };
    /*!
     * @brief Defines the bandtype for the IIR filter design.
     * @ingroup rtseis_utils_design_iir
     */
    enum Bandtype
    {
        LOWPASS,  /*!< Lowpass filter. */
        HIGHPASS, /*!< Highpass filter. */
        BANDPASS, /*!< Bandpass filter. */
        BANDSTOP  /*!< Bandstop filter. */
    };
    /*! 
     * @brief Pairing strategy when converting a ZPK filter to an 
     *        SOS filter.
     * @ingroup rtseis_utils_design_iir
     */
    enum Pairing
    {
        NEAREST = 0,  /*!< This attempts to minimize the peak gain. */
        KEEP_ODD = 1  /*!< This attempts to minimize the peak gain
                           subject to the constraint that odd-order
                           systems should retain one section as first order. */
    };
    /*!
     * @brief Convenience function to design a digital or analog filter from
     *        an analog prototype.
     * @param[in] n        Order of filter to design.
     * @param[in] W        A scalar or length-2 array defining the critical
     *                     frequencies.  For digital filters, these are
     *                     normalized frequencies in the range [0,1] where 1
     *                     is the Nyquist frequency in pi radians/sample
     *                     (thereby making Wn half-cycles per sample.)   For
     *                     analog filters, Wn is the angular frequency in rad/s.
     * @param[in] rp       For Chebyshev I filters this specifies the maximum
     *                     ripple in the passband in decibels.
     * @param[in] rs       For Chebyshev II filters this specifies the minimum
     *                     attenutation in the stop band in decibels.
     * @param[in] btype    The type of filter, e.g., lowpass, highpass,
     *                     bandpass, or bandstop.
     * @param[in] ftype    The IIR filter protoytpe, e.g., Butterworth, Bessel,
     *                     Chebyshev1, or Chebyshev2.
     * @param[out] ba      The corresponding IIR filter specified in terms
     *                     of numerator and denominator coefficients.
     * @param[in] lanlaog  If true then this designs an analog filter.  The
     *                     default is a digital filter.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int iirfilter(const int n, const double *W,
                  const double rp, const double rs,
                  const Bandtype btype,
                  const Prototype ftype,
                  BA &ba,
                  const bool lanalog = false);
    /*!
     * @brief Convenience function to design a digital or analog filter from
     *        an analog prototype.
     * @param[in] n        Order of filter to design.
     * @param[in] W        A scalar or length-2 array defining the critical
     *                     frequencies.  For digital filters, these are
     *                     normalized frequencies in the range [0,1] where 1
     *                     is the Nyquist frequency in pi radians/sample
     *                     (thereby making Wn half-cycles per sample.)   For
     *                     analog filters, Wn is the angular frequency in rad/s.
     * @param[in] rp       For Chebyshev I filters this specifies the maximum
     *                     ripple in the passband in decibels.
     * @param[in] rs       For Chebyshev II filters this specifies the minimum
     *                     attenutation in the stop band in decibels.
     * @param[in] btype    The type of filter, e.g., lowpass, highpass,
     *                     bandpass, or bandstop.
     * @param[in] ftype    The IIR filter protoytpe, e.g., Butterworth, Bessel,
     *                     Chebyshev1, or Chebyshev2.
     * @param[out] zpk     The corresponding IIR filter specified in terms
     *                     of zeros, poles, and gain.
     * @param[in] lanlaog  If true then this designs an analog filter.  The
     *                     default is a digital filter.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int iirfilter(const int n, const double *W,
                  const double rp, const double rs,
                  const Bandtype btype,
                  const Prototype ftype,
                  ZPK &zpk,
                  const bool lanalog = false);
    /*!
     * @brief Convenience function to design a digital or analog filter from
     *        an analog prototype.
     * @param[in] n        Order of filter to design.
     * @param[in] W        A scalar or length-2 array defining the critical
     *                     frequencies.  For digital filters, these are
     *                     normalized frequencies in the range [0,1] where 1
     *                     is the Nyquist frequency in pi radians/sample
     *                     (thereby making Wn half-cycles per sample.)   For
     *                     analog filters, Wn is the angular frequency in rad/s.
     * @param[in] rp       For Chebyshev I filters this specifies the maximum
     *                     ripple in the passband in decibels.
     * @param[in] rs       For Chebyshev II filters this specifies the minimum
     *                     attenutation in the stop band in decibels.
     * @param[in] btype    The type of filter, e.g., lowpass, highpass,
     *                     bandpass, or bandstop.
     * @param[in] ftype    The IIR filter protoytpe, e.g., Butterworth, Bessel,
     *                     Chebyshev1, or Chebyshev2.
     * @param[out] sos     The corresponding IIR filter as a cascaded series of
     *                     second order sections.
     * @param[in] lanlaog  If true then this designs an analog filter.  The
     *                     default is a digital filter.
     * @param[in] pairing  Defines the pairing strategy.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int iirfilter(const int n, const double *W,
                   const double rp, const double rs,
                   const Bandtype btype,
                   const Prototype ftype,
                   SOS &sos,
                   const bool lanalog,
                   const Pairing pairing = Pairing::NEAREST);
    /*!
     * @brief Convert a filter specified as zeros, poles, and a gain to
     *        a filter consisting of cascaded second order sections.
     * @param[in] zpk      ZPK filter to convert to second order sections.
     * @param[out] sos     The corresponding filter stored as cascaded
     *                     section order sections.
     * @param[in] pairint  The pairing strategy.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int zpk2sos(const ZPK &zpk, SOS &sos,
                const Pairing pairing = Pairing::NEAREST);
    /*! 
     * @brief Computes the polynomial transfer function from a pole-zero
     *        representation.
     * @param[in] zpk  Input transfer function specified in terms of zeros,
     *                 poles, and gain.
     * @param[out] ba  Corresponding transfer function specified in terms of
     *                 numerator and denominator coefficients.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int zpk2tf(const ZPK &zpk, BA &ba);
    /*!
     * @brief Converts a lowpass filter prototype to a different cutoff
     *        frequency.
     * @param[in] zpkIn    Input lowpass filter prototype to convert.
     * @param[in] w0       Desired cutoff frequency (rad/s).
     * @param[out] zpkOut  The corresponding lowpass filter.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int zpklp2lp(const ZPK &zpkIn, const double w0, ZPK &zpkOut);
    /*!
     * @brief Converts a lowpass filter prototype to a highpss filter.
     * @param[in] zpkIn    Input lowpass filter prototype to convert.
     * @param[in] w0       Desired cutoff frequency (rad/s).
     * @param[out] zpkOut  The corresponding highpass filter.
     * @result 0 indicates success. 
     * @ingroup rtseis_utils_design_iir
     */
    int zpklp2hp(const ZPK &zpkIn, const double w0, ZPK &zpkOut);
    /*!
     * @brief Transforms a lowpass filter prototype to a bandpass filter.
     *        The passband width is defined by [w0,w1] = [w0,w0+bw] in rad/s.
     * @param[in] zpkIn    Input lowpass filter prototype to convert.
     * @param[in] w0       Desired cutoff frequency in rad/s.
     * @param[in] bw       The desired passband width in rad/s.
     * @param[out] zpkOut  The corresponding bandpass filter.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int zpklp2bp(const ZPK &zpkIn, const double w0,
                 const double bw, ZPK &zpkOut);
    /*!
     * @brief Transforms a lowpass filter prototype to a bandstop filter.
     *        The stopband width is defined by [w0,w1] = [w0,w0+bw] in rad/s.
     * @param[in] zpkIn    Input lowpass filter prototype to convert.
     * @param[in] w0       Desired cutoff frequency in rad/s.
     * @param[in] bw       The desired stopband width in rad/s.
     * @param[out] zpkOut  The corresponding stopband filter.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int zpklp2bs(const ZPK &zpkIn, const double w0,
                 const double bw, ZPK &zpkOut);
    /*!
     * @brief Converts an analog filter to a digital filter using hte bilinear
     *        transform.  This works by transforming a set of poles and zeros
     *        from the s-plane to the digital z-plane using Tustin's method,
     *        which substitutes \f$ s = \frac{z-1}{z+1} \f$ thereby maintaining
     *        the shape of the frequency response.
     *
     * @param[in] zpk     The analog zeros, poles, and gain to transform to the
     *                    z-plane.
     * @param[in] fs      The sampling rate in Hz.  Note, that pre-warping is
     *                    performed in this function.
     * @param[out] zpkbl  The bilinear-transformed variant of zpk.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_iir
     */
    int zpkbilinear(const ZPK zpk, const double fs, ZPK &zpkbl);
}; /* End IIR */

}; /* End FilterDesign */
}; /* End Utils */

}; /* End RTseis */

#endif
