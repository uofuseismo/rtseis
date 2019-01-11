#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include <vector>
#include <cmath>
#include "rtseis/utils/design.hpp"
#include "rtseis/utils/polynomial.hpp"
#include "rtseis/utils/ipps.hpp"
#include "rtseis/log.h"
#include <ipps.h>

using namespace RTSeis::Utils::FilterDesign;

/*!
 * @defgroup rtseis_utils_design_response Response
 * @brief Utility functions for looking at the response of a filter.
 *        This code is originally from ISTI's ISCL and has been modified
 *        to conform with C++.  Function names have also been changed to
 *        conform with rtseis's naming conventions.
 * @copyright ISTI distributed under the Apache 2 license.
 * @ingroup rtseis_utils_design
 */

/*!
 * @brief Tabulates the frequency response of an analog filter.
 * @param[in] ba     The transfer function defining the analog filter.
 * @param[in] w      The angular frequencies (rad/s) at which to tabulate the
 *                   response.
 * @param[out] h     The frequency repsonse, \f$ H(i \omega) \f$, tabulated at
 *                   the angular frequencies.  This has dimension [w.size()].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design_response
 */
int Response::freqs(const BA ba, const std::vector<double> w,
                    std::vector<std::complex<double>> &h)
{
    size_t nw = w.size();
    h.resize(0);
    if (nw == 0){return 0;}
    // Get a handle on the numerator and denominator coefficients
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients();
    bool lzero = true;
    for (size_t ia=0; ia<a.size(); ia++)
    {
        if (a[ia] != 0)
        {
            lzero = false;
            break;
        }
    } 
    if (lzero)
    {
        RTSEIS_ERRMSG("%s", "a is entirely 0; division by zero");
        return -1;
    }
    // Compute s = i \omega
    std::vector<std::complex<double>> s;
    s.resize(nw);
    for (size_t i=0; i<nw; i++)
    {
        s[i] = std::complex<double> (0, w[i]); // s = i omega
    }
    // Evaluate the numerator and denominator polynoimals
    std::vector<std::complex<double>> bz;
    bz.resize(b.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<b.size(); i++){bz[i] = std::complex<double> (b[i], 0);}
    std::vector<std::complex<double>> hsNum;
    int ierr = Math::Polynomial::polyval(bz, s, hsNum);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hsNum");
        return -1;
    }
    std::vector<std::complex<double>> hsDen;
    std::vector<std::complex<double>> az;
    az.resize(a.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<a.size(); i++){az[i] = std::complex<double> (a[i], 0);}
    ierr = Math::Polynomial::polyval(az, s, hsDen);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hsDen");
        return -1;
    }
    // Compute the transfer function
    ierr = IPPS::Div(hsDen, hsNum, h);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Division failed");
        return -1;
    }
    return 0;
}
/*!
 * @brief Computes the complex frequency response H(z) of a digital filter
 *
 *        \f[
 *           H(z)
 *         = \frac{B(z)}{A(z)}
 *         = \frac{b[0] + b[1]z^{-1} + b[2]z^{-2} + ... + b[nb-1]z^{-nb}}
 *                {a[0] + a[1]z^{-1} + a[2]z^{-2} + ... + a[na-1]z^{-na}}
 *        \f]
 * 
 *        given the numerator and denominator coefficients at normalized
 *        angular frequencies, w.
 *
 * @param[in] ba   The transfer function defining the digital filter.
 * @param[in] w    The normalized angular frequencies at which to evaluate
 *                 the transfer function.  Values should be in the range
 *                 \f$ \omega \in [0, \pi] \f$ where \f$ 0 \f$ is the zero
 *                 frequency and \f$ \pi \f$ is the Nyquist frequency.
 * @param[out] h   The frequency repsonse, \f$ H(z) \f$, tabulated at the
 *                 normalized angular frequencies.  This has dimension
 *                 [w.size()].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design_response
 */
int Response::freqz(const BA ba, const std::vector<double> w,
                    std::vector<std::complex<double>> &h) 
{
    size_t nw = w.size();
    h.resize(0);
    if (nw == 0){return 0;} 
    // Get a handle on the numerator and denominator coefficients
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients();
    bool lzero = true;
    for (size_t ia=0; ia<a.size(); ia++)
    {   
        if (a[ia] != 0)
        {
            lzero = false;
            break;
        }
    }   
    if (lzero)
    {   
        RTSEIS_ERRMSG("%s", "a is entirely 0; division by zero");
        return -1; 
    }
    // Copy numerators and denominators in reverse order for consistency with
    // polyval and then compute z.  The polyval convention requires us to
    // compute z=e^{-i\omega} to compensate.  Normally, this is would
    // be z^{i \omega}.
    size_t nb = b.size();
    std::vector<std::complex<double>> bz;
    bz.resize(nb);
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<nb; i++)
    {
        bz[i] = b[nb-1-i];
    }
    size_t na = a.size();
    std::vector<std::complex<double>> az;
    az.resize(a.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<na; i++)
    {
        az[i] = a[na-1-i];
    }
    std::vector<std::complex<double>> z;
    z.resize(nw);
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<nw; i++)
    {
        z[i] = std::exp(std::complex<double> (0, -w[i]));
    }
    // Evaluate the numerator and denominator polynoimals
    std::vector<std::complex<double>> hzNum;
    int ierr = Math::Polynomial::polyval(bz, z, hzNum);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hzNum");
        return -1; 
    }
    std::vector<std::complex<double>> hzDen;
    ierr = Math::Polynomial::polyval(az, z, hzDen);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hzDen");
        return -1;
    }
    // Compute the transfer function: H = Num/Den
    ierr = IPPS::Div(hzDen, hzNum, h); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Division failed");
        return -1;
    }
    return 0;
}
