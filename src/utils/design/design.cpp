#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include <vector>
#include <cmath>
#include "rtseis/utils/design.h"
#include "rtseis/utils/polynomial.h"
#include "rtseis/log.h"


/*!
 * @defgroup rtseis_utils_design Design
 * @brief Utility functions for filter design.  This code is originally
 *        from ISTI's ISCL and has been modified to conform with C++.
 *        Function names have also been changed to conform with rtseis's
 *        naming conventions.
 * @copyright ISTI distributed under the Apache 2 license.
 * @ingroup rtseis_utils_design
 */
Design::Design(void)
{
    return;
}

Design::~Design(void)
{
    return;
}
/*!
 * @brief Tabulates the frequency response of an analog filter.
 * @param[in] ba     The transfer function defining the analog filter.
 * @param[in] w      The angular frequencies (rad/s) at which to tabulate the
 *                   response.
 * @param[out] h     The frequency repsonse, \f$ H(i \omega) \f$, tabulated at
 *                   the angular frequencies.  This has dimension [w.size()].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design
 */
int Design::freqs(const BA ba, const std::vector<double> w,
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
    #pragma ivdep
    for (size_t i=0; i<b.size(); i++){bz[i] = std::complex<double> (b[i], 0);}
    std::vector<std::complex<double>> hsNum;
    Polynomial polynomial;
    int ierr = polynomial.polyval(bz, s, hsNum);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hsNum");
        return -1;
    }
    std::vector<std::complex<double>> hsDen;
    std::vector<std::complex<double>> az;
    az.resize(a.size());
    #pragma ivdep
    for (size_t i=0; i<a.size(); i++){az[i] = std::complex<double> (a[i], 0);}
    ierr = polynomial.polyval(az, s, hsDen);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hsDen");
        return -1;
    }
    // Compute the transfer function
    h.resize(nw);
    #pragma ivdep
    for (size_t i=0; i<nw; i++)
    {
        h[i] = hsNum[i]/hsDen[i];
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
 * @ingroup rtseis_utils_design
 */
int Design::freqz(const BA ba, const std::vector<double> w,
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
    #pragma ivdep
    for (size_t i=0; i<nb; i++)
    {
        bz[i] = b[nb-1-i];
    }
    size_t na = a.size();
    std::vector<std::complex<double>> az;
    az.resize(a.size());
    #pragma ivdep
    for (size_t i=0; i<na; i++)
    {
        az[i] = a[na-1-i];
    }
    std::vector<std::complex<double>> z;
    z.resize(nw);
    #pragma ivdep
    for (size_t i=0; i<nw; i++)
    {
        z[i] = std::exp(std::complex<double> (0, -w[i]));
    }
    // Evaluate the numerator and denominator polynoimals
    std::vector<std::complex<double>> hzNum;
    Polynomial polynomial;
    int ierr = polynomial.polyval(bz, z, hzNum);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hzNum");
        return -1; 
    }
    std::vector<std::complex<double>> hzDen;
    ierr = polynomial.polyval(az, z, hzDen);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hzDen");
        return -1;
    }
    // Compute the transfer function
    h.resize(nw);
    #pragma ivdep
    for (size_t i=0; i<nw; i++)
    {
        h[i] = hzNum[i]/hzDen[i];
    }
    return 0;
}
