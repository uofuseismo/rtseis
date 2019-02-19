#ifndef RTSEIS_UTILS_DESIGN_RESPONSE_HPP
#define RTSEIS_UTILS_DESIGN_RESPONSE_HPP 1
#include <complex>
#include <vector>
#include "rtseis/config.h"

namespace RTSeis
{
class BA;
namespace Utilities
{
namespace FilterDesign
{
/*!
 * @defgroup rtseis_utils_design_response Response
 * @brief Utility functions for analyzing the response of a filter.
 *        This code is originally from ISTI's ISCL and has been modified
 *        to conform with C++.  Function names have also been changed to
 *        conform with rtseis's naming conventions.
 * @copyright ISTI distributed under the Apache 2 license.
 * @ingroup rtseis_utils_design
 */
namespace Response
{
    /*!
     * @brief Computes the complex frequency response H(s) of an analog filter
     *
     *        \f[
     *           H(s)
     *         = \frac{B(s)}{A(s)}
     *         = \frac{b[0]s^{n_b} + b[1]s^{n_b-1} + ... + b[nb-1]s + b[nb-1]}
     *                {a[0]s^{n_a} + a[1]s^{n_a-1} + ... + a[na-1]s + a[na-1]}
     *        \f]
     *
     *        given the numerator ande denomiantor coefficients at
     *        angular frequencies, w.
     *
     * @param[in] ba     The transfer function defining the analog filter.
     * @param[in] w      The angular frequencies (rad/s) at which to
     *                   tabulate the response.
     * @param[out] h     The frequency repsonse, \f$ H(i \omega) \f$, tabulated 
     *                   at the angular frequencies.  This has dimension
     *                   [w.size()].
     * @result 0 indicates success.
     * @ingroup rtseis_utils_design_response
     */
    int freqs(const BA &ba, const std::vector<double> &w,
              std::vector<std::complex<double>> &h);
    /*!
     * @brief Computes the complex frequency response H(z) of a digital filter
     *
     *        \f[
     *           H(z)
     *         = \frac{B(z)}{A(z)}
     *         = \frac{b[0] + b[1]z^{-1} + b[2]z^{-2} + ... + b[nb-1]z^{-n_b}}
     *                {a[0] + a[1]z^{-1} + a[2]z^{-2} + ... + a[na-1]z^{-n_a}}
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
    int freqz(const BA &ba, const std::vector<double> &w,
              std::vector<std::complex<double>> &h);
    int freqz(const BA &ba, const std::vector<double> &w,
              std::vector<std::complex<double>> &h);
    /*!
     * @brief Computes the group delay of a filter.  The group delay
     *        is a measure of the average delay of the filter as a function
     *        of frequency.  
     *        \f[
     *            \tau_g(\omega) =-\frac{d \theta(\omega}}{d \omega}
     *        \f$.
     *        where \f$ \theta(\omega) \f$ is the argument of the transfer
     *        function \f$ H(z) = H(e^{i\omega}} \f$.
     * @param[in] ba   The transfer function defining the digital filter.
     * @param[in] w    The normalized frequencies at which to evaluate the group
     *                 delay.  In this case, 0 is the zero frequency and
     *                 1 is the Nyquist frequency.
     * @param[out] gd  The group delay in samples.  This will have dimension
     *                 [w.size()].
     * @result 0 indicates success.
     */
     int groupDelay(const BA &ba, const std::vector<double> &w,
                    std::vector<double> &gd);


}; /* End Response */

}; /* End FilterDesign */
}; /* End Utils */

}; /* End RTseis */

#endif
