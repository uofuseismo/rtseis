#ifndef RTSEIS_UTILS_DESIGN_RESPONSE_HPP
#define RTSEIS_UTILS_DESIGN_RESPONSE_HPP 1
#include <complex>
#include <vector>

namespace RTSeis
{
namespace Utilities
{
namespace FilterRepresentations
{
class BA;
};
namespace FilterDesign
{
/*!
 * @defgroup rtseis_utils_design_response Response
 * @brief Utility functions for analyzing the response of a filter.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filterDesign
 */
namespace Response
{
/*!
 * @brief Computes the complex frequency response H(s) of an analog filter
 *        \f[
 *           H(s)
 *         = \frac{B(s)}{A(s)}
 *         = \frac{b[0]s^{n_b} + b[1]s^{n_b-1} + ... + b[nb-1]s + b[nb-1]}
 *                {a[0]s^{n_a} + a[1]s^{n_a-1} + ... + a[na-1]s + a[na-1]}
 *        \f]
 *        given the numerator ande denomiantor coefficients at
 *        angular frequencies, w.
 *
 * @param[in] ba  The transfer function defining the analog filter.
 * @param[in] w   The angular frequencies (rad/s) at which to
 *                tabulate the response.
 * @result The frequency repsonse, \f$ H(i \omega) \f$, tabulated 
 *         at the angular frequencies.  This has dimension [w.size()].
 * @throws std::invalid_argument if there are no numerator or denominator
 *         coefficients or if all of the denominator coefficients are 0.
 * @ingroup rtseis_utils_design_response
 */
std::vector<std::complex<double>>
freqs(const FilterRepresentations::BA &ba,
      const std::vector<double> &w);

/*!
 * @brief Computes the complex frequency response H(z) of a digital filter
 *
 *        \f[
 *           H(z)
 *         = \frac{B(z)}{A(z)}
 *         = \frac{b[0] + b[1]z^{-1} + b[2]z^{-2} + ... + b[nb-1]z^{-n_b}}
 *                {a[0] + a[1]z^{-1} + a[2]z^{-2} + ... + a[na-1]z^{-n_a}}
 *        \f]
 *        given the numerator and denominator coefficients at normalized
 *        angular frequencies, w.
 *
 * @param[in] ba  The transfer function defining the digital filter.
 * @param[in] w   The normalized angular frequencies at which to evaluate
 *                the transfer function.  Values should be in the range
 *                \f$ \omega \in [0, \pi] \f$ where \f$ 0 \f$ is the zero
 *                frequency and \f$ \pi \f$ is the Nyquist frequency.
 * @result The frequency repsonse, \f$ H(z) \f$, tabulated at the
 *         normalized angular frequencies.  This has dimension [w.size()].
 * @throws std::invalid_argument if there are no numerator or denominator
 *         coefficients or if all of the denominator coefficients are 0.
 * @ingroup rtseis_utils_design_response
 */
std::vector<std::complex<double>>
freqz(const FilterRepresentations::BA &ba,
      const std::vector<double> &w);
 
/*!
 * @brief Computes the group delay of a filter.  The group delay
 *        is a measure of the average delay of the filter as a function
 *        of frequency.  
 *        \f[
 *            \tau_g(\omega) =-\frac{d \theta(\omega)}{d \omega}
 *        \f].
 *        where \f$ \theta(\omega) \f$ is the argument of the transfer
 *        function \f$ H(z) = H(e^{i\omega}) \f$.
 * @param[in] ba   The transfer function defining the digital filter.
 * @param[in] w    The normalized frequencies at which to evaluate the group
 *                 delay.  In this case, 0 is the zero frequency and
 *                 1 is the Nyquist frequency.
 * @result The group delay in samples.  This will have dimension [w.size()].
 * @result 0 indicates success.
 * @throws std::invalid_argument if there are no numerator or denominator
 *         coefficients.
 * @ingroup rtseis_utils_design_response
 */
std::vector<double>
groupDelay(const FilterRepresentations::BA &ba,
           const std::vector<double> &w);


}; /* End Response */

}; /* End FilterDesign */
}; /* End Utils */

}; /* End RTseis */

#endif
