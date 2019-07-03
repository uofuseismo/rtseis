#ifndef RTSEIS_UTILS_WINDOWFUNCTIONS_HPP
#define RTSEIS_UTILS_WINDOWFUNCTIONS_HPP 1
#include <vector>

namespace RTSeis
{
namespace Utilities
{
namespace WindowFunctions
{

/*!
 * @defgroup rtseis_utils_windowFunctions Window Functions
 * @brief Utilities for generating window functions.
 * @{
 */
/*!
 * @brief Creates a Hamming window which is defined as
 *        \f$ 
 *            w_n = 0.54 - 0.46 \cos \left ( \frac{2 \pi n}{L - 1} \right )
 *        \f$
 * @param[in] len      The window length.  This must be positive.
 * @result The Hamming window of length len.
 * @throws std::invalid_argument if len is not positive.
 */
//std::vector<double> hamming(const int len);
//std::vector<float> hamming(const int len);
/*!
 * @brief Creates a Hamming window which is defined as
 *        \f$ 
 *            w_n = 0.54 - 0.46 \cos \left ( \frac{2 \pi n}{L - 1} \right )
 *        \f$
 * @param[in] len      The window length.  This must be positive.
 * @param[out] window  The Hamming window of length len.
 * @throws std::invalid_argument if len is not positive or window is NULL.
 */
void hamming(const int len, double *window[]);
void hamming(const int len, float *window[]);

/*!
 * @brief Creates a Hann window which is defined as
 *        \f$
 *            w_n = 0.5 - 0.5 \cos \left ( \frac{2 \pi n}{L - 1} \right )
 *        \f$
 * @param[in] len      The window length.  This must be positive.
 * @param[out] window  The Hann window of length len.
 * @throws std::invalid_argment if len is not positive or window is NULL.
 */
void hann(const int len, double *window[]);
void hann(const int len, float *window[]);

/*!
 * @brief Creates a Blackman window which is defined as
 *        \f$ 
 *           w_n
 *         = 0.42
 *         - 0.5 \cos \left ( \frac{2\pi n}{L-1} \right )
 *         + 0.08 \cos \left ( \frac{4 \pi n}{L-1} \right )
 *        \f$
 * @param[in] len      The window length.  This must be positive.
 * @param[out] window  The Blackman window of length len.
 * @throws std::invalid_argument if len is not positive or window is NULL.
 */
void blackman(const int len, double *window[]);
void blackman(const int len, float *window[]);

/*!
 * @brief Creates a sine window which is defined as
 *        \f$
 *           w_n = \sin
 *                 \left (
 *                   \frac{\pi n}{L - 1}
 *                 \right )
 *        \f$.
 * @param[in] len      The window length.  This must be positive.
 * @param[out] window  The sine window of length len.
 * @throws std::invalid_argument if len is not postiive or window is NULL.
 */
void sine(const int len, double *window[]);
void sine(const int len, float *window[]);

/*!
 * @brief Creates a Bartlett window which is defined as
 *        \f$
 *           w_n = \left \{
 *                   \begin{array}{lr}
 *                     \frac{2n}{L-1}   & 0 \le n \le \frac{L-1}{2}  \\
 *                     2-\frac{2n}{L-1} & \frac{L-1}{2} < n \le L -1 
 *                   \end{array}
 *                 \right . 
 *        \f$
 * @param[in] len      The window length.  This must be positive.
 * @param[out] window  The Bartlett window of length len.
 * @throw std::invalid_argument if len is not positive or window is NULL.
 */
void bartlett(const int len, double *window[]);
void bartlett(const int len, float *window[]);

/*!
 * @brief Creates a Kaiser window which is defined as
 *        \f$
 *           w_n
 *          =\frac{ I_0 \left (
 *                  \beta \sqrt{ 1 - \left (\frac{n - N/2}{N/2} \right )^2 }
 *                      \right ) }
 *                { I_0 \left ( \beta \right ) } 
 *        \f$ where \f$ N = L - 1 \f$ and \f$ L \f$ is the length.
 * @param[in] len      The window length.  This must be positive.
 * @param[in] beta     An optional parameter that controls the sidelobe 
 *                     attenuation of the Fourier transform window.  The
 *                     default is \f$ \beta = 0.5 \f$.
 * @param[out] window  The Kaiser window of length len.
 * @throws std::invalid_argument if len is negative, beta is negative, or
 *         window is NULL. 
 * @bug For whatever reason IPP does not adequately compute Bessel functions.
 *      If using a C++14 or lesser compiler then the Kaiser window accuracy
 *      will only be valid to 6 digits.  C++17 will compute a high-accuracy
 *      Bessel function using an intrinsic.
 */
void kaiser(const int len, double *window[], const double beta = 0.5);
void kaiser(const int len, float *window[], const float beta = 0.5f);

/*!
 * @}
 */


} /* End WindowFunctions */
} /* End Utilities */
} /* End RTSeis */

#endif
