#ifndef RTSEIS_UTILS_POLYNOMIAL_HPP
#define RTSEIS_UTILS_POLYNOMIAL_HPP 1
#include <complex>
#include <vector>

namespace RTSeis
{
/*!
 * @defgroup rtseis_utils Utilities
 * @brief Herein lies much of the primitives on which the high-level
 *        processing modules are built.  This includes basic mathematical
 *        operations, filter design, filter implementations, and transform
 *        implementations. 
 * @copyright Ben Baker distributed under the MIT license.
 */
namespace Utilities
{
/*!
 * @defgroup rtseis_utils_math Math Utilities
 * @brief Utility functions for performing common mathematical operations.
 * @ingroup rtseis_utils
 */
namespace Math 
{

/*!
 * @defgroup rtseis_utils_math_polynomial Polynomial
 * @brief Utility functions for polynomial handling.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_math
 */
namespace Polynomial
{
    /*!
     * @brief Computes the roots of a polynomial:
     *        \f[
     *             q(x) = c_0 x^p + c_1 x^{p-1} + \cdots + c_{p+1}
     *        \f]
     *        where \f$ p \f$ is the polynomial order
     *        and \f$ coeffs[0] = c_0 \f$.
     * @param[in] coeffs   The coefficients of the polynomial whose order 
     *                     is defined above.  Note, the coeffs[0] cannot be
     *                     0 and coeffs must have length at least 2.
     * @param[out] roots   The roots of the polynomial.  This has dimension
     *                     [coeffs.size() - 1].
     * @result 0 indicates success.
     * @ingroup rtseis_utils_math_polynomial
     */
    int roots(const std::vector<double> &coeffs,
              std::vector<std::complex<double>> &roots);
    /*!
     * @{
     * @brief Returns a polynomial whose roots are given by p.
     * @param[in] p    The polynomial roots.  
     * @param[out] y   The polynomial coefficients corresponding to the
     *                 roots of the given polynomial.  This has dimension
     *                 [p.size()+1] and is ordered so that the last coefficient
     *                 is the constant term and the first coefficient scales
     *                 the highest order polynomial.
     * @result 0 indicates success.
     * @ingroup rtseis_utils_math_polynomial
     */
    int poly(const std::vector<std::complex<double>> &p,
             std::vector<std::complex<double>> &y);
    /*!
     * @copydoc poly
     * @ingroup rtseis_utils_math_polynomial
     */
    int poly(const std::vector<double> &p,
             std::vector<double> &y);
    /*! @} */
    /*!
     * @{
     * @brief Evaluates the polynomial
     *        \f[
     *            p(x) = p_{n_{order}}
     *                 + x p_{n_{order}-1}
     *                 + \cdots
     *                 + x^{n_{order}} p_0
     *        \f]
     *        at points \f$ x_j, j=1,2,...,n_x \f$.
     * @param[in] p    The polynomial coefficients ordered such that the
     *                 highest order coefficient comes first.  This has
     *                 dimension [order+1].
     * @param[in] x    The points at which to evaluate the polynomial.  This
     *                 has dimension [x.size()].
     * @param[out] y   \f$ y = p(x) \f$ evaluated at each \f$ x_i \f$.  This
     *                 has dimension [x.size()].
     * @result 0 indicates success.
     * @ingroup rtseis_utils_math_polynomial 
     */
    int polyval(const std::vector<double> &p,
                const std::vector<double> &x,
                std::vector<double> &y);
    /*!
     * @copydoc polyval
     * @ingroup rtseis_utils_math_polynomial
     */
    int polyval(const std::vector<std::complex<double>> &p,
                const std::vector<std::complex<double>> &x,
                std::vector<std::complex<double>> &y);
    /*! @} */
}; /* End polynomial. */

}; /* End math. */
}; /* End utils. */
}; /* End rtseis. */

#endif
