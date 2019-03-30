#ifndef RTSEIS_UTILITIES_VECTORMATH_HPP
#define RTSEIS_UTILITIES_VECTORMATH_HPP 1
#include <complex>
#include <vector>

namespace RTSeis
{
namespace Utilities
{
namespace Math
{
/*!
 * @defgroup rtseis_utils_math_vm Vector Math
 * @brief Utilities for vectorized math operations.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_math
 */
namespace VectorMath
{

/*!
 * @brief Divides element-wise res = num/den.
 * @param[in] den   The array of denominators.
 * @param[in] num   The array of numerators.
 * @param[out] res  The result, res=num/den.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_math_vm
 */
int divide(const std::vector<std::complex<double>> &den,
           const std::vector<std::complex<double>> &num,
           std::vector<std::complex<double>> &res);

/*!
 * @brief Computes the real part of a complex vector.
 * @param[in] z   The complex array of which to take the real.
 * @param[out] r  The real part of the z.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_math_vm
 */
int real(const std::vector<std::complex<double>> &z,
         std::vector<double> &r);

/*!
 * @defgroup rtseis_utils_math_vm_copysign Copysign
 * @brief Computes the copysign of a float or double precision array s.t.
 *        \f[
 *           y =
 *           \left \{
 *             \begin{array}{rr}
 *               -1 & x \le -0 \\
 *               +1 & x \ge +0 \\
 *             \end{array}     
 *           \right .
 *        \f]
 * @param[in] n   The length of the input array.
 * @param[in] x   The array of which to compute the sign bit.
 *                This has dimension [n].
 * @param[out] y  The sign bit of the input array.  This has dimension [n].
 * @ingroup rtseis_utils_math_vm
 */
template<typename T> int copysign(const int n, const T x[], T y[]);
/*!
 * @brief Computes the copysign of a float or double precision array s.t.
 *        \f[
 *           y =
 *           \left \{
 *             \begin{array}{rr}
 *               -1 & x \le -0 \\
 *               +1 & x \ge +0 \\
 *             \end{array}     
 *           \right .
 *        \f]
 * @param[in] x   The array of which to compute the sign bit.
 * @param[out] y  The sign bit of the input array.  This has
 *                dimension [x.size()].
 * @ingroup rtseis_utils_math_vm_copysign
 */
template<typename T> int copysign(const std::vector<T> &x, std::vector<T> &y);

/*!
 * @defgroup rtseis_utilties_math_vm_isSorted isSorted
 * @ingroup rtseis_utils_math_vm
 * @brief Determines if the vector is sorted in increasing order.
 * @param[in] x  The array to determine if is sorted or not.
 * @retval True indicates that x is sorted in increasing order.
 * @retval False indicates that x is not sorted in increasing order.
 * @ingroup rtseis_utils_math_vm
 */
template<typename T> bool isSorted(const std::vector<T> &x);

};
}; // End Math
}; // End Utilities
}; // End RTSeis

#endif
