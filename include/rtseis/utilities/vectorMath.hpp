#ifndef RTSEIS_UTILITIES_VECTORMATH_HPP
#define RTSEIS_UTILITIES_VECTORMATH_HPP 1
#include <complex>
#include <vector>
#include "rtseis/config.h"

namespace RTSeis
{
namespace Utilities
{
namespace Math
{
namespace Vector
{

/*!
 * @brief Divides element-wise res = num/den.
 * @param[in] den   The array of denominators.
 * @param[in] num   The array of numerators.
 * @param[out] res  The result, res=num/den.
 * @result 0 indicates success.
 */
int divide(const std::vector<std::complex<double>> &den,
           const std::vector<std::complex<double>> &num,
           std::vector<std::complex<double>> &res);

/*!
 * @brief Computes the real part of a complex vector.
 * @param[in] z   The complex array of which to take the real.
 * @param[out] r  The real part of the z.
 * @result 0 indicates success.
 */
int real(const std::vector<std::complex<double>> &z,
         std::vector<double> &r);


};
}; // End Math
}; // End Utilities
}; // End RTSeis

#endif
