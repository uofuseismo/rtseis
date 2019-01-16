#ifndef RTSEIS_UTILS_POLYNOMIAL_HPP
#define RTSEIS_UTILS_POLYNOMIAL_HPP 1
#include "rtseis/config.h"
#include <complex>
#include <vector>

namespace RTSeis
{
namespace Utils
{
namespace Math 
{

namespace Polynomial
{
    int roots(const std::vector<double> coeffs,
              std::vector<std::complex<double>> &roots);
    int poly(const std::vector<std::complex<double>> p,
             std::vector<std::complex<double>> &y);
    int poly(const std::vector<double> p, std::vector<double> &y);
    int polyval(const std::vector<double> p,
                const std::vector<double> x,
                std::vector<double> &y);
    int polyval(const std::vector<std::complex<double>> p,
                const std::vector<std::complex<double>> x,
                std::vector<std::complex<double>> &y);
}; /* End polynomial. */

}; /* End math. */
}; /* End utils. */
}; /* End rtseis. */

#endif
