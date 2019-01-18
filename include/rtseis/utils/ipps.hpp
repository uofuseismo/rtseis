#ifndef RTSEIS_UTILS_IPPS_HPP
#define RTSEIS_UTILS_IPPS_HPP 1
#include "rtseis/config.h"
#include "rtseis/enums.h"
#include <complex>
#include <vector>
#include <memory>

namespace RTSeis
{
namespace Utils
{
namespace IPPS
{
    /* Divides res = num/den. */
    int Div(std::vector<std::complex<double>> den,
            std::vector<std::complex<double>> num,
            std::vector<std::complex<double>> &res);
    /* Multiplies x = val*x. */
    int MulC(const std::complex<double> val,
             std::vector<std::complex<double>> x);

};
};
};

#endif
