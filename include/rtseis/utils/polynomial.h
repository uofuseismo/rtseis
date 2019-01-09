#ifndef RTSEIS_UTILS_ZPK_HPP
#define RTSEIS_UTILS_ZPK_HPP 1
#include "rtseis/config.h"
#include <stdio.h>
#include <complex>
#include <vector>

#ifdef __cplusplus
class Polynomial 
{
    public:
        Polynomial(void);
        ~Polynomial(void); 
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
};
#endif

#endif
