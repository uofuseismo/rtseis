#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include "rtseis/utils/ipps.hpp"
#include "rtseis/log.h"
#include <ipps.h>

/*!
 * @defgroup rtseis_utils_ipps IPP C++ Interfaces
 * @brief Generic interfaes to the IPP library.
 * @ingroup rtseis_utils
 */

/*!
 * @brief Divides res = num/den.
 * @param[in] den   The denominator.
 * @param[in] num   The numerator.
 * @param[out] res  The result of num/den.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_ipps
 */
int RTSeis::Utils::IPPS::Div(
    std::vector<std::complex<double>> den,
    std::vector<std::complex<double>> num,
    std::vector<std::complex<double>> &res)
{
    int n = static_cast<int> (den.size());
    res.resize(n);
    Ipp64fc *pSrc1 = static_cast<Ipp64fc *> (static_cast<void *> (den.data()));
    Ipp64fc *pSrc2 = static_cast<Ipp64fc *> (static_cast<void *> (num.data()));
    Ipp64fc *pDst  = static_cast<Ipp64fc *> (static_cast<void *> (res.data()));
    IppStatus status = ippsDiv_64fc(pSrc1, pSrc2, pDst, n);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Division failed");
        return -1;
    }
    return 0;
}
/*!
 * @brief Multiplies a vector by a scalar: x = val*x inplace.
 * @param[in] val   Scalar to multiply x by.
 * @param[inout] x  On input this is a vector.  On exit this the scaled vector
 *                  so that x = val*x.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_ipps
 */
int RTSeis::Utils::IPPS::MulC(
    const std::complex<double> val, std::vector<std::complex<double>> x)
{
    int n = static_cast<int> (x.size());
    if (n == 0){return 0;}
    Ipp64fc *pSrcDst = static_cast<Ipp64fc *> (static_cast<void *> (x.data()));
    Ipp64fc zval;
    zval.re = std::real(val);
    zval.im = std::imag(val);
    IppStatus status = ippsMulC_64fc_I(zval, pSrcDst, n);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Division failed");
        return -1;
    }
    return 0;
}
