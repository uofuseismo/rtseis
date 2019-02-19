#include <cstdio>
#include <cstdlib>
#include <complex>
#include <vector>
#include <functional>
#include <algorithm>
#ifdef __INTEL_COMPILER
#include <pstl/algorithm>
#include <pstl/execution>
#endif
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/utilities/vectorMath.hpp"
#include <ipps.h>

using namespace RTSeis::Utilities::Math::Vector;

int RTSeis::Utilities::Math::Vector::divide(
    const std::vector<std::complex<double>> &den,
    const std::vector<std::complex<double>> &num,
    std::vector<std::complex<double>> &res)
{
    int nDen = static_cast<int> (den.size());
    int nNum = static_cast<int> (num.size());
    if (nNum != nDen)
    {
        RTSEIS_ERRMSG("num has length = %d while den has length = %d",
                      nNum, nDen);
        return -1;
    }
    res.resize(nNum);
    if (nNum <= 0){return 0;}
#ifdef __INTEL_COMPILER
    std::transform(pstl::execution::unseq,
                   num.data(), num.data()+nNum, den.data(), res.data(),
                   std::divides< std::complex<double> > ());
#else
    const Ipp64fc *pSrc1 = static_cast<const Ipp64fc *>
                           (static_cast<const void *> (den.data()));
    const Ipp64fc *pSrc2 = static_cast<const Ipp64fc *>
                           (static_cast<const void *> (num.data()));
    Ipp64fc *pDst  = static_cast<Ipp64fc *> (static_cast<void *> (res.data()));
    IppStatus status = ippsDiv_64fc(pSrc1, pSrc2, pDst, nNum); 
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Division failed");
        return -1;
    }
#endif
    return 0;
}

int RTSeis::Utilities::Math::Vector::real(
    const std::vector<std::complex<double>> &z,
    std::vector<double> &r)
{
    int n = static_cast<int> (z.size());
    r.resize(n);
    if (n <= 0){return 0;}
    const Ipp64fc *pSrc = static_cast<const Ipp64fc *>
                          (static_cast<const void *> (z.data()));
    Ipp64f *pDst  = r.data();
    IppStatus status = ippsReal_64fc(pSrc, pDst, n); 
    if (status != ippStsNoErr)
    {   
        RTSEIS_ERRMSG("%s", "Division failed");
        return -1; 
    }   
    return 0;
}
