#include <iostream>
#include <complex>
#include <vector>
#include <functional>
#include <algorithm>
#ifdef __INTEL_COMPILER
#include <pstl/algorithm>
#include <pstl/execution>
#endif
#include <ipps.h>
#include "rtseis/utilities/math/vectorMath.hpp"

int RTSeis::Utilities::Math::VectorMath::divide(
    const std::vector<std::complex<double>> &den,
    const std::vector<std::complex<double>> &num,
    std::vector<std::complex<double>> &res)
{
    int nDen = static_cast<int> (den.size());
    int nNum = static_cast<int> (num.size());
    if (nNum != nDen)
    {
        std::cerr << "Error: num has length = "
                  << nNum << " while den has length = "
                  << nDen << std::endl;
        return -1;
    }
    res.resize(nNum);
    if (nNum <= 0){return 0;}
#ifdef __INTEL_COMPILER
    std::transform(pstl::execution::unseq,
                   num.data(), num.data()+nNum, den.data(), res.data(),
                   std::divides< std::complex<double> > ());
#else
    auto *pSrc1 = static_cast<const Ipp64fc *>
                  (static_cast<const void *> (den.data()));
    auto *pSrc2 = static_cast<const Ipp64fc *>
                  (static_cast<const void *> (num.data()));
    auto *pDst  = static_cast<Ipp64fc *> (static_cast<void *> (res.data()));
    IppStatus status = ippsDiv_64fc(pSrc1, pSrc2, pDst, nNum); 
    if (status != ippStsNoErr)
    {
        std::cerr << "Division failed" << std::endl;
        return -1;
    }
#endif
    return 0;
}

//============================================================================//

int RTSeis::Utilities::Math::VectorMath::real(
    const std::vector<std::complex<double>> &z,
    std::vector<double> &r)
{
    int n = static_cast<int> (z.size());
    r.resize(n);
    if (n <= 0){return 0;}
    auto *pSrc = static_cast<const Ipp64fc *>
                 (static_cast<const void *> (z.data()));
    Ipp64f *pDst  = r.data();
    IppStatus status = ippsReal_64fc(pSrc, pDst, n); 
    if (status != ippStsNoErr)
    {   
        std::cerr << "Division failed" << std::endl;
        return -1; 
    }   
    return 0;
}

//============================================================================//
template<typename T> int RTSeis::Utilities::Math::VectorMath::copysign(
    const  std::vector<T> &x, std::vector<T> &y)
{
    int nx = static_cast<int> (x.size());
    y.resize(nx);
    if (nx <= 0){return 0;}
    int ierr = copysign(nx, x.data(), y.data());
    if (ierr != 0){y.resize(0);}
    return ierr;
}

template<typename T> int RTSeis::Utilities::Math::VectorMath::copysign(
    const int n, const T x[], T y[])
{
    if (n <= 0){return 0;}
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){std::cerr << "x is NULL" << std::endl;}
        if (y == nullptr){std::cerr << "y is NULL" << std::endl;}
        return -1;
    }
    constexpr T one = 1;
    #pragma omp simd
    for (int i=0; i<n; i++){y[i] = std::copysign(one, x[i]);}
    return 0;
}
// Instantiate the templates in the library
template int RTSeis::Utilities::Math::VectorMath::copysign<double>(
    const std::vector<double> &x, std::vector<double> &y);
template int RTSeis::Utilities::Math::VectorMath::copysign<double>(
    const int n, const double x[], double y[]);
template int RTSeis::Utilities::Math::VectorMath::copysign<float>(
    const std::vector<float> &x, std::vector<float> &y);
template int RTSeis::Utilities::Math::VectorMath::copysign<float>(
    const int n, const float x[], float y[]);

//============================================================================//

template<typename T> bool RTSeis::Utilities::Math::VectorMath::isSorted(
    const int npts, const T x[])
{
    bool lsorted;
#ifdef __INTEL_COMPILER
    //lsorted = std::is_sorted(pstl::execution::unseq, x, x+npts,
    //                           std::less_equal<T>());
    lsorted = std::is_sorted(x, x+npts);
#else
    lsorted = std::is_sorted(x, x+npts);
#endif
    return lsorted;
}
template<typename T> bool RTSeis::Utilities::Math::VectorMath::isSorted(
    const std::vector<T> &x)
{
    bool lsorted;
#ifdef __INTEL_COMPILER
    //lsorted = std::is_sorted(pstl::execution::unseq, x.begin(), x.end(),
    //                           std::less_equal<T>());
    lsorted = std::is_sorted(x.begin(), x.end());//, std::less_equal<T>());
#else
    lsorted = std::is_sorted(x.begin(), x.end());
#endif
    return lsorted;
}
template bool RTSeis::Utilities::Math::VectorMath::isSorted<double>(
    const int npts, const double x[]);
template bool RTSeis::Utilities::Math::VectorMath::isSorted<float>(
    const int npts, const float x[]);
template bool RTSeis::Utilities::Math::VectorMath::isSorted<int>(
    const int npts, const int x[]);

template bool RTSeis::Utilities::Math::VectorMath::isSorted<double>(
    const std::vector<double> &x);
template bool RTSeis::Utilities::Math::VectorMath::isSorted<float>(
    const std::vector<float> &x);
template bool RTSeis::Utilities::Math::VectorMath::isSorted<int>(
    const std::vector<int> &x);
