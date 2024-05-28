#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>
#include <mkl.h>
#include <mkl_vsl.h>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/utilities/math/convolve.hpp"
#include "private/convolve.hpp"

using namespace  RTSeis::Utilities::Math;


//static IppEnum getImplementation(const Convolve::Implementation implementation);
//static std::pair<int,int> computeTrimIndices(const Convolve::Mode mode, const int n1, const int n2);

namespace
{
/// @brief Utility function to get appropriate IPP implementaiton.
/// @param[in] implementatation  The desired implementation.
/// @result The corresponding IPP implementation enum.
/// @ingroup rtseis_utils_convolve
IppEnum getImplementation(const Convolve::Implementation implementation)
{
    if (implementation == Convolve::Implementation::FFT)
    {
        return ippAlgFFT;
    }
    else if (implementation == Convolve::Implementation::DIRECT)
    {
        return ippAlgDirect;
    }
    else
    {
        return ippAlgAuto;
    }
}
MKL_INT getMKLImplementationConvolution(
    const Convolve::Implementation implementation)
{
    if (implementation == Convolve::Implementation::FFT)
    {
        return VSL_CONV_MODE_FFT;
    }   
    else if (implementation == Convolve::Implementation::DIRECT)
    {   
        return VSL_CONV_MODE_DIRECT;
    }   
    else
    {   
        return VSL_CONV_MODE_AUTO;
    }   
}
}

std::vector<double>
Convolve::convolve(const std::vector<double> &a,
                   const std::vector<double> &b, 
                   const Convolve::Mode mode,
                   const Convolve::Implementation implementation)
{
    std::vector<double> c;
    int src1Len = static_cast<size_t> (a.size());
    int src2Len = static_cast<size_t> (b.size());
    c.resize(0);
    if (src1Len < 1 || src2Len < 1)
    {
        if (src1Len < 1){throw std::invalid_argument("No points in a");}
        throw std::invalid_argument("No points in b");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int len = indexes.second - indexes.first;
    int nc;
    c.resize(len);
    double *cdata = c.data();
    convolve(src1Len, a.data(), src2Len, b.data(),
             len, &nc, &cdata, mode, implementation);
#ifdef DEBUG
    assert(nc == static_cast<int> (c.size()));
#endif
    return c;
}

std::vector<std::complex<double>>
Convolve::convolve(const std::vector<std::complex<double>> &a, 
                   const std::vector<std::complex<double>> &b, 
                   const Convolve::Mode mode,
                   const Convolve::Implementation implementation)
{
    std::vector<std::complex<double>> c;
    int src1Len = static_cast<size_t> (a.size());
    int src2Len = static_cast<size_t> (b.size());
    c.resize(0);
    if (src1Len < 1 || src2Len < 1)
    {
        if (src1Len < 1){throw std::invalid_argument("No points in a");}
        throw std::invalid_argument("No points in b");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int len = indexes.second - indexes.first;
    int nc;
    c.resize(len, 0);
    std::complex<double> *cdata = c.data();
    convolve(src1Len, a.data(), src2Len, b.data(),
             len, &nc, &cdata, mode, implementation);
#ifdef DEBUG
    assert(nc == static_cast<int> (c.size()));
#endif
    return c;
}

void Convolve::convolve(const int src1Len, const double a[],
                        const int src2Len, const double b[],
                        const int maxc, int *nc, double *cIn[],
                        const Convolve::Mode mode,
                        const Convolve::Implementation implementation)
{
    // Check the inputs
    *nc = 0;
    if (src1Len < 1 || a == nullptr || src2Len < 1 || b == nullptr)
    {
        if (src1Len < 1){throw std::invalid_argument("No points in a");}
        if (src2Len < 1){throw std::invalid_argument("No points in b");}
        if (a == nullptr){throw std::invalid_argument("a is NULL");}
        throw std::invalid_argument("b is NULL");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int fullLen = indexes.second - indexes.first;
    double *c = *cIn;
    if (maxc < fullLen || c == nullptr)
    {
       if (maxc < fullLen)
       {
           throw std::invalid_argument("maxc = " + std::to_string(maxc)
                                     + " must be at "
                                     + std::to_string(fullLen));
       }
       throw std::invalid_argument("c is NULL");
    }
    int len = src1Len + src2Len - 1;
    // Figure out the buffer size
    int bufSize = 0;
    IppEnum funCfg = getImplementation(implementation);
    IppStatus status = ippsConvolveGetBufferSize(src1Len, src2Len, ipp64f,
                                                 funCfg, &bufSize);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("Failed to compute buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    // Perform the convolution
    const double *pSrc1 = a;
    const double *pSrc2 = b;
    double *pDst = nullptr;
    if (mode == Mode::FULL)
    {
        pDst = c;
    }
    else
    {
        pDst = ippsMalloc_64f(len);
    }
    status = ippsConvolve_64f(pSrc1, src1Len, pSrc2, src2Len, pDst,
                              funCfg, pBuffer);
    ippsFree(pBuffer);
    if (status != ippStsNoErr)
    {
        if (mode != Mode::FULL){ippsFree(pDst);}
        throw std::runtime_error("Failed to compute convolution");
    }
    // The full convolution is desired
    *nc = fullLen;
    if (mode == Mode::FULL){return;}
    // Trim the full convolution
    int i1 = indexes.first;
    int i2 = indexes.second;
    len = i2 - i1;
    ippsCopy_64f(&pDst[i1], c, len);
    ippsFree(pDst);
}

void Convolve::convolve(const int src1Len, const std::complex<double> a[],
                        const int src2Len, const std::complex<double> b[],
                        const int maxc, int *nc, std::complex<double> *cIn[],
                        const Convolve::Mode mode,
                        const Convolve::Implementation implementation)
{
    // Check the inputs
    *nc = 0;
    if (src1Len < 1 || a == nullptr || src2Len < 1 || b == nullptr)
    {
        if (src1Len < 1){throw std::invalid_argument("No points in a");}
        if (src2Len < 1){throw std::invalid_argument("No points in b");}
        if (a == nullptr){throw std::invalid_argument("a is NULL");}
        throw std::invalid_argument("b is NULL");
    }
    // Get indices
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int outputLen = indexes.second - indexes.first;
    // Check output size is okay
    if (maxc < outputLen || *cIn == nullptr)
    {
       if (maxc < outputLen)
       {
           throw std::invalid_argument("maxc = " + std::to_string(maxc)
                                     + " must be at "
                                     + std::to_string(outputLen));
       }
       throw std::invalid_argument("c is NULL");
    }
    int fullLen = src1Len + src2Len - 1;
    // Setup for VSL
    VSLConvTaskPtr task;
    auto mklImplementation = getMKLImplementationConvolution(implementation);
    auto status = vslzConvNewTask1D(&task, mklImplementation,
                                    src1Len, src2Len, fullLen);
    if (status != VSL_STATUS_OK)
    {
        throw std::runtime_error("Failed to create convolution task");
    }
    auto aMKL = reinterpret_cast<const MKL_Complex16 *> (a);
    auto bMKL = reinterpret_cast<const MKL_Complex16 *> (b);
    MKL_Complex16 *cMKL{nullptr};
    if (mode == Mode::FULL)
    {
        cMKL = reinterpret_cast<MKL_Complex16 *> (*cIn);
    }
    else
    {
        cMKL = reinterpret_cast<MKL_Complex16 *>
               (mkl_calloc(fullLen, sizeof(MKL_Complex16), 64));
    }
    status = vslzConvExec1D(task, aMKL, 1, bMKL, 1, cMKL, 1);
    if (status != 0)
    {
        if (mode != Mode::FULL){mkl_free(cMKL);}
        throw std::runtime_error("Failed to compute convolution");
    }
    auto error = vslConvDeleteTask(&task);
    if (error != 0)
    {
        if (mode != Mode::FULL){mkl_free(cMKL);}
        throw std::runtime_error("Failed to delete task");
    }
    // The full convolution is desired
    *nc = outputLen;
    if (mode == Mode::FULL){return;}
    // Trim the full convolution
    int i1 = indexes.first;
    int i2 = indexes.second;
    auto pDst = reinterpret_cast<const std::complex<double> *> (cMKL);
    auto cOut = *cIn;
    std::copy(pDst + i1, pDst + i2, cOut + 0);
    mkl_free(cMKL);
}

//============================================================================//
//                                     Correlate                              //
//============================================================================//

std::vector<double>
Convolve::correlate(const std::vector<double> &a,
                    const std::vector<double> &b, 
                    const Convolve::Mode mode,
                    const Convolve::Implementation implementation)
{
    std::vector<double> c;
    int src1Len = static_cast<size_t> (a.size());
    int src2Len = static_cast<size_t> (b.size());
    c.resize(0);
    if (src1Len < 1 || src2Len < 1)
    {
        if (src1Len < 1){throw std::invalid_argument("No points in a");}
        throw std::invalid_argument("no points in b");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int len = indexes.second - indexes.first;
    int nc;
    c.resize(len);
    double *cdata = c.data();
    correlate(src1Len, a.data(), src2Len, b.data(),
              len, &nc, &cdata, mode, implementation);
#ifdef DEBUG
    assert(nc == static_cast<int> (c.size()));
#endif
    return c;
}

void Convolve::correlate(const int src1Len, const double a[],
                         const int src2Len, const double b[],
                         const int maxc, int *nc, double *cIn[],
                         const Convolve::Mode mode,
                         const Convolve::Implementation implementation)
{
    // Check the inputs
    *nc = 0;
    if (src1Len < 1 || a == nullptr || src2Len < 1 || b == nullptr)
    {
        if (src1Len < 1){throw std::invalid_argument("No points in a");}
        if (src2Len < 1){throw std::invalid_argument("No points in b");}
        if (a == nullptr){throw std::invalid_argument("a is NULL");}
        throw std::invalid_argument("b is NULL");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int fullLen = indexes.second - indexes.first;
    double *c = *cIn;
    if (maxc < fullLen || c == nullptr)
    {
       if (maxc < fullLen)
       {
           throw std::invalid_argument("maxc = " + std::to_string(maxc)
                                     + " must be at "
                                     + std::to_string(fullLen));
       }
       throw std::invalid_argument("c is NULL");
    }
    int len = src1Len + src2Len - 1;
    // Figure out the buffer size
    IppEnum funCfg = getImplementation(implementation) | ippsNormNone;
    int bufSize = 0;
    const int lowLag =-src2Len + 1; //(std::max(src1Len, src2Len) - 1);
    IppStatus status = ippsCrossCorrNormGetBufferSize(src2Len, src1Len, len,
                                                      lowLag, ipp64f, funCfg,
                                                      &bufSize);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("Failed to compute buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    // Perform the correlation
    const double *pSrc1 = a;
    const double *pSrc2 = b;
    double *pDst = nullptr;
    if (mode == Mode::FULL)
    {
        pDst = c;
    }
    else
    {
        pDst = ippsMalloc_64f(len);
    }
    // Perform the correlation noting that IPP uses a formula whose convention
    // is reverse from Matlab's
    status = ippsCrossCorrNorm_64f(pSrc2, src2Len, pSrc1, src1Len,
                                   pDst, len, lowLag, funCfg,
                                   pBuffer);
    ippsFree(pBuffer);
    if (status != ippStsNoErr)
    {
        if (mode != Mode::FULL){ippsFree(pDst);}
        throw std::runtime_error("Failed to compute correlation");
    }
    // The full correlation is desired
    *nc = fullLen;
    if (mode == Mode::FULL){return;}
    // Trim the full correlation
    int i1 = indexes.first;
    int i2 = indexes.second;
    len = i2 - i1;
    ippsCopy_64f(&pDst[i1], c, len);
    ippsFree(pDst);
}

//============================================================================//

std::vector<double>
Convolve::autocorrelate(const std::vector<double> &a, 
                        const Convolve::Mode mode,
                        const Convolve::Implementation implementation)
{
    std::vector<double> c;
    int src1Len = static_cast<size_t> (a.size());
    c.resize(0);
    if (src1Len < 1){throw std::invalid_argument("No points in a");}
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src1Len);
    int len = indexes.second - indexes.first;
    int nc;
    c.resize(len);
    double *cdata = c.data();
    autocorrelate(src1Len, a.data(),
                  len, &nc, &cdata, mode, implementation);
#ifdef DEBUG
    assert(nc == static_cast<int> (c.size()));
#endif
    return c;
}

void Convolve::autocorrelate(const int src1Len, const double a[],
                             const int maxc, int *nc, double *cIn[],
                             const Convolve::Mode mode,
                             const Convolve::Implementation implementation)
{
    // Check the inputs
    *nc = 0;
    if (src1Len < 1 || a == nullptr)
    {
        if (src1Len < 1){throw std::invalid_argument("No points in a");}
        throw std::invalid_argument("a is NULL");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src1Len);
    int fullLen = indexes.second - indexes.first;
    double *c = *cIn;
    if (maxc < fullLen || c == nullptr)
    {
       if (maxc < fullLen)
       {
           throw std::invalid_argument("maxc = " + std::to_string(maxc)
                                     + " must be at "
                                     + std::to_string(fullLen));
       }
       throw std::invalid_argument("c is NULL");
    }
    // Figure out the buffer size
    int len = (2*src1Len - 1)/2 + 1;
    IppEnum funCfg = getImplementation(implementation) | ippsNormNone;
    int bufSize = 0;
    IppStatus status = ippsAutoCorrNormGetBufferSize(src1Len, len,
                                                     ipp64f, funCfg, &bufSize);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("Failed to compute buffer size");
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    // Perform the autocorrelation
    const double *pSrc1 = a;
    double *pDst = nullptr;
    pDst = ippsMalloc_64f(len);
    status = ippsAutoCorrNorm_64f(pSrc1, src1Len, pDst, len, funCfg, pBuffer);
    ippsFree(pBuffer);
    if (status != ippStsNoErr)
    {
        ippsFree(pDst);
        throw std::runtime_error("Failed to compute autocorrelation");
    }
    // Only the first half is computed
    if (mode == Mode::FULL)
    {
        ippsFlip_64f(pDst, c, len);
        ippsCopy_64f(&pDst[1], &c[len], len-1); 
        ippsFree(pDst);
    }
    else
    {
        Ipp64f *scratch = ippsMalloc_64f(2*src1Len - 1);
        ippsFlip_64f(pDst, scratch, len);
        ippsCopy_64f(&pDst[1], &scratch[len], len-1); 
        ippsFree(pDst);
        // Trim the full autocorrelation
        int i1 = indexes.first;
        int i2 = indexes.second;
        len = i2 - i1; 
        ippsCopy_64f(&scratch[i1], c, len);
        ippsFree(scratch);
    }
    *nc = fullLen;
}

//============================================================================//

int Convolve::computeConvolutionLength(const int n1, const int n2,
                                       const Convolve::Mode mode)
{
    std::pair<int, int> indexes = computeTrimIndices(mode, n1, n2);
    int len = indexes.second - indexes.first; 
    return len;
}

/*!
 * @brief Utility function to get appropriate IPP implementaiton.
 * @param[in] implementatation  The desired implementation.
 * @result The corresponding IPP implementation enum.
 * @ingroup rtseis_utils_convolve
 */
/*
IppEnum getImplementation(const Convolve::Implementation implementation)
{
    if (implementation == Convolve::Implementation::FFT)
    {   
        return ippAlgFFT;
    }   
    else if (implementation == Convolve::Implementation::DIRECT)
    {   
        return ippAlgDirect;
    }   
    else
    {   
        return ippAlgAuto;
    }   
}
*/

/*!
 * @brief Computes the start and stop indices of the convolution or correlation.
 *        The resulting copy would go from conv(res.first:res.second) where the
 *        upper limit, res.second, is non-inclusive.
 * @param[in] mode   Defines the convolution or correlation mode.
 * @param[in] n1     Length of the first array.
 * @param[in] n2     Length of the second array.
 * @result The start and stop index from which to copy the full convolution.
 * @ingroup rtseis_utils_convolve
 */
/*
std::pair<int,int> computeTrimIndices(const Convolve::Mode mode,
                                      const int n1,  const int n2)
{
    int lc;
    int nLeft = 0;
    int nRight = 0;
    // Full
    if (mode == Convolve::Mode::FULL)
    {
        lc = n1 + n2 - 1; // Length of full convolution
        nLeft = 0;
        nRight = lc; 
    }
    // Valid
    else if (mode == Convolve::Mode::VALID)
    {
        lc = std::max(n1,n2) - std::min(n1,n2) + 1;
        if (n2 > n1)
        {
            nLeft = n1/2 + 1;
            nRight = nLeft + lc; 
        }
        else
        {
            nLeft = n2/2 + 1;
            nRight = nLeft + lc; 
        }
    }
    // Same
    else if (mode == Convolve::Mode::SAME)
    {
        lc = std::max(n1, n2);
        if (n1 < n2)
        {
            nLeft = n1/2;
            nRight = nLeft + lc;
        }
        else
        {
            nLeft = n2/2;
            nRight = nLeft + lc;
        }
    }
    std::pair<int,int> result(nLeft, nRight);
    return result;
}
*/
