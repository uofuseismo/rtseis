#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cassert>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "private/throw.hpp"
#include "rtseis/utilities/math/convolve.hpp"
#include "rtseis/log.h"

using namespace::RTSeis::Utilities::Math;


static IppEnum getImplementation(const Convolve::Implementation implementation);
static std::pair<int,int> computeTrimIndices(const Convolve::Mode mode, const int n1, const int n2);

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
        if (src1Len < 1){RTSEIS_THROW_IA("%s", "No points in a");}
        RTSEIS_THROW_IA("%s", "No points in b");
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
        if (src1Len < 1){RTSEIS_THROW_IA("%s", "No points in a");}
        if (src2Len < 1){RTSEIS_THROW_IA("%s", "No points in b");}
        if (a == nullptr){RTSEIS_THROW_IA("%s", "a is NULL");}
        RTSEIS_THROW_IA("%s", "b is NULL");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int fullLen = indexes.second - indexes.first;
    double *c = *cIn;
    if (maxc < fullLen || c == nullptr)
    {
       if (maxc < fullLen)
       {
           RTSEIS_THROW_IA("maxc = %d must be at least %d", maxc, fullLen);
       }
       if (c == nullptr){RTSEIS_THROW_IA("%s", "c is NULL");}
       RTSEIS_THROW_IA("%s", "Invalid arguments");
    }
    int len = src1Len + src2Len - 1;
    // Figure out the buffer size
    int bufSize = 0;
    IppEnum funCfg = getImplementation(implementation);
    IppStatus status = ippsConvolveGetBufferSize(src1Len, src2Len, ipp64f,
                                                 funCfg, &bufSize);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute buffer size");
        return;
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
        RTSEIS_ERRMSG("%s", "Failed to compute convolution");
        if (mode != Mode::FULL){ippsFree(pDst);}
        return;
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
        if (src1Len < 1){RTSEIS_THROW_IA("%s", "No points in a");}
        RTSEIS_THROW_IA("%s", "No points in b");
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
        if (src1Len < 1){RTSEIS_THROW_IA("%s", "No points in a");}
        if (src2Len < 1){RTSEIS_THROW_IA("%s", "No points in b");}
        if (a == nullptr){RTSEIS_THROW_IA("%s", "a is NULL");}
        if (b == nullptr){RTSEIS_THROW_IA("%s", "b is NULL");}
        RTSEIS_THROW_IA("%s", "Invalid lengths");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int fullLen = indexes.second - indexes.first;
    double *c = *cIn;
    if (maxc < fullLen || c == nullptr)
    {
       if (maxc < fullLen)
       {
           RTSEIS_THROW_IA("maxc = %d must be at least %d", maxc, fullLen);
       }
       if (c == nullptr){RTSEIS_THROW_IA("%s", "c is NULL");}
       RTSEIS_THROW_IA("%s", "Invalid arguments");
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
        RTSEIS_ERRMSG("%s", "Failed to compute buffer size");
        return; 
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
        RTSEIS_ERRMSG("%s", "Failed to compute correlation");
        if (mode != Mode::FULL){ippsFree(pDst);}
        return;
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
    if (src1Len < 1){RTSEIS_THROW_IA("%s", "No points in a");}
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
        if (src1Len < 1){RTSEIS_THROW_IA("%s", "No points in a");}
        if (a == nullptr){RTSEIS_THROW_IA("%s", "a is NULL");}
        RTSEIS_THROW_IA("%s", "Invalid lengths");
    }
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src1Len);
    int fullLen = indexes.second - indexes.first;
    double *c = *cIn;
    if (maxc < fullLen || c == nullptr)
    {
       if (maxc < fullLen)
       {
           RTSEIS_THROW_IA("maxc = %d must be at least %d", maxc, fullLen);
       }
       RTSEIS_THROW_IA("%s", "c is NULL");
    }
    // Figure out the buffer size
    int len = (2*src1Len - 1)/2 + 1;
    IppEnum funCfg = getImplementation(implementation) | ippsNormNone;
    int bufSize = 0;
    IppStatus status = ippsAutoCorrNormGetBufferSize(src1Len, len,
                                                     ipp64f, funCfg, &bufSize);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute buffer size");
        return;
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
        RTSEIS_ERRMSG("%s", "Failed to compute autocorrelation");
        ippsFree(pDst);
        return;
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
