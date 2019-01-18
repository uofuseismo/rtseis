#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include <cmath>
#include <vector>
#include "rtseis/utils/design.hpp"
#include "rtseis/utils/convolve.hpp"
#include "rtseis/log.h"
#include <ipps.h>

using namespace::RTSeis::Utils::Math;

/*!
 * @defgroup rtseis_utils_convolve Convolution and Correlation
 * @brief Utility functions for convolution and correlation.
 *        This code is originally from ISTI's ISCL and has been
 *        modified to conform with C++.  Function names have also been
 *        changed to conform with rtseis's naming convention.
 * @copyright ISTI distributed under the Apache 2 license.
 * @ingroup rtseis_utils
 */

static IppEnum getImplementation(const Convolve::Implementation implementation);
static std::pair<int,int> computeTrimIndices(
    const Convolve::Mode mode,
    const int n1, const int n2);

/*!
 * @brief Computes the convolution \f$ c[k] = \sum_n a[n] b[n-k] \f$.
 * @param[in] a               First array in convolution.  This has length [m].
 * @param[in] b               Second array in convolution.  This has length [n].
 * @param[out] c              The resulting convolution.
 * @param[in] mode            Defines the convolution mode which can be
 *                            FULL, VALID, or SAME.  The default is FULL.
 * @param[in] implementation  Defines the implementation type.  This can
 *                            be AUTO, DIRECT, or FFT.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_convolve
 */
int Convolve::convolve(const std::vector<double> a,
                       const std::vector<double> b,
                       std::vector<double> &c,
                       const Convolve::Mode mode,
                       const Convolve::Implementation implementation)
{
    IppEnum funCfg = getImplementation(implementation);
    int src1Len = static_cast<size_t> (a.size());
    int src2Len = static_cast<size_t> (b.size());
    c.resize(0);
    if (src1Len < 1 || src2Len < 1)
    {
        if (src1Len < 1){RTSEIS_ERRMSG("%s", "No points in a");}
        if (src2Len < 1){RTSEIS_ERRMSG("%s", "No points in b");}
        return -1;
    }
    int len = src1Len + src2Len - 1;
    // Figure out the buffer size
    int bufSize = 0;
    IppStatus status = ippsConvolveGetBufferSize(src1Len, src2Len, ipp64f,
                                                 funCfg, &bufSize);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute buffer size");
        return -1;
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    // Perform the convolution
    const double *pSrc1 = a.data();
    const double *pSrc2 = b.data();
    double *pDst = nullptr;
    if (mode == Mode::FULL)
    {
        c.resize(len);
        pDst = c.data();
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
        c.resize(0);
        return -1;
    }
    // The full convolution is desired
    if (mode == Mode::FULL){return 0;}
    // Trim the full convolution
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int i1 = indexes.first;
    int i2 = indexes.second;
    len = i2 - i1;
    c.resize(len);
#ifdef _INTEL_COMPILER
    #pragma ivdep
#endif
    std::copy(&pDst[i1], &pDst[i1+i2], c.begin());
    ippsFree(pDst);
    return 0;
}
/*!
 * @brief Computes the correlation \f$ c[k] = \sum_n a[n] b[n+k] \f$.
 * @param[in] a               First array in correlation.  This has length [m].
 * @param[in] b               Second array in correlation.  This has length [n].
 * @param[out] c              The resulting correlation.
 * @param[in] mode            Defines the correlation mode which can be
 *                            FULL, VALID, or SAME.  The default is FULL.
 * @param[in] implementation  Defines the implementation type.  This can
 *                            be AUTO, DIRECT, or FFT.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_convolve
 */
int Convolve::correlate(const std::vector<double> a,
                        const std::vector<double> b,
                        std::vector<double> &c, 
                        const Convolve::Mode mode,
                        const Convolve::Implementation implementation)
{
    IppEnum funCfg = getImplementation(implementation) | ippsNormNone;
    int src1Len = static_cast<size_t> (a.size());
    int src2Len = static_cast<size_t> (b.size());
    c.resize(0);
    if (src1Len < 1 || src2Len < 1)
    {   
        if (src1Len < 1){RTSEIS_ERRMSG("%s", "No points in a");}
        if (src2Len < 1){RTSEIS_ERRMSG("%s", "No points in b");}
        return -1; 
    }   
    int len = src1Len + src2Len - 1;
    // Figure out the buffer size
    int bufSize = 0;
    const int lowLag =-src2Len + 1; //(std::max(src1Len, src2Len) - 1);
    IppStatus status = ippsCrossCorrNormGetBufferSize(src1Len, src2Len, len,
                                                      lowLag, ipp64f, funCfg,
                                                      &bufSize);
    if (status != ippStsNoErr)
    {   
        RTSEIS_ERRMSG("%s", "Failed to compute buffer size");
        return -1; 
    }   
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    // Perform the correlation
    const double *pSrc1 = a.data();
    const double *pSrc2 = b.data();
    double *pDst = nullptr;
    if (mode == Mode::FULL)
    {   
        c.resize(len);
        pDst = c.data();
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
        c.resize(0);
        return -1;
    }
    // The full correlation is desired
    if (mode == Mode::FULL){return 0;}
    // Trim the full correlation
    std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src2Len);
    int i1 = indexes.first;
    int i2 = indexes.second;
    len = i2 - i1;
    c.resize(len);
#ifdef _INTEL_COMPILER
    #pragma ivdep
#endif
    std::copy(&pDst[i1], &pDst[i1]+i2, c.begin());
    ippsFree(pDst);
    return 0;
}
/*!
 * @brief Computes the autocorrelation \f$ c[k] = \sum_n a[n] a[n+k] \f$.
 * @ingroup rtseis_utils_convolve
 * @param[in] a               Array to autocorrelation.  This has length [m].
 * @param[out] c              The resulting autocorrelation.
 * @param[in] mode            Defines the correlation mode which can be
 *                            FULL, VALID, or SAME.  The default is FULL.
 * @param[in] implementation  Defines the implementation type.  This can
 *                            be AUTO, DIRECT, or FFT.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_convolve
 */
int Convolve::autocorrelate(const std::vector<double> a,
                            std::vector<double> &c, 
                            const Convolve::Mode mode,
                            const Convolve::Implementation implementation)
{
    IppEnum funCfg = getImplementation(implementation) | ippsNormNone;
    int src1Len = static_cast<size_t> (a.size());
    c.resize(0);
    if (src1Len < 1)
    {
        if (src1Len < 1){RTSEIS_ERRMSG("%s", "No points in a");}
        return -1;
    }
    int len = (2*src1Len - 1)/2 + 1;
    // Figure out the buffer size
    int bufSize = 0;
    IppStatus status = ippsAutoCorrNormGetBufferSize(src1Len, len,
                                                     ipp64f, funCfg, &bufSize);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute buffer size");
        return -1;
    }
    Ipp8u *pBuffer = ippsMalloc_8u(bufSize);
    // Perform the autocorrelation
    const double *pSrc1 = a.data();
    double *pDst = nullptr;
    pDst = ippsMalloc_64f(len);
    status = ippsAutoCorrNorm_64f(pSrc1, src1Len, pDst, len, funCfg, pBuffer);
    ippsFree(pBuffer);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute autocorrelation");
        c.resize(0);
        return -1;
    }
    // Only the first half is computed
    if (mode == Mode::FULL)
    {
        c.resize(2*src1Len - 1);
        double *cPtr = c.data();
        ippsFlip_64f(pDst, cPtr, len);
        ippsCopy_64f(&pDst[1], &cPtr[len], len-1); 
        ippsFree(pDst);
        return 0;
    }
    else
    {
        Ipp64f *scratch = ippsMalloc_64f(2*src1Len - 1);
        ippsFlip_64f(pDst, scratch, len);
        ippsCopy_64f(&pDst[1], &scratch[len], len-1); 
        ippsFree(pDst);
        // Trim the full autocorrelation
        std::pair<int,int> indexes = computeTrimIndices(mode, src1Len, src1Len);
        int i1 = indexes.first;
        int i2 = indexes.second;
        len = i2 - i1; 
        c.resize(len);
#ifdef _INTEL_COMPILER
        #pragma ivdep
#endif
        std::copy(&scratch[i1], &scratch[i1]+i2, c.begin());
        ippsFree(scratch);
    }
    return 0;
}
/*!
 * @brief Utility function to get appropriate IPP implementaiton.
 * @param[in] implementatation  The desired implementation.
 * @result The corresponding IPP implementation enum.
 * @ingroup rtseis_utils_convolve
 */
static IppEnum getImplementation(const Convolve::Implementation implementation)
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
static std::pair<int,int> computeTrimIndices(
    const Convolve::Mode mode,
    const int n1, const int n2)
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
        nLeft = 0;
        nRight = 0;
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
