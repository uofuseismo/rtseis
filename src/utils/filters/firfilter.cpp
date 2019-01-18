#include <stdio.h>
#include <stdlib.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#define IPPS_CORE_SRC 1
#include "rtseis/log.h"
#include "rtseis/utils/filters.hpp"

using namespace RTSeis::Utils::Filters;

/*!
 * @brief Default constructor.
 * @ingroup rtseis_utils_filters_fir
 */
FIRFilter::FIRFilter(void)
{
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_utils_filters_fir
 */
FIRFilter::~FIRFilter(void)
{
    clear();
    return;
}
/*!
 * @brief Initializes the FIR filter for IPP.
 * @param[in] nb              Number of numerator coefficients.
 * @param[in] b               Numerator coefficients.  This is an array of
 *                            dimension [nb].
 * @param[in] lisRealTime     If true then the filter is for real-time
 *                            processing.
 * @param[in] lisRealTime     Otherwise, the filter is for post-processing.
 * @param[in] precision       Precision of FIR filter. 
 * @param[in] implementation  Defines the implementation.  This can specify
 *                            Implementation::DIRECT form for direct form
 *                            implementation or Implementation::FFT for an
 *                            FFT-based overlap-add implementation which 
 *                            can be advantageous when \f$ \log_2 L < N \f$
 *                            where \f$ L \f$ is the signal length and \f N \f$
 *                            the number of filter taps.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_fir
 */
int FIRFilter::initialize(const int nb, const double b[],
                              const bool lisRealTime,
                              const enum rtseisPrecision_enum precision,
                              Implementation implementation)
{
    clear();
    // Checks
    if (nb < 1 || b == NULL)
    {
        if (nb < 1){RTSEIS_ERRMSG("%s", "No b coefficients\n");}
        if (b == NULL){RTSEIS_ERRMSG("%s", "b is NULL\n");}
        return -1; 
    }
    // Figure out sizes and save some basic info
    tapsLen_ = nb;
    order_ = nb - 1;
    nbDly_ = std::max(128, order_+1);
    tapsRef_ = ippsMalloc_64f(nb);
    ippsCopy_64f(b, tapsRef_, nb);
    zi_ = ippsMalloc_64f(std::max(1, order_));
    ippsZero_64f(zi_, std::max(1, order_));
    // Determine the algorithm type
    IppAlgType algType = ippAlgDirect;
    if (implementation == Implementation::FFT){algType = ippAlgFFT;}
    // Initialize FIR filter
    if (precision == RTSEIS_DOUBLE)
    {
        Ipp64f *dlysrc = ippsMalloc_64f(nbDly_);
        Ipp64f *dlydst = ippsMalloc_64f(nbDly_);
        ippsZero_64f(dlysrc, nbDly_);
        ippsZero_64f(dlydst, nbDly_);
        Ipp64f *pTaps = ippsMalloc_64f(tapsLen_);
        ippsZero_64f(pTaps, tapsLen_);
        ippsCopy_64f(b, pTaps, nb);
        IppStatus status = ippsFIRSRGetSize(tapsLen_, ipp64f,
                                            &specSize_, &bufferSize_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting state size");
            return -1;
        }
        IppsFIRSpec_64f *pFIRspec
            = reinterpret_cast<IppsFIRSpec_64f *> (ippsMalloc_8u(specSize_));
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        status = ippsFIRSRInit_64f(pTaps, tapsLen_, algType, pFIRspec); 
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error initializing state structure");
            return -1;
        }
        // Set pointers 
        dlysrc64_ = dlysrc;
        dlydst64_ = dlydst;
        pTaps64_ = pTaps;
        pFIRspec64_ = pFIRspec;
        pBuf_ = pBuf;
    }
    else if (precision == RTSEIS_FLOAT)
    {
        Ipp32f *dlysrc = ippsMalloc_32f(nbDly_);
        Ipp32f *dlydst = ippsMalloc_32f(nbDly_);
        ippsZero_32f(dlysrc, nbDly_);
        ippsZero_32f(dlydst, nbDly_);
        Ipp32f *pTaps = ippsMalloc_32f(tapsLen_);
        ippsZero_32f(pTaps, tapsLen_);
        ippsConvert_64f32f(b, pTaps, nb);
        IppStatus status = ippsFIRSRGetSize(tapsLen_, ipp32f,
                                            &specSize_, &bufferSize_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error getting state size");
            return -1; 
        }
        IppsFIRSpec_32f *pFIRspec
            = reinterpret_cast<IppsFIRSpec_32f *> (ippsMalloc_8u(specSize_));
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        status = ippsFIRSRInit_32f(pTaps, tapsLen_, algType, pFIRspec); 
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error initializing state structure");
            return -1; 
        }
        // Set pointers 
        dlysrc32_ = dlysrc;
        dlydst32_ = dlydst;
        pTaps32_ = pTaps;
        pFIRspec32_ = pFIRspec;
        pBuf_ = pBuf;
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        clear();
        return -1;
    }
    implementation_ = implementation;
    lrt_ = lisRealTime;
    linit_ = true;
    return 0;
}
/*!
 * @brief Clears the module and resets all parameters.
 * @ingroup rtseis_utils_filters_fir
 */
void FIRFilter::clear(void)
{
    if (isInitialized())
    {
        if (precision_ == RTSEIS_DOUBLE)
        {
             IppsFIRSpec_64f *pSpec = static_cast<IppsFIRSpec_64f *> (pFIRspec64_);//.get();
             ippsFree(pSpec);
             Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc64_);
             ippsFree(dlysrc);
             Ipp64f *pTaps  = static_cast<Ipp64f *> (pTaps64_);
             ippsFree(pTaps);
        }
        else
        {
             IppsFIRSpec_32f *pSpec = static_cast<IppsFIRSpec_32f *> (pFIRspec32_);
             ippsFree(pSpec);
             Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc32_);
             ippsFree(dlysrc);
             Ipp32f *pTaps  = static_cast<Ipp32f *> (pTaps32_);
             ippsFree(pTaps);
        }
        Ipp8u *pBuf = static_cast<Ipp8u *> (pBuf_);
        ippsFree(pBuf);
        ippsFree(tapsRef_);
        ippsFree(zi_);
    }
    pFIRspec64_ = nullptr;
    pTaps64_ = nullptr;
    dlysrc64_ = nullptr;
    pFIRspec64_ = nullptr;
    pTaps32_ = nullptr;
    dlysrc64_ = nullptr;
    tapsLen_ = 0;
    nbDly_ = 0;
    bufferSize_ = 0;
    specSize_ = 0;
    order_ = 0;
    implementation_ = Implementation::DIRECT;
    precision_ = RTSEIS_DOUBLE;
    lrt_ = false;
    linit_ = false;
    return;
}

int FIRFilter::apply(const int n,
                         const double x[],
                         double y[])
{
    if (n <= 0){return 0;} // Nothing to do
    if (x == NULL || y == NULL)
    {
        if (x == NULL){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == NULL){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    if (precision_ == RTSEIS_DOUBLE)
    {
        if (isRealTime())
        {
            Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc64_);
            Ipp64f *dlydst = static_cast<Ipp64f *> (dlydst64_);

            ippsCopy_64f(dlydst, dlysrc, order_);
        }
        else
        {

        }
    }
    else
    {

    }
    return 0;
}
