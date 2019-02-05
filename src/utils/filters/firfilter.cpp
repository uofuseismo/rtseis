#include <stdio.h>
#include <stdlib.h>
#include <ipps.h>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/utils/filters.hpp"

using namespace RTSeis::Utils::Filters;

/*!
 * @defgroup rtseis_utils_filters_fir FIR Filtering
 * @brief This is the core implementation for FIR filtering.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
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
 * @brief Copy constructor.
 * @param[in] fir  Class from which to initialize.
 * @ingroup rtseis_utils_filters_sos
 */
FIRFilter::FIRFilter(const FIRFilter &fir)
{
    *this = fir;
    return;
}
/*!
 * @brief Copy operator.
 * @param[in] fir   FIR filter class to copy.
 * @result A deep copy of the input class.
 * @ingroup rtseis_utils_filters_sos
 */
FIRFilter& FIRFilter::operator=(const FIRFilter &fir)
{
    if (&fir == this){return *this;}
    clear();
    if (!fir.linit_){return *this;}
    // Reinitialize the filter
    initialize(fir.tapsLen_, fir.tapsRef_, fir.isRealTime(),
               fir.getPrecision(), fir.implementation_);
    // Now copy the filter states
    if (bufferSize_ > 0)
    {
        Ipp8u *pBufIn = static_cast<Ipp8u *> (fir.pBuf_);
        Ipp8u *pBufOut = static_cast<Ipp8u *> (pBuf_);
        ippsCopy_8u(pBufIn, pBufOut, bufferSize_);
    }
    // Copy the initial conditions
    if (order_ > 0){ippsCopy_64f(fir.zi_, zi_, order_);}
    if (nwork_ > 0)
    {
        if (isDoublePrecision())
        {
            Ipp64f *dlysrcIn  = static_cast<Ipp64f *> (fir.dlysrc_);
            Ipp64f *dlysrcOut = static_cast<Ipp64f *> (dlysrc_);
            ippsCopy_64f(dlysrcIn, dlysrcOut, nwork_); 
            Ipp64f *dlydstIn  = static_cast<Ipp64f *> (fir.dlydst_);
            Ipp64f *dlydstOut = static_cast<Ipp64f *> (dlydst_);
            ippsCopy_64f(dlydstIn, dlydstOut, nwork_);
        }
        else
        {
            Ipp32f *dlysrcIn  = static_cast<Ipp32f *> (fir.dlysrc_);
            Ipp32f *dlysrcOut = static_cast<Ipp32f *> (dlysrc_);
            ippsCopy_32f(dlysrcIn, dlysrcOut, nwork_); 
            Ipp32f *dlydstIn  = static_cast<Ipp32f *> (fir.dlydst_);
            Ipp32f *dlydstOut = static_cast<Ipp32f *> (dlydst_);
            ippsCopy_32f(dlydstIn, dlydstOut, nwork_);
        }
    }
    return *this;
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
                              const RTSeis::Precision precision,
                              Implementation implementation)
{
    clear();
    // Checks
    if (nb < 1 || b == nullptr)
    {
        if (nb < 1){RTSEIS_ERRMSG("%s", "No b coefficients\n");}
        if (b == nullptr){RTSEIS_ERRMSG("%s", "b is NULL\n");}
        return -1; 
    }
    // Figure out sizes and save some basic info
    tapsLen_ = nb;
    order_ = nb - 1;
    nwork_ = std::max(128, order_+1);
    tapsRef_ = ippsMalloc_64f(nb);
    ippsCopy_64f(b, tapsRef_, nb);
    zi_ = ippsMalloc_64f(std::max(1, order_));
    ippsZero_64f(zi_, std::max(1, order_));
    // Determine the algorithm type
    IppAlgType algType = ippAlgDirect;
    if (implementation == Implementation::FFT){algType = ippAlgFFT;}
    // Initialize FIR filter
    if (precision == RTSeis::Precision::DOUBLE)
    {
        Ipp64f *dlysrc = ippsMalloc_64f(nwork_);
        Ipp64f *dlydst = ippsMalloc_64f(nwork_);
        ippsZero_64f(dlysrc, nwork_);
        ippsZero_64f(dlydst, nwork_);
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
        IppsFIRSpec_64f *pFIRSpec
            = reinterpret_cast<IppsFIRSpec_64f *> (ippsMalloc_8u(specSize_));
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        status = ippsFIRSRInit_64f(pTaps, tapsLen_, algType, pFIRSpec); 
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error initializing state structure");
            return -1;
        }
        // Set pointers 
        dlysrc_ = dlysrc;
        dlydst_ = dlydst;
        pTaps_ = pTaps;
        pFIRSpec_ = pFIRSpec;
        pBuf_ = pBuf;
    }
    else if (precision == RTSeis::Precision::FLOAT)
    {
        Ipp32f *dlysrc = ippsMalloc_32f(nwork_);
        Ipp32f *dlydst = ippsMalloc_32f(nwork_);
        ippsZero_32f(dlysrc, nwork_);
        ippsZero_32f(dlydst, nwork_);
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
        IppsFIRSpec_32f *pFIRSpec
            = reinterpret_cast<IppsFIRSpec_32f *> (ippsMalloc_8u(specSize_));
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        status = ippsFIRSRInit_32f(pTaps, tapsLen_, algType, pFIRSpec); 
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Error initializing state structure");
            return -1; 
        }
        // Set pointers 
        dlysrc_ = dlysrc;
        dlydst_ = dlydst;
        pTaps_ = pTaps;
        pFIRSpec_ = pFIRSpec;
        pBuf_ = pBuf;
    }
    else
    {
        RTSEIS_ERRMSG("%s", "Invalid precision");
        clear();
        return -1;
    }
    implementation_ = implementation;
    setPrecision(precision);
    toggleRealTime(lisRealTime);
    linit_ = true;
    return 0;
}
/*!
 * @brief Clears the module and resets all parameters.
 * @ingroup rtseis_utils_filters_fir
 */
void FIRFilter::clear(void)
{
    if (pFIRSpec_ != nullptr){ippsFree(pFIRSpec_);}
    if (pTaps_ != nullptr){ippsFree(pTaps_);}
    if (dlysrc_ != nullptr){ippsFree(dlysrc_);}
    if (dlydst_ != nullptr){ippsFree(dlydst_);}
    if (pBuf_ != nullptr){ippsFree(pBuf_);}
    if (tapsRef_ != nullptr){ippsFree(tapsRef_);}
    if (zi_ != nullptr){ippsFree(zi_);}
    pFIRSpec_ = nullptr;
    pTaps_ = nullptr;
    dlysrc_ = nullptr;
    dlydst_ = nullptr;
    pBuf_ = nullptr;
    tapsRef_ = nullptr;
    zi_ = nullptr;
    setPrecision(RTSeis::Precision::DOUBLE);
    toggleRealTime(false);
    tapsLen_ = 0;
    nwork_ = 0;
    bufferSize_ = 0;
    specSize_ = 0;
    order_ = 0;
    implementation_ = Implementation::DIRECT;
    linit_ = false;
    return;
}
/*!
 * @brief Sets the initial conditions for the filter.  This should be called
 *        prior to filter application as it will reset the filter.
 * @param[in] nz   The FIR filter initial conditions.
 *                 This should be equal to getInitialConditionLength().
 * @param[in] zi   The initial conditions.  This has dimension [nz].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_fir
 */
int FIRFilter::setInitialConditions(const int nz, const double zi[])
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    resetInitialConditions();
    int nzRef = getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef){RTSEIS_ERRMSG("nz=%d should equal %d", nz, nzRef);}
        if (zi == nullptr){RTSEIS_ERRMSG("%s", "zi is NULL");}
        return -1;
    }
    if (nzRef > 0)
    {
        ippsCopy_64f(zi, zi_, nzRef);
        if (isDoublePrecision())
        {
            Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
            ippsCopy_64f(zi, dlysrc, nzRef);
        }
        else
        {
            Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
            ippsConvert_64f32f(zi, dlysrc, nzRef);
        }
    }
    return 0;
}
/*!
 * @brief Resets the initial conditions on the source delay line to the
 *        default initial conditions or the initial conditions set 
 *        when FIRFilter::setInitialConditions() was called.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_fir
 */
int FIRFilter::resetInitialConditions(void)
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    if (order_ > 0)
    {
        if (isDoublePrecision())
        {
            Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
            ippsCopy_64f(zi_, dlysrc, order_);
        }
        else
        {
            Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
            ippsConvert_64f32f(zi_, dlysrc, order_);
        }
    }
    return 0;
}
/*!
 * @brief Applies the FIR filter to the data.
 * @param[in] n   Number of points in signals.
 * @param[in] x   Signal to filter.  This has dimension [n].
 * @param[out] y  The filtered signal.  This has dimension [n].
 * @ingroup rtseis_utils_filters_fir
 */
int FIRFilter::apply(const int n, const double x[], double y[])
{
    if (n <= 0){return 0;} // Nothing to do
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    if (isFloatPrecision())
    {
        Ipp32f *x32 = ippsMalloc_32f(n);
        Ipp32f *y32 = ippsMalloc_32f(n);
        ippsConvert_64f32f(x, x32, n); 
        int ierr = apply(n, x32, y32);
        ippsFree(x32);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to apply filter");
            ippsFree(y32);
            return -1; 
        }
        ippsConvert_32f64f(y32, y, n); 
        ippsFree(y32);
        return 0;
    }
    IppsFIRSpec_64f *pFIRSpec = static_cast<IppsFIRSpec_64f *> (pFIRSpec_);
    Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
    Ipp64f *dlydst = static_cast<Ipp64f *> (dlydst_);
    Ipp8u *pBuf = static_cast<Ipp8u *> (pBuf_);
    IppStatus status = ippsFIRSR_64f(x, y, n, pFIRSpec, dlysrc, dlydst, pBuf);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply FIR filter");
        return -1;
    }
    if (isRealTime() && order_ > 0){ippsCopy_64f(dlydst, dlysrc, order_);}
    return 0;
}
/*!
 * @brief Applies the FIR filter to the data.
 * @param[in] n   Number of points in signals.
 * @param[in] x   Signal to filter.  This has dimension [n].
 * @param[out] y  The filtered signal.  This has dimension [n].
 * @ingroup rtseis_utils_filters_fir
 */
int FIRFilter::apply(const int n, const float x[], float y[])
{
    if (n <= 0){return 0;} // Nothing to do
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    if (isDoublePrecision())
    {
        Ipp64f *x64 = ippsMalloc_64f(n);
        Ipp64f *y64 = ippsMalloc_64f(n);
        ippsConvert_32f64f(x, x64, n); 
        int ierr = apply(n, x64, y64);
        ippsFree(x64);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to apply filter");
            ippsFree(y64);
            return -1; 
        }
        ippsConvert_64f32f(y64, y, n); 
        ippsFree(y64);
        return 0;
    }
    IppsFIRSpec_32f *pFIRSpec = static_cast<IppsFIRSpec_32f *> (pFIRSpec_);
    Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
    Ipp32f *dlydst = static_cast<Ipp32f *> (dlydst_);
    Ipp8u *pBuf = static_cast<Ipp8u *> (pBuf_);
    IppStatus status = ippsFIRSR_32f(x, y, n, pFIRSpec, dlysrc, dlydst, pBuf);
    if (status != ippStsNoErr)
    {   
        RTSEIS_ERRMSG("%s", "Failed to apply FIR filter");
        return -1;
    }
    if (isRealTime() && order_ > 0){ippsCopy_32f(dlydst, dlysrc, order_);}
    return 0;
}
/*!
 * @brief Returns the length of the initial conditions.
 * @result The length of initial conditions array.
 * @ingroup rtseis_utils_filters_fir
 */
int FIRFilter::getInitialConditionLength(void) const
{
    if (!linit_)
    {
        RTSEIS_ERRMSG("%s", "Class is not yet initialized");
        return -1;
    }
    return order_;
}
/*!
 * @brief Returns a copy of the initial conditions.
 * @param[in] nz   The length of zi.  This must be at least 
 *                 getInitialConditionLength().
 * @param[out] zi  A copy of the initial conditions.  This has dimension
 *                 [nz] however only the first getInitialConditionLength()
 *                 elements are copied.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_filters_fir
 */
int FIRFilter::getInitialConditions(const int nz, double zi[]) const
{
    if (!linit_)
    {   
        RTSEIS_ERRMSG("%s", "Class is not yet initialized");
        return -1;
    }
    int nzLen = getInitialConditionLength();
    if (nzLen > nz)
    {
        RTSEIS_ERRMSG("nz = %d must be at least %d", nz, nzLen);
        return -1;
    }
    if (nzLen > 0)
    {
        if (zi == nullptr)
        {
            RTSEIS_ERRMSG("%s", "zi is NULL");
            return -1;
        }
        ippsCopy_64f(zi_, zi, nzLen);
    }
    return 0;
}
