#include <stdio.h>
#include <stdlib.h>
#include <ipps.h>
#include <cmath>
#define RTSEIS_LOGGING 1
#define IPPS_CORE_SRC 1
#include "rtseis/utils/filters.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utils::Filters;

IIRIIRFilter::IIRIIRFilter(void)
{
    clear();
    return;
}

IIRIIRFilter::~IIRIIRFilter(void)
{
    clear();
    return;
}

IIRIIRFilter::IIRIIRFilter(const IIRIIRFilter &iiriir)
{
    *this = iiriir;
    return;
}

IIRIIRFilter& IIRIIRFilter::operator=(const IIRIIRFilter &iiriir)
{
    if (&iiriir == this){return *this;}
    clear();
    if (!iiriir.isInitialized_){return *this;}
    initialize(iiriir.nbRef_, iiriir.bRef_,
               iiriir.naRef_, iiriir.aRef_,
               iiriir.getPrecision());
    // Now copy the filter states
    if (bufferSize_ > 0)
    {
        Ipp8u *pBufIn = static_cast<Ipp8u *> (iiriir.pBuf_);
        Ipp8u *pBufOut = static_cast<Ipp8u *> (pBuf_);
        ippsCopy_8u(pBufIn, pBufOut, bufferSize_);
    }
    // Copy the initial conditions
    if (nwork_ > 0)
    {
        ippsCopy_64f(iiriir.zi_, zi_, nwork_);
        if (isDoublePrecision())
        {
            Ipp64f *dlysrcIn  = static_cast<Ipp64f *> (iiriir.dlysrc_);
            Ipp64f *dlysrcOut = static_cast<Ipp64f *> (dlysrc_);
            ippsCopy_64f(dlysrcIn, dlysrcOut, nwork_);
        }
        else
        {
            Ipp32f *dlysrcIn  = static_cast<Ipp32f *> (iiriir.dlysrc_);
            Ipp32f *dlysrcOut = static_cast<Ipp32f *> (dlysrc_);
            ippsCopy_32f(dlysrcIn, dlysrcOut, nwork_);
        }
    }
    lhaveZI_ = iiriir.lhaveZI_;
    return *this;
}

void IIRIIRFilter::clear(void)
{
    if (pTaps_ != nullptr){ippsFree(pTaps_);}
    if (dlysrc_ != nullptr){ippsFree(dlysrc_);}
    if (pBuf_ != nullptr){ippsFree(pBuf_);}
    if (bRef_ != nullptr){ippsFree(bRef_);}
    if (aRef_ != nullptr){ippsFree(aRef_);}
    if (zi_ != nullptr){ippsFree(zi_);}
    setPrecision(RTSeis::Precision::DOUBLE);
    toggleRealTime(false);
    pState_ = nullptr;
    pTaps_ = nullptr;
    dlysrc_ = nullptr;
    pBuf_ = nullptr;
    bRef_ = nullptr;
    aRef_ = nullptr;
    zi_ = nullptr;
    nwork_ = 0;
    bufferSize_ = 0;
    order_ = 0;
    nbRef_ = 0;
    naRef_ = 0;
    lhaveZI_ = false;
    isInitialized_ = false;
    return;
}

int IIRIIRFilter::initialize(const int nb, const double b[],
                             const int na, const double a[],
                             const RTSeis::Precision precision)
{
    clear();
    // Check inputs
    if (nb < 1 || na < 1 || b == nullptr || a == nullptr)
    {
        if (nb < 1){RTSEIS_ERRMSG("%s", "No b coefficients");}
        if (na < 1){RTSEIS_ERRMSG("%s", "No a coefficients");}
        if (b == nullptr){RTSEIS_ERRMSG("%s", "a is NULL");}
        if (a == nullptr){RTSEIS_ERRMSG("%s", "b is NULL");}
    }
    // Copy basics
    nbRef_ = nb;
    naRef_ = na;
    order_ = std::max(nb, na) - 1;
    nwork_ = std::max(32, order_);
    bRef_ = ippsMalloc_64f(nb);
    ippsCopy_64f(b, bRef_, nb);
    aRef_ = ippsMalloc_64f(nb); 
    ippsCopy_64f(a, aRef_, na);
    zi_ = ippsMalloc_64f(nwork_);
    ippsZero_64f(zi_, nwork_);
    if (precision == RTSeis::Precision::DOUBLE)
    {
        // Workspace query
        IppStatus status = ippsIIRIIRGetStateSize_64f(order_, &bufferSize_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to get buffer size");
            clear();
            return -1;
        }
        // Set the workspace
        Ipp64f *pTaps = ippsMalloc_64f(2*(order_+1));
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        Ipp64f *dlysrc = ippsMalloc_64f(nwork_);
        ippsZero_64f(dlysrc, nwork_);
        // Copy the filter taps 
        ippsZero_64f(pTaps, 2*(order_+1));
        ippsCopy_64f(b, &pTaps[0],        nb);
        ippsCopy_64f(a, &pTaps[order_+1], na);
        // Initialize the IIR state
        IppsIIRState_64f *pIIRState = nullptr;
        status = ippsIIRIIRInit_64f(&pIIRState, pTaps, order_, NULL, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to initialize IIRIIR state");
            clear();
            return -1;
        }
        // Set pointers
        pState_ = pIIRState;
        pTaps_ = pTaps;
        pBuf_ = pBuf;
        dlysrc_ = dlysrc;
    }
    else
    {
        // Workspace query 
        IppStatus status = ippsIIRIIRGetStateSize_32f(order_, &bufferSize_);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to get buffer size");
            clear();
            return -1;
        }
        // Set the workspace
        Ipp32f *pTaps = ippsMalloc_32f(2*(order_+1));
        Ipp8u *pBuf = ippsMalloc_8u(bufferSize_);
        Ipp32f *dlysrc = ippsMalloc_32f(nwork_);
        ippsZero_32f(dlysrc, nwork_);
        // Copy the filter taps 
        ippsZero_32f(pTaps, 2*(order_+1));
        ippsConvert_64f32f(b, &pTaps[0],        nb);
        ippsConvert_64f32f(a, &pTaps[order_+1], na);
        // Initialize the IIR state 
        IppsIIRState_32f *pIIRState = nullptr;
        status = ippsIIRIIRInit_32f(&pIIRState, pTaps, order_, NULL, pBuf);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to initialize IIRIIR state");
            clear();
            return -1;
        }
        // Set pointers
        pState_ = pIIRState;
        pTaps_ = pTaps;
        pBuf_ = pBuf;
        dlysrc_ = dlysrc;
    }
    setPrecision(precision);
    toggleRealTime(false);
    lhaveZI_ = false;
    isInitialized_ = true;
    return 0;
}

int IIRIIRFilter::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    lhaveZI_ = false;
    resetInitialConditions();
    int nzRef = getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef){RTSEIS_ERRMSG("nz=%d should equal %d", nz, nzRef);}
        if (zi == nullptr){RTSEIS_ERRMSG("%s", "zi is NULL");}
        return -1;
    }
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
    lhaveZI_ = true;
    return 0;
}

int IIRIIRFilter::apply(const int n, const double x[], double y[])
{
    if (n <= 0){return 0;}
    if (!isInitialized_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
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
    // Get handle to filter state
    IppsIIRState_64f *pState = static_cast<IppsIIRState_64f *> (pState_);
    // Set a delay line if the user desires it.  Note, the initialization
    // sets a NULL delay line.
    IppStatus status;
    if (lhaveZI_)
    {
        Ipp64f *dlysrc = static_cast<Ipp64f *> (dlysrc_);
        status = ippsIIRIIRSetDlyLine_64f(pState, dlysrc);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to set delay line");
            return -1; 
        }
    }
    // Apply the filter
    status = ippsIIRIIR_64f(x, y, n, pState);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
    }
    // Undo the action of setting a delay line
    if (lhaveZI_){ippsIIRIIRSetDlyLine_64f(pState, NULL);}
    return 0;
}

int IIRIIRFilter::apply(const int n, const float x[], float y[])
{
    if (n <= 0){return 0;} 
    if (!isInitialized_)
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1;
    }
    if (isFloatPrecision())
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
    // Get handle to filter state
    IppsIIRState_32f *pState = static_cast<IppsIIRState_32f *> (pState_);
    // Set a delay line if the user desires it.  Note, the initialization sets
    // a NULL delay line.
    IppStatus status;
    if (lhaveZI_)
    {
        Ipp32f *dlysrc = static_cast<Ipp32f *> (dlysrc_);
        status = ippsIIRIIRSetDlyLine_32f(pState, dlysrc);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to set delay line");
            return -1;
        }
    }
    // Apply the filter
    status = ippsIIRIIR_32f(x, y, n, pState);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
    }
    // Undo the action of setting a delay line
    if (lhaveZI_){ippsIIRIIRSetDlyLine_32f(pState, NULL);}
    return 0;
}

int IIRIIRFilter::resetInitialConditions(void)
{
    if (!isInitialized_)
    {
        RTSEIS_ERRMSG("%s", "Class is not initialized");
        return -1; 
    }
    return 0;
}

int IIRIIRFilter::getInitialConditionLength(void) const
{
    if (!isInitialized_)
    {
        RTSEIS_ERRMSG("%s", "Class is not initialized");
        return -1;
    }
    return order_;
}

int IIRIIRFilter::getFilterOrder(void) const
{
    return order_;
}
