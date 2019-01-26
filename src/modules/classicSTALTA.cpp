#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <vector>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/classicSTALTA.hpp"
#include "rtseis/utils/filters.hpp"
#include <ipps.h>

using namespace RTSeis::Modules;

ClassicSTALTAParameters::ClassicSTALTAParameters(void)
{
    clear();
}

ClassicSTALTAParameters::ClassicSTALTAParameters(
    const ClassicSTALTAParameters &parameters)
{
    *this = parameters;
    return;
}

ClassicSTALTAParameters&
ClassicSTALTAParameters::operator=(const ClassicSTALTAParameters &parameters)
{
    if (&parameters == this){return *this;}
    clear();
    if (!parameters.isInitialized()){return *this;}
    precision_ = parameters.precision_;
    isRealTime_ = parameters.isRealTime_;
    nsta_ = parameters.nsta_;
    nlta_ = parameters.nlta_;
    chunkSize_ = parameters.chunkSize_;
    isInitialized_ = parameters.isInitialized_;
    return *this;
}

ClassicSTALTAParameters::ClassicSTALTAParameters(
    const int nsta, const int nlta,
    const bool lrt,
    const RTSeis::Precision prec) :
    chunkSize_(1024)
{
    size_t chunk = chunkSize_;
    ClassicSTALTAParameters(nsta, nlta, chunk, lrt, prec);
    return;
}

ClassicSTALTAParameters::ClassicSTALTAParameters(
    const double staWin, const double ltaWin, const double dt,
    const bool lrt,
    const RTSeis::Precision prec) :
    chunkSize_(1024)
{
    size_t chunk = chunkSize_;
    ClassicSTALTAParameters(staWin, ltaWin, dt, chunk, lrt, prec);
    return;
}

ClassicSTALTAParameters::ClassicSTALTAParameters(
    const int nsta, const int nlta,
    const size_t chunkSize,
    const bool lrt,
    const RTSeis::Precision prec)
{
    clear();
    // Check parameters
    if (nsta < 1)
    {   
        RTSEIS_ERRMSG("STA window=%d samples must be at least 1", nsta);
        return;
    }   
    if (nlta <= nsta)
    {   
        RTSEIS_ERRMSG("LTA window=%d samples must be greater than %d",
                      nlta, nsta);
        return;
    }   
    if (chunkSize < 1)
    {
        RTSEIS_ERRMSG("%s", "Chunksize must be positive");
        return;
    }
    // Initialize
    precision_ = prec;
    isRealTime_ = lrt;
    nsta_ = nsta;
    nlta_ = nlta;
    chunkSize_ = chunkSize;
    isInitialized_ = true;
    return;
}



ClassicSTALTAParameters::ClassicSTALTAParameters(
    const double staWin, const double ltaWin, const double dt, 
    const size_t chunkSize,
    const bool lrt,
    const RTSeis::Precision prec)
{
    clear();
    // Check parameters
    if (chunkSize < 1)
    {
        RTSEIS_ERRMSG("%s", "Chunksize must be positive");
        return;
    }
    if (dt <= 0)
    {
        RTSEIS_ERRMSG("dt=%lf must be postiive", dt) 
        return;
    }
    if (staWin < dt) 
    {
        RTSEIS_ERRMSG("STA window length=%lf (s) must be at least %lf (s) ",
                      staWin, dt);
        return;
    }
    if (ltaWin <= staWin + dt) 
    {
        RTSEIS_ERRMSG("LTA window length=%lf (s) must be at least %lf (s) ",
                      ltaWin, staWin + dt);
        return;
    }
    int nsta = static_cast<int> (staWin/dt + 0.5);
    int nlta = static_cast<int> (ltaWin/dt + 0.5);
    if (nlta <= nsta)
    {
        RTSEIS_ERRMSG("%s", "Algorithmic failure");
        return;
    }
    // Initialize
    precision_ = prec;
    isRealTime_ = lrt;
    nsta_ = nsta;
    nlta_ = nlta;
    chunkSize_ = chunkSize;
    isInitialized_ = true;
    return;
}

ClassicSTALTAParameters::~ClassicSTALTAParameters(void)
{
    clear();
}

void ClassicSTALTAParameters::clear(void)
{
    nsta_ = 0;
    nlta_ = 0;
    precision_ = defaultPrecision_;
    isRealTime_ = false;
    isInitialized_ = false;
    return;
}

size_t ClassicSTALTAParameters::getChunkSize(void) const
{
    return chunkSize_;
}

int ClassicSTALTAParameters::getLongTermWindowSize(void) const
{
    return nlta_;
}

int ClassicSTALTAParameters::getShortTermWindowSize(void) const
{
    return nsta_;
}

bool ClassicSTALTAParameters::isInitialized(void) const
{
    return isInitialized_;
}

bool ClassicSTALTAParameters::isRealTime(void) const
{
    return isRealTime_;
}

RTSeis::Precision ClassicSTALTAParameters::getPrecision(void) const
{
    return precision_;
}

//============================================================================//
//                                 End Parameters                             //
//============================================================================//

ClassicSTALTA::ClassicSTALTA(void)
{
    clear();
}

ClassicSTALTA::ClassicSTALTA(const ClassicSTALTA &cstalta)
{
    *this = cstalta;
}

ClassicSTALTA::ClassicSTALTA(const ClassicSTALTAParameters &parameters)
{
    const bool lrt = true;
    clear();
    if (!parameters.isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Parameters are not initialized");
        return;
    }
    parms_ = parameters;
    // Set the short-term averaging filter coefficients
    int nsta = parms_.getShortTermWindowSize();
    Ipp64f *xsta = ippsMalloc_64f(nsta);
    double xdiv = 1.0/static_cast<double> (nsta);
    ippsSet_64f(xdiv, xsta, nsta);
    int ierr = firNum_.initialize(nsta, xsta, lrt, parms_.getPrecision(),
                    RTSeis::Utils::Filters::FIRFilter::Implementation::DIRECT);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to set numerator FIR filter");
        ippsFree(xsta);
        clear();
        return;
    } 
    // Set the initial conditions to 0
    ippsSet_64f(0, xsta, nsta);
    ierr = firNum_.setInitialConditions(nsta-1, xsta); 
    ippsFree(xsta);
    if (ierr != 0)
    {   
        RTSEIS_ERRMSG("%s", "Failed to set numerator initial conditions");
        clear();
        return;
    }
    // Set the long-term averaging filter coefficients
    int nlta = parms_.getLongTermWindowSize();
    Ipp64f *xlta = ippsMalloc_64f(nlta);
    xdiv = 1.0/static_cast<double> (nlta);
    ippsSet_64f(xdiv, xlta, nlta);
    ierr = firDen_.initialize(nlta, xlta, lrt, parms_.getPrecision(),
                  RTSeis::Utils::Filters::FIRFilter::Implementation::DIRECT);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize denominator");
        ippsFree(xlta);
        clear();
        return;
    }
    // Set the initial conditions to something large
    double xset = DBL_MAX/static_cast<double> (nlta)/4.0; 
    if (parms_.getPrecision() == RTSeis::Precision::FLOAT)
    {
        xset = FLT_MAX/static_cast<double> (nlta)/4.0;
    }
    ippsSet_64f(xset, xlta, nlta); 
    ierr = firDen_.setInitialConditions(nlta-1, xlta); 
    ippsFree(xlta);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to set numerator initial conditions");
        clear();
        return;
    }
    // Set the workspace
    int chunkSize = static_cast<int> (parms_.getChunkSize());
    if (parms_.getPrecision() == RTSeis::Precision::DOUBLE)
    {
        x2_   = ippsMalloc_64f(chunkSize);
        ynum_ = ippsMalloc_64f(chunkSize);
        yden_ = ippsMalloc_64f(chunkSize);
    }
    else
    {
        x2_   = ippsMalloc_32f(chunkSize);
        ynum_ = ippsMalloc_32f(chunkSize);
        yden_ = ippsMalloc_32f(chunkSize);
    }
    isInitialized_ = true;
    return;
}

ClassicSTALTA& ClassicSTALTA::operator=(const ClassicSTALTA &cstalta)
{
    if (&cstalta == this){return *this;}
    clear();
    if (!cstalta.isInitialized()){return *this;}
    firNum_ = cstalta.firNum_;
    firDen_ = cstalta.firDen_;
    parms_ = cstalta.parms_;
    int chunkSize = static_cast<int> (parms_.getChunkSize());
    if (parms_.getPrecision() == RTSeis::Precision::DOUBLE)
    {
        Ipp64f *x2In  = static_cast<Ipp64f *> (cstalta.x2_);
        Ipp64f *x2Out = ippsMalloc_64f(chunkSize);
        ippsCopy_64f(x2In, x2Out, chunkSize);

        Ipp64f *ynumIn  = static_cast<Ipp64f *> (cstalta.ynum_);
        Ipp64f *ynumOut = ippsMalloc_64f(chunkSize);
        ippsCopy_64f(ynumIn, ynumOut, chunkSize);

        Ipp64f *ydenIn  = static_cast<Ipp64f *> (cstalta.yden_);
        Ipp64f *ydenOut = ippsMalloc_64f(chunkSize);
        ippsCopy_64f(ydenIn, ydenOut, chunkSize);
    }
    else
    {
        Ipp32f *x2In  = static_cast<Ipp32f *> (cstalta.x2_);
        Ipp32f *x2Out = ippsMalloc_32f(chunkSize);
        ippsCopy_32f(x2In, x2Out, chunkSize);

        Ipp32f *ynumIn  = static_cast<Ipp32f *> (cstalta.ynum_);
        Ipp32f *ynumOut = ippsMalloc_32f(chunkSize);
        ippsCopy_32f(ynumIn, ynumOut, chunkSize);

        Ipp32f *ydenIn  = static_cast<Ipp32f *> (cstalta.yden_);
        Ipp32f *ydenOut = ippsMalloc_32f(chunkSize);
        ippsCopy_32f(ydenIn, ydenOut, chunkSize);
    }
    isInitialized_ = cstalta.isInitialized_;
    return *this;
}

ClassicSTALTA::~ClassicSTALTA(void)
{
    clear();
}

void ClassicSTALTA::clear(void)
{
    if (x2_ != nullptr){ippsFree(x2_);}
    if (ynum_ != nullptr){ippsFree(ynum_);}
    if (yden_ != nullptr){ippsFree(yden_);}
    x2_ = nullptr;
    ynum_ = nullptr;
    yden_ = nullptr;
    parms_.clear();
    firNum_.clear();
    firDen_.clear();
    isInitialized_ = false;
    return;
}

int ClassicSTALTA::setInitialConditions(const int nzNum, const double zNum[],
                                        const int nzDen, const double zDen[])
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module not initialized");
        return -1;
    }
    resetInitialConditions();
    if (nzNum != firNum_.getInitialConditionLength())
    {
        RTSEIS_ERRMSG("nzNum=%d must equal %d",
                      nzNum, firNum_.getInitialConditionLength())
        return -1;
    }
    if (nzNum > 0 && zNum == nullptr)
    {
        RTSEIS_ERRMSG("%s", "zNum is NULL");
        return -1;
    }
    if (nzDen != firDen_.getInitialConditionLength())
    {
        RTSEIS_ERRMSG("nzDen=%d must equal %d",
                      nzNum, firDen_.getInitialConditionLength())
        return -1;
    }
    if (nzDen > 0 && zDen == nullptr)
    {
        RTSEIS_ERRMSG("%s", "zDen is NULL");
        return -1;
    }
    // Get a copy of the initial conditions in case a mistake occurs
    int nwork = std::max(nzNum, nzDen);
    std::vector<double> ziNum(nwork);
    std::vector<double> ziDen(nwork);
    firNum_.getInitialConditions(nwork, ziNum.data());
    // Set the initial conditions
    int ierr;
    ierr = firNum_.setInitialConditions(nzNum, zNum);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to set initial numerator conditions");
        if (nzNum > 0){firNum_.setInitialConditions(nzNum, ziNum.data());}
        return -1;
    }
    ierr = firDen_.setInitialConditions(nzDen, zDen);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to set intiial denominator conditions");
        if (nzNum > 0){firNum_.setInitialConditions(nzNum, ziNum.data());}
        if (nzDen > 0){firDen_.setInitialConditions(nzDen, ziDen.data());}
        return -1;
    }
    return 0;
}

int ClassicSTALTA::resetInitialConditions(void)
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module not initialized");
        return -1; 
    }
    firNum_.resetInitialConditions();
    firDen_.resetInitialConditions();
    return 0;
}

int ClassicSTALTA::getNumeratorInitialConditionLength(void) const
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module not initialized");
        return -1;
    }
    int nz = firNum_.getInitialConditionLength();
    if (nz < 0){RTSEIS_ERRMSG("%s", "Internal error");}
    return nz;
}

int ClassicSTALTA::getDenominatorInitialConditionLength(void) const
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module not initialized");
        return -1;
    }
    int nz = firDen_.getInitialConditionLength();
    if (nz < 0){RTSEIS_ERRMSG("%s", "Internal error");}
    return nz; 
}

bool ClassicSTALTA::isInitialized(void) const
{
    return isInitialized_;
}

int ClassicSTALTA::apply(const int nx, const double x[], double y[])
{
    if (nx <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module not initialized");
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1;
    }
    int chunkSize = static_cast<int> (parms_.getChunkSize());
    Ipp64f *x2   = static_cast<Ipp64f *> (x2_);
    Ipp64f *ynum = static_cast<Ipp64f *> (ynum_);
    Ipp64f *yden = static_cast<Ipp64f *> (yden_);
    for (int i=0; i<nx; i=i+chunkSize)
    {
        int nloc = std::min(chunkSize, nx - i);
        // Compute the squared signal
        IppStatus status = ippsSqr_64f(&x[i], x2, nloc); 
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to square signal");
            return -1;
        }
        // Compute the numerator average 
        int ierr = firNum_.apply(nloc, x2, ynum);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to filter numerator");
            return -1;
        }
        // Compute the denominator average
        ierr = firDen_.apply(nloc, x2, yden);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to filter denominator");
            return -1;
        }
        // Pointwise division
        status = ippsDiv_64f(yden, ynum, &y[i], nloc);
        if (status != ippStsNoErr)
        {
            // Division by zero error can be handled.  This means the
            // the numerator must be 0 as well - so we force the division
            // 0/0 = 0.  The head is that a dead signal won't trigger.
            if (status == ippStsDivByZero)
            {
                RTSEIS_WARNMSG("%s", "Division by zero detected");
                // if |yden| < DBL_MIN then y = 0
                status = ippsThreshold_LTAbsVal_64f(yden, y, nloc,
                                                    DBL_MIN, 0.0);
            }
            // Unrecoverable error
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error computing ynum/yden");
                return -1;
            }
        }
    }
    // Reset the initial conditions for post-processing
    if (!parms_.isRealTime())
    {
        resetInitialConditions();
    }
    return 0;
}

int ClassicSTALTA::apply(const int nx, const float x[], float y[])
{
    if (nx <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module not initialized");
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1;
    }
    int chunkSize = static_cast<int> (parms_.getChunkSize());
    Ipp32f *x2   = static_cast<Ipp32f *> (x2_);
    Ipp32f *ynum = static_cast<Ipp32f *> (ynum_);
    Ipp32f *yden = static_cast<Ipp32f *> (yden_);
    for (int i=0; i<nx; i=i+chunkSize)
    {
        int nloc = std::min(chunkSize, nx - i);
        // Compute the squared signal
        IppStatus status = ippsSqr_32f(&x[i], x2, nloc);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to square signal");
            return -1;
        }
        // Compute the numerator average 
        int ierr = firNum_.apply(nloc, x2, ynum);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to filter numerator");
            return -1;
        }
        // Compute the denominator average
        ierr = firDen_.apply(nloc, x2, yden);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to filter denominator");
            return -1;
        }
        // Pointwise division
        status = ippsDiv_32f(yden, ynum, &y[i], nloc);
        if (status != ippStsNoErr)
        {
            // Division by zero error can be handled.  This means the
            // the numerator must be 0 as well - so we force the division
            // 0/0 = 0.  The head is that a dead signal won't trigger.
            if (status == ippStsDivByZero)
            {
                RTSEIS_WARNMSG("%s", "Division by zero detected");
                // if |yden| < DBL_MIN then y = 0
                status = ippsThreshold_LTAbsVal_32f(yden, y, nloc,
                                                    FLT_MIN, 0.0f);
            }
            // Unrecoverable error
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error computing ynum/yden");
                return -1;
            }
        }
    }
    // Reset the initial conditions for post-processing
    if (!parms_.isRealTime())
    {
        resetInitialConditions();
    }
    return 0;
}
