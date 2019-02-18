#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <vector>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/modules/classicSTALTA.hpp"
#include "rtseis/utilities/filters.hpp"
#include <ipps.h>

using namespace RTSeis::Modules;

class ClassicSTALTA::ClassicSTALTAImpl
{
    public:
        /// Default constructor
        ClassicSTALTAImpl(void)
        {
            return;
        }
        /// Copy constructor
        ClassicSTALTAImpl(const ClassicSTALTAImpl &stalta)
        { 
            *this = stalta;
            return;
        }
        /// Classic destructor
        ~ClassicSTALTAImpl(void) 
        {
            clear();
            return;
        }
        /// Copy assignment operator
        ClassicSTALTAImpl& operator=(const ClassicSTALTAImpl &stalta)
        {
            if (&stalta == this){return *this;} 
            clear();
            if (!stalta.linit_){return *this;}
            firNum_ = stalta.firNum_;
            firDen_ = stalta.firDen_;
            nsta_ = stalta.nsta_;
            nlta_ = stalta.nlta_;
            chunkSize_ = stalta.chunkSize_;
            mode_ = stalta.mode_;
            precision_ = stalta.precision_;
            linit_ = stalta.linit_;
            if (chunkSize_ > 0)
            {
                if (precision_ == RTSeis::Precision::DOUBLE)
                {
                    x264f_   = ippsMalloc_64f(chunkSize_);
                    ippsCopy_64f(stalta.x264f_, x264f_, chunkSize_);
                    ynum64f_ = ippsMalloc_64f(chunkSize_);
                    ippsCopy_64f(stalta.ynum64f_, ynum64f_, chunkSize_);
                    yden64f_ = ippsMalloc_64f(chunkSize_);
                    ippsCopy_64f(stalta.yden64f_, yden64f_, chunkSize_);
                }
                else
                {
                    x232f_   = ippsMalloc_32f(chunkSize_);
                    ippsCopy_32f(stalta.x232f_, x232f_, chunkSize_);
                    ynum32f_ = ippsMalloc_32f(chunkSize_);
                    ippsCopy_32f(stalta.ynum32f_, ynum32f_, chunkSize_);
                    yden32f_ = ippsMalloc_32f(chunkSize_);
                    ippsCopy_32f(stalta.yden32f_, yden32f_, chunkSize_);
                }
            }
            return *this;
        }
        /// Releases memory on the module
        void clear(void)
        {
            firNum_.clear();
            firDen_.clear();
            if (x264f_ != nullptr){ippsFree(x264f_);}
            if (ynum64f_ != nullptr){ippsFree(ynum64f_);}
            if (yden64f_ != nullptr){ippsFree(yden64f_);}
            if (x232f_ != nullptr){ippsFree(x232f_);}
            if (ynum32f_ != nullptr){ippsFree(ynum32f_);}
            if (yden32f_ != nullptr){ippsFree(yden32f_);}
            x264f_ = nullptr;
            ynum64f_ = nullptr;
            yden64f_ = nullptr;
            x232f_ = nullptr;
            ynum32f_ = nullptr;
            yden32f_ = nullptr;
            chunkSize_ = 1024;
            nsta_ = 0;
            nlta_ = 0;
            mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
            precision_ = RTSeis::Precision::DOUBLE;
            linit_ = false;
            return;
        }
        //--------------------------------------------------------------------//
        int initialize(const int nsta, const int nlta,
                       const int chunkSize,
                       const RTSeis::ProcessingMode mode,
                       const RTSeis::Precision precision)
        {
            clear();
            // Set constants
            RTSeis::ProcessingMode modeRT = RTSeis::ProcessingMode::REAL_TIME;
            nsta_ = nsta;
            nlta_ = nlta;
            chunkSize_ = chunkSize;
            // Set the short-term averaging filter coefficients
            Ipp64f *xsta = ippsMalloc_64f(nsta_);
            double xdiv = 1.0/static_cast<double> (nsta_);
            ippsSet_64f(xdiv, xsta, nsta_);
            int ierr = firNum_.initialize(nsta_, xsta, modeRT, precision,
                 RTSeis::Utilities::Filters::FIRFilter::Implementation::DIRECT);
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Failed to set numerator FIR filter");
                ippsFree(xsta);
                clear();
                return -1;
            }
            // Set the initial conditions to 0
            ippsSet_64f(0, xsta, nsta_);
            ierr = firNum_.setInitialConditions(nsta_-1, xsta); 
            ippsFree(xsta);
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Failed to set numerator ics");
                clear();
                return -1;
            }
            // Set the long-term averaging filter coefficients
            Ipp64f *xlta = ippsMalloc_64f(nlta_);
            xdiv = 1.0/static_cast<double> (nlta_);
            ippsSet_64f(xdiv, xlta, nlta_);
            ierr = firDen_.initialize(nlta_, xlta, modeRT, precision,
                 RTSeis::Utilities::Filters::FIRFilter::Implementation::DIRECT);
            if (ierr != 0)
            {
                 RTSEIS_ERRMSG("%s", "Failed to initialize denominator");
                 ippsFree(xlta);
                 clear();
                 return -1;
            }
            // Set the initial conditions to something large
            double xset = DBL_MAX/static_cast<double> (nlta_)/4.0; 
            if (precision == RTSeis::Precision::FLOAT)
            {
                 xset = FLT_MAX/static_cast<double> (nlta_)/4.0;
            }
            ippsSet_64f(xset, xlta, nlta_);
            ierr = firDen_.setInitialConditions(nlta_-1, xlta);
            ippsFree(xlta);
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Failed to set denominator ics");
                clear();
                return -1;
            }
            // Initialize workspace
            if (precision == RTSeis::Precision::DOUBLE)
            {
                x264f_   = ippsMalloc_64f(chunkSize_);
                ynum64f_ = ippsMalloc_64f(chunkSize_);
                yden64f_ = ippsMalloc_64f(chunkSize_);
            } 
            else
            {
                x232f_   = ippsMalloc_32f(chunkSize_);
                ynum32f_ = ippsMalloc_32f(chunkSize_);
                yden32f_ = ippsMalloc_32f(chunkSize_);
            }
            mode_ = mode;
            precision_ = precision;
            linit_ = true;
            return 0;
        } 
        /// Determines if the module is initialized
        bool isInitialized(void) const
        {
            return linit_;
        }
        /// Gets length of the numerator initial conditons
        int getNumeratorInitialConditionLength(void) const
        {
            return firNum_.getInitialConditionLength();
        }
        /// Gets length of the denominator initial conditons
        int getDenominatorInitialConditionLength(void) const
        {
            return firDen_.getInitialConditionLength();
        }
        /// Sets the initial conditions
        int setInitialConditions(const int nzNum, const double zNum[],
                                 const int nzDen, const double zDen[])
        {
            resetInitialConditions(); 
            // Set the initial conditions
            int nzNumRef = getNumeratorInitialConditionLength();
            int nzDenRef = getDenominatorInitialConditionLength();
            if (nzNumRef != nzNum){RTSEIS_ERRMSG("%s", "Shouldn't happen");}
            if (nzDenRef != nzDen){RTSEIS_ERRMSG("%s", "Shouldn't happen");}
            // Set numerator initial conditions
            int ierr = firNum_.setInitialConditions(nzNumRef, zNum);
            if (ierr != 0)
            {
                 RTSEIS_ERRMSG("%s", "Failed to set numerator ics");
                 return -1;
            }
            ierr = firDen_.setInitialConditions(nzDenRef, zDen);
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Failed to set denominator ics");
                return -1;
            }
            return 0;
        }
        /// Resets the initial conditions
        int resetInitialConditions(void)
        {
            firNum_.resetInitialConditions();
            firDen_.resetInitialConditions();
            return 0;
        }
        /// Applies the STA/LTA
        int apply(const int nx, const double x[], double y[])
        {
            if (nx <= 0){return 0;} // Nothing to do
            for (int i=0; i<nx; i=i+chunkSize_)
            {
                int nloc = std::min(chunkSize_, nx - i);
                // Compute the squared signal
                ippsSqr_64f(&x[i], x264f_, nloc);
                // Compute the numerator average 
                int ierr = firNum_.apply(nloc, x264f_, ynum64f_);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to filter numerator");
                    return -1;
                }
                // Compute the denominator average
                ierr = firDen_.apply(nloc, x264f_, yden64f_);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to filter denominator");
                    return -1;
                }
                // Pointwise division
                IppStatus status = ippsDiv_64f(yden64f_, ynum64f_, &y[i], nloc);
                if (status != ippStsNoErr)
                {
                    // Division by zero error can be handled.  This means the
                    // the numerator must be 0 as well - so we force the division
                    // 0/0 = 0.  The head is that a dead signal won't trigger.
                    if (status == ippStsDivByZero)
                    {
                       RTSEIS_WARNMSG("%s", "Division by zero detected");
                       // if |yden| < DBL_MIN then y = 0
                       status = ippsThreshold_LTAbsVal_64f(yden64f_, &y[i],
                                                           nloc,
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
            if (mode_ == RTSeis::ProcessingMode::POST_PROCESSING)
            {
                resetInitialConditions();
            }
            return 0;
        }
        /// Applies the STA/LTA
        int apply(const int nx, const float x[], float y[])
        {
            if (nx <= 0){return 0;} // Nothing to do
            for (int i=0; i<nx; i=i+chunkSize_)
            {
                int nloc = std::min(chunkSize_, nx - i);
                // Compute the squared signal
                ippsSqr_32f(&x[i], x232f_, nloc);
                // Compute the numerator average 
                int ierr = firNum_.apply(nloc, x232f_, ynum32f_);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to filter numerator");
                    return -1;
                }
                // Compute the denominator average
                ierr = firDen_.apply(nloc, x232f_, yden32f_);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to filter denominator");
                    return -1;
                }
                // Pointwise division
                IppStatus status = ippsDiv_32f(yden32f_, ynum32f_, &y[i], nloc);
                if (status != ippStsNoErr)
                {
                    // Division by zero error can be handled.  This means the
                    // the numerator must be 0 as well - so we force the division
                    // 0/0 = 0.  The head is that a dead signal won't trigger.
                    if (status == ippStsDivByZero)
                    {
                        RTSEIS_WARNMSG("%s", "Division by zero detected");
                        // if |yden| < DBL_MIN then y = 0
                        status = ippsThreshold_LTAbsVal_32f(yden32f_, &y[i],
                                                            nloc,
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
            if (mode_ == RTSeis::ProcessingMode::POST_PROCESSING)
            {
                resetInitialConditions();
            }
            return 0;
        }
    private:
        /// Tabulates the numerator short-term average
        RTSeis::Utilities::Filters::FIRFilter firNum_;
        /// Tabulates the denominator long-term average
        RTSeis::Utilities::Filters::FIRFilter firDen_;
        /// The characteristic function.  This has dimension [chunkSize_].
        Ipp64f *x264f_ = nullptr;
        /// Workspace array for holding numerator.  This has
        /// dimension [chunkSize_].
        Ipp64f *ynum64f_ = nullptr;
        /// Workspace array for hodling denominator.  This has
        /// dimension [chunkSize_].
        Ipp64f *yden64f_ = nullptr;
        /// The characteristic function.  This has dimension [chunkSize_].
        Ipp32f *x232f_ = nullptr;
        /// Workspace array for holding numerator.  This has
        /// dimension [chunkSize_].
        Ipp32f *ynum32f_ = nullptr;
        /// Workspace array for hodling denominator.  This has
        /// dimension [chunkSize_].
        Ipp32f *yden32f_ = nullptr;
        /// The number of points in the STA window
        int nsta_ = 0;
        /// The number of points in the LTA window
        int nlta_ = 0;
        /// Workspace for numerator and denominators 
        int chunkSize_ = 1024;
        /// The processing mode
        RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        /// The precision of th emodule
        RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
        /// Flag indicating the class is initialized
        bool linit_ = false;
};

//============================================================================//

ClassicSTALTAParameters::ClassicSTALTAParameters(void)
{
    return;
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
    precision_ = parameters.precision_;
    processingMode_ = parameters.processingMode_;
    nsta_ = parameters.nsta_;
    nlta_ = parameters.nlta_;
    chunkSize_ = parameters.chunkSize_;
    isValid_ = parameters.isValid_;
    return *this;
}

ClassicSTALTAParameters::ClassicSTALTAParameters(
    const int nsta, const int nlta,
    const RTSeis::ProcessingMode mode,
    const RTSeis::Precision prec) :
    chunkSize_(1024)
{
    // Set the long-term and short-term parameters
    int ierr = setShortTermAndLongTermWindowSize(nsta, nlta);
    if (ierr != 0)
    {
        clear();
        return;
    }
    setProcessingMode(mode);
    precision_ = prec;
    // Validate
    validate_();
    return;
}

ClassicSTALTAParameters::ClassicSTALTAParameters(
    const double staWin, const double ltaWin, const double dt,
    const RTSeis::ProcessingMode mode,
    const RTSeis::Precision prec) :
    chunkSize_(1024)
{
    // Check parameters
    int ierr = setShortTermAndLongTermWindowSize(staWin, ltaWin, dt);
    if (ierr != 0)
    {
        clear();
        return;
    }
    setProcessingMode(mode);
    precision_ = prec;
    // Validate
    validate_();
    return;
}

ClassicSTALTAParameters::ClassicSTALTAParameters(
    const int nsta, const int nlta,
    const size_t chunkSize,
    const RTSeis::ProcessingMode mode,
    const RTSeis::Precision prec)
{
    // Set the long-term and short-term parameters
    int ierr = setShortTermAndLongTermWindowSize(nsta, nlta);
    if (ierr != 0)
    {
        clear();
        return;
    }
    ierr = setChunkSize(chunkSize);
    if (ierr != 0)
    {
        clear();
        return;
    }
    setProcessingMode(mode);
    precision_ = prec;
    // Validate
    validate_();
    return;
}

ClassicSTALTAParameters::ClassicSTALTAParameters(
    const double staWin, const double ltaWin, const double dt, 
    const size_t chunkSize,
    const RTSeis::ProcessingMode mode,
    const RTSeis::Precision prec)
{
    // Check parameters
    int ierr = setShortTermAndLongTermWindowSize(staWin, ltaWin, dt);
    if (ierr != 0)
    {
        clear();
        return;
    }
    ierr = setChunkSize(chunkSize);
    if (ierr != 0)
    {
        clear();
        return;
    }
    setProcessingMode(mode);
    precision_ = prec;
    // Validate
    validate_();
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
    chunkSize_ = 1024;
    precision_ = defaultPrecision_;
    processingMode_ = RTSeis::ProcessingMode::POST_PROCESSING;
    isValid_ = false;
    return;
}

int ClassicSTALTAParameters::setChunkSize(const size_t chunkSize)
{
    if (chunkSize < 1)
    {
        RTSEIS_ERRMSG("%s", "Chunksize must be positive");
        return -1;
    }
    chunkSize_ = chunkSize;
    validate_();
    return 0;
}

size_t ClassicSTALTAParameters::getChunkSize(void) const
{
    return chunkSize_;
}

int ClassicSTALTAParameters::setShortTermAndLongTermWindowSize(
    const int nsta, const int nlta)
{
    if (nsta < 1)
    {
        RTSEIS_ERRMSG("STA window=%d samples must be at least 1", nsta);
        return -1;
    }
    if (nlta <= nsta)
    {
        RTSEIS_ERRMSG("LTA window=%d samples must be greater than %d",
                      nlta, nsta);
        return -1;
    }
    nsta_ = nsta;
    nlta_ = nlta;
    validate_();
    return 0;
}

int ClassicSTALTAParameters::setShortTermAndLongTermWindowSize(
    const double staWin, const double ltaWin, const double dt)
{
    if (dt <= 0)
    {
        RTSEIS_ERRMSG("dt=%lf must be postiive", dt);
        return -1;
    }
    if (staWin < 0)
    {
        RTSEIS_ERRMSG("STA window length=%lf (s) must be at least %lf (s) ",
                      staWin, dt);
        return -1;
    }
    if (ltaWin <= staWin + dt/2)
    {
        RTSEIS_ERRMSG("LTA window length=%lf (s) must be at least %lf (s) ",
                      ltaWin, staWin + dt);
        return -1;
    }
    int nsta = static_cast<int> (staWin/dt + 0.5) + 1;
    int nlta = static_cast<int> (ltaWin/dt + 0.5) + 1;
    if (nlta <= nsta)
    {
        RTSEIS_ERRMSG("%s", "Algorithmic failure");
        return -1;
    }
    nsta_ = nsta;
    nlta_ = nlta;
    validate_();
    return 0;
}

int ClassicSTALTAParameters::getLongTermWindowSize(void) const
{
    return nlta_;
}

int ClassicSTALTAParameters::getShortTermWindowSize(void) const
{
    return nsta_;
}

void ClassicSTALTAParameters::setProcessingMode(
    const RTSeis::ProcessingMode mode)
{
    processingMode_ = mode;
    validate_();
    return;
}

RTSeis::ProcessingMode ClassicSTALTAParameters::getProcessingMode(void) const
{
    return processingMode_;
}

RTSeis::Precision ClassicSTALTAParameters::getPrecision(void) const
{
    return precision_;
}

bool ClassicSTALTAParameters::isValid(void) const
{
    return isValid_;
}

void ClassicSTALTAParameters::validate_(void)
{
    isValid_ = false;
    if (chunkSize_ < 1){return;}
    if (nsta_ < 1){return;}
    if (nlta_ <= nsta_){return;}
    if (getPrecision() != RTSeis::Precision::DOUBLE &&
        getPrecision() != RTSeis::Precision::FLOAT){return;}
    isValid_ = true;
    return;
}

//============================================================================//
//                                 End Parameters                             //
//============================================================================//

ClassicSTALTA::ClassicSTALTA(void) :
    pSTALTA_(new ClassicSTALTAImpl())
{
    clear();
}

ClassicSTALTA::ClassicSTALTA(const ClassicSTALTA &cstalta)
{
    *this = cstalta;
}

ClassicSTALTA::~ClassicSTALTA(void)
{
    clear();
    return;
}

void ClassicSTALTA::clear(void)
{
    pSTALTA_->clear();
    return;
}

ClassicSTALTA::ClassicSTALTA(const ClassicSTALTAParameters &parameters) :
    pSTALTA_(new ClassicSTALTAImpl())
{
    clear();
    if (!parameters.isValid())
    {
        RTSEIS_ERRMSG("%s", "Parameters are not valid");
        return;
    }
    int nsta = parameters.getShortTermWindowSize();
    int nlta = parameters.getLongTermWindowSize();
    int chunkSize = static_cast<int> (parameters.getChunkSize());
    RTSeis::Precision precision = parameters.getPrecision();
    RTSeis::ProcessingMode mode = parameters.getProcessingMode();
    int ierr = pSTALTA_->initialize(nsta, nlta, chunkSize, mode, precision);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "FAiled to initialize module");
        return;
    }
    return;
}

ClassicSTALTA& ClassicSTALTA::operator=(const ClassicSTALTA &cstalta)
{
    if (&cstalta == this){return *this;}
    if (pSTALTA_){pSTALTA_->clear();}
    pSTALTA_ = std::unique_ptr<ClassicSTALTAImpl>
              (new ClassicSTALTAImpl(*cstalta.pSTALTA_));
    return *this;
}

int ClassicSTALTA::setInitialConditions(const int nzNum, const double zNum[],
                                        const int nzDen, const double zDen[])
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module not initialized");
        return -1;
    }
    int nzNumRef = pSTALTA_->getNumeratorInitialConditionLength();
    int nzDenRef = pSTALTA_->getDenominatorInitialConditionLength(); 
    if (nzNumRef != nzNum)
    {
        RTSEIS_ERRMSG("nzNum=%d must equal %d", nzNum, nzNumRef);
        return -1;
    }
    if (nzDenRef != nzDen)
    {
        RTSEIS_ERRMSG("nzDen=%d must equal %d", nzDen, nzDenRef);
        return -1;
    }
    if (nzNum > 0 && zNum == nullptr)
    {
        RTSEIS_ERRMSG("%s", "zNum is NULL");
        return -1;
    }
    if (nzDen > 0 && zDen == nullptr)
    {
        RTSEIS_ERRMSG("%s", "zDen is NULL");
        return -1;
    }
    int ierr = pSTALTA_->setInitialConditions(nzNum, zNum, nzDen, zDen);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to set initial conditions");
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
    pSTALTA_->resetInitialConditions();
    return 0;
}

int ClassicSTALTA::getNumeratorInitialConditionLength(void) const
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module not initialized");
        return -1;
    }
    int nz = pSTALTA_->getNumeratorInitialConditionLength();
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
    int nz = pSTALTA_->getDenominatorInitialConditionLength();
    if (nz < 0){RTSEIS_ERRMSG("%s", "Internal error");}
    return nz; 
}

bool ClassicSTALTA::isInitialized(void) const
{
    return pSTALTA_->isInitialized();
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
    int ierr = pSTALTA_->apply(nx, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
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
    int ierr = pSTALTA_->apply(nx, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
    }
    return 0;
}
