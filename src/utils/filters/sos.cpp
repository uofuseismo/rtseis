#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/utils/filters.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utils::Filters;

class SOSFilter::SOSFilterImpl
{
    public:
        /// Default constructor
        SOSFilterImpl(void)
        {
            return;
        }
        /// Copy constructor
        SOSFilterImpl(const SOSFilterImpl &sos)
        {
            *this = sos;
        }
        /// (Deep) copy operator
        SOSFilterImpl& operator=(const SOSFilterImpl &sos)
        {
            if (&sos == this){return *this;}
            if (!sos.linit_){return *this;}
            // Reinitialize the filter
            initialize(sos.nsections_, sos.bsRef_, sos.asRef_, 
                       sos.mode_, sos.precision_);
            // Now copy the filter states
            if (bufferSize_ > 0)
            {
                ippsCopy_8u(sos.pBuf_, pBuf_, bufferSize_);
            }
            // Copy the initial conditions
            if (nsections_ > 0){ippsCopy_64f(sos.zi_, zi_, 2*nsections_);}
            // And the delay lines
            if (nwork_ > 0)
            {
                if (precision_ == RTSeis::Precision::DOUBLE)
                {
                    ippsCopy_64f(sos.dlySrc64f_, dlySrc64f_, nwork_);
                    ippsCopy_64f(sos.dlyDst64f_, dlyDst64f_, nwork_);
                }
                else
                {
                    ippsCopy_32f(sos.dlySrc32f_, dlySrc32f_, nwork_); 
                    ippsCopy_32f(sos.dlyDst32f_, dlyDst32f_, nwork_);
                }
            }
            return *this;
        }
        /// Default constructor
        ~SOSFilterImpl(void)
        {
            clear();
            return;
        }
        /// Clears the memory off the module
        void clear(void)
        {
            if (pTaps64f_ != nullptr){ippsFree(pTaps64f_);}
            if (dlySrc64f_ != nullptr){ippsFree(dlySrc64f_);}
            if (dlyDst64f_ != nullptr){ippsFree(dlyDst64f_);}
            if (pTaps32f_ != nullptr){ippsFree(pTaps32f_);}
            if (dlySrc32f_ != nullptr){ippsFree(dlySrc32f_);}
            if (dlyDst32f_ != nullptr){ippsFree(dlyDst64f_);}
            if (pBuf_ != nullptr){ippsFree(pBuf_);}
            if (bsRef_ != nullptr){ippsFree(bsRef_);}
            if (asRef_ != nullptr){ippsFree(asRef_);} 
            if (zi_ != nullptr){ippsFree(zi_);}
            pState64f_ = nullptr;
            pTaps64f_ = nullptr;
            dlySrc64f_ = nullptr;
            dlyDst64f_ = nullptr;
            pState32f_ = nullptr;
            pTaps32f_ = nullptr;
            dlySrc32f_ = nullptr;
            dlyDst32f_ = nullptr; 
            pBuf_ = nullptr;
            bsRef_ = nullptr;
            asRef_ = nullptr;
            zi_ = nullptr;
            nsections_ = 0;
            tapsLen_ = 0;
            nwork_ = 0;
            bufferSize_ = 0;
            mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
            precision_ = RTSeis::Precision::DOUBLE;
            linit_ = false;
            return;
        }
        //====================================================================//
        int initialize(const int ns,
                       const double bs[],
                       const double as[],
                       const RTSeis::ProcessingMode mode,
                       const RTSeis::Precision precision)
        {
            clear();
            // Figure out sizes and copy the inputs
            nsections_ = ns;
            tapsLen_ = 6*nsections_;
            nwork_ = std::max(128, 2*nsections_);
            bsRef_ = ippsMalloc_64f(3*nsections_);
            ippsCopy_64f(bs, bsRef_, 3*nsections_);
            asRef_ = ippsMalloc_64f(3*nsections_);
            ippsCopy_64f(as, asRef_, 3*nsections_);
            zi_ = ippsMalloc_64f(2*nsections_);
            ippsZero_64f(zi_, 2*nsections_);
            IppStatus status;
            if (precision == RTSeis::Precision::DOUBLE)
            {
                status = ippsIIRGetStateSize_BiQuad_64f(nsections_,
                                                        &bufferSize_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Failed to get state size");
                    clear();
                    return -1;
                }
                pBuf_ = ippsMalloc_8u(bufferSize_);
                pTaps64f_ = ippsMalloc_64f(tapsLen_);
                for (int i=0; i<nsections_; i++)
                {
                    pTaps64f_[6*i+0] = bs[3*i+0];
                    pTaps64f_[6*i+1] = bs[3*i+1];
                    pTaps64f_[6*i+2] = bs[3*i+2];
                    pTaps64f_[6*i+3] = as[3*i+0];
                    pTaps64f_[6*i+4] = as[3*i+1];
                    pTaps64f_[6*i+5] = as[3*i+2];
                }
                dlySrc64f_ = ippsMalloc_64f(nwork_);
                ippsZero_64f(dlySrc64f_, nwork_);
                dlyDst64f_ = ippsMalloc_64f(nwork_);
                ippsZero_64f(dlyDst64f_, nwork_);
                status = ippsIIRInit_BiQuad_64f(&pState64f_, pTaps64f_,
                                                nsections_,
                                                dlySrc64f_, pBuf_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Failed to initialized biquad filter");
                    clear();
                    return -1;
                }
            }
            else
            {
                status = ippsIIRGetStateSize_BiQuad_32f(nsections_,
                                                        &bufferSize_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Failed to get state size");
                    clear();
                    return -1; 
                }
                pBuf_ = ippsMalloc_8u(bufferSize_);
                pTaps32f_ = ippsMalloc_32f(tapsLen_);
                for (int i=0; i<nsections_; i++)
                {
                    pTaps32f_[6*i+0] = static_cast<float> (bs[3*i+0]);
                    pTaps32f_[6*i+1] = static_cast<float> (bs[3*i+1]);
                    pTaps32f_[6*i+2] = static_cast<float> (bs[3*i+2]);
                    pTaps32f_[6*i+3] = static_cast<float> (as[3*i+0]);
                    pTaps32f_[6*i+4] = static_cast<float> (as[3*i+1]);
                    pTaps32f_[6*i+5] = static_cast<float> (as[3*i+2]);
                }
                dlySrc32f_ = ippsMalloc_32f(nwork_);
                ippsZero_32f(dlySrc32f_, nwork_);
                dlyDst32f_ = ippsMalloc_32f(nwork_);
                ippsZero_32f(dlyDst32f_, nwork_);
                status = ippsIIRInit_BiQuad_32f(&pState32f_, pTaps32f_,
                                                nsections_,
                                                dlySrc32f_, pBuf_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Failed to initialized biquad filter");
                    clear();
                    return -1;
                }
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
        /// Determines the length of the initial conditions
        int getInitialConditionLength(void) const
        {
            return 2*nsections_;
        }
        /// Gets the number of sections
        int getNumberOfSections(void) const
        {
            return nsections_;
        }
        /// Sets the initial conditions
        int setInitialConditions(const int nz, const double zi[])
        {
            resetInitialConditions();
            int nzRef = getInitialConditionLength();
            if (nz != nzRef){RTSEIS_WARNMSG("%s", "Shouldn't be here");}
            ippsCopy_64f(zi, zi_, nzRef);
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(zi_, dlySrc64f_, nzRef);
            }
            else
            {
                ippsConvert_64f32f(zi_, dlySrc32f_, nzRef);
            }
            return 0;
        }
        /// Resets the initial conditions
        int resetInitialConditions(void)
        {
            if (precision_ == RTSeis::Precision::DOUBLE)
            {
                ippsCopy_64f(zi_, dlySrc64f_, 2*nsections_);
            }
            else
            {
                ippsConvert_64f32f(zi_, dlySrc32f_, 2*nsections_);
            }
            return 0; 
        } 
        /// Applies the filter
        int apply(const int n, const double x[], double y[])
        {
            if (n <= 0){return 0;}
            if (precision_ == RTSeis::Precision::FLOAT)
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
            // Get a pointer to the filter state and set the initial conditions
            IppStatus status = ippsIIRSetDlyLine_64f(pState64f_, dlySrc64f_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to set delay line");
                return -1;
            }
            // Apply the filters
            status = ippsIIR_64f(x, y, n, pState64f_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to apply filter");
                return -1;
            }
            if (mode_ == RTSeis::ProcessingMode::REAL_TIME)
            {
                status = ippsIIRGetDlyLine_64f(pState64f_, dlyDst64f_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Failed to get delay line");
                    return -1;
                }
                ippsCopy_64f(dlyDst64f_, dlySrc64f_, 2*nsections_);
            }
            return 0;
        }
        /// Applies the filter
        int apply(const int n, const float x[], float y[])
        {
            if (n <= 0){return 0;} 
            if (precision_ == RTSeis::Precision::DOUBLE)
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
            // Get a pointer to the filter state and set the initial conditions
            IppStatus status = ippsIIRSetDlyLine_32f(pState32f_, dlySrc32f_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to set delay line");
                return -1;
            }
            // Apply the filters
            status = ippsIIR_32f(x, y, n, pState32f_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to apply filter");
                return -1;
            }
            if (mode_ == RTSeis::ProcessingMode::REAL_TIME)
            {
                status = ippsIIRGetDlyLine_32f(pState32f_, dlyDst32f_);
                if (status != ippStsNoErr)
                {
                    RTSEIS_ERRMSG("%s", "Failed to get delay line");
                    return -1;
                }
                ippsCopy_32f(dlyDst32f_, dlySrc32f_, 2*nsections_);
            }
            return 0;
        }
    private:
        IppsIIRState_64f *pState64f_ = nullptr;
        /// Filter taps.  This has dimension [tapsLen_].
        Ipp64f *pTaps64f_ = nullptr;
        /// Initial conditions. This has dimension [nwork_].
        Ipp64f *dlySrc64f_ = nullptr;
        /// Final conditions.  This has dimension [nwork_].
        Ipp64f *dlyDst64f_ = nullptr;
        /// Final conditions.  This has dimension [nwork_].
        /// Handle on filter state. 
        IppsIIRState_32f *pState32f_ = nullptr;
        /// Filter taps.  This has dimension [tapsLen_].
        Ipp32f *pTaps32f_ = nullptr;
        /// Initial conditions. This has dimension [nwork_].
        Ipp32f *dlySrc32f_ = nullptr;
        /// Final conditions.  This has dimension [nwork_].
        Ipp32f *dlyDst32f_ = nullptr;
        /// The workspace buffer.
        Ipp8u *pBuf_ = nullptr;
        /// A copy of the numerator filter coefficients.  This has
        /// dimension [3 x nsections_].
        double *bsRef_ = nullptr;
        /// A copy of the denominator filter coefficients.  This has
        /// dimension [3 x nsections_].
        double *asRef_ = nullptr; 
        /// A copy of the initial conditions.  This has dimension
        /// [2 x nsections_].
        double *zi_ = nullptr;
        /// The number of sections.
        int nsections_ = 0;
        /// The number of filter taps.  This equals 6*nsections_.
        int tapsLen_ = 0;
        /// Workspace for the delay lines.
        int nwork_ = 0;
        /// Size of workspace buffer.
        int bufferSize_ = 0;
        /// By default the module does post-procesing.
        RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        /// The default module implementation.
        RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
        /// Flag indicating the module is intiialized
        bool linit_ = false;
};

//============================================================================//

SOSFilter::SOSFilter(void) :
    pSOS_(new SOSFilterImpl())
{
    return;
}

SOSFilter::~SOSFilter(void)
{
    clear();
    return;
}

void SOSFilter::clear(void)
{
    pSOS_->clear();
    return;
}

SOSFilter::SOSFilter(const SOSFilter &sos)
{
    *this = sos;
    return;
}

SOSFilter& SOSFilter::operator=(const SOSFilter &sos)
{
    if (&sos == this){return *this;}
    if (pSOS_){pSOS_->clear();}
    pSOS_ = std::unique_ptr<SOSFilterImpl> (new SOSFilterImpl(*sos.pSOS_));
    return *this;
}

int SOSFilter::initialize(const int ns,
                          const double bs[],
                          const double as[],
                          const RTSeis::ProcessingMode mode,
                          const RTSeis::Precision precision)
{
    clear();
    // Checks
    if (ns < 1 || bs == NULL || as == NULL)
    {
        if (ns < 1){RTSEIS_ERRMSG("%s", "No sections\n");}
        if (bs == NULL){RTSEIS_ERRMSG("%s", "bs is NULL\n");}
        if (as == NULL){RTSEIS_ERRMSG("%s", "as is NULL\n");}
        return -1;
    }
    // Verify the highest order coefficients make sense
    for (int i=0; i<ns; i++)
    {
        if (bs[3*i] == 0.0)
        {
            RTSEIS_ERRMSG("Leading bs coefficient of section %d is zero", i);
            return -1;
        }
        if (as[3*i] == 0.0)
        {
            RTSEIS_ERRMSG("Leading as coefficient of section %d is zero", i);
            return -1;
        }
    }
    int ierr = pSOS_->initialize(ns, bs, as, mode, precision);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize sos filter");
        clear();
        return -1;
    }
    return 0;
}

int SOSFilter::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    resetInitialConditions();
    int nzRef = pSOS_->getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef){RTSEIS_ERRMSG("nz=%d should equal %d", nz, nzRef);}
        if (zi == nullptr){RTSEIS_ERRMSG("%s", "zi is NULL");}
        return -1;
    }
    pSOS_->setInitialConditions(nz, zi);
    return 0;
}

int SOSFilter::resetInitialConditions(void)
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    pSOS_->resetInitialConditions();
    return 0;
}

int SOSFilter::apply(const int n, const double x[], double y[]) 
{
    if (n <= 0){return 0;}
    if (!isInitialized())
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
    int ierr = pSOS_->apply(n, x, y); 
    if (ierr != 0)
    {   
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1; 
    }   
    return 0;
}

int SOSFilter::apply(const int n, const float x[], float y[])
{
    if (n <= 0){return 0;}
    if (!isInitialized())
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
    int ierr = pSOS_->apply(n, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return -1;
    }
    return 0;
}

int SOSFilter::getInitialConditionLength(void) const
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    return pSOS_->getInitialConditionLength();
}

int SOSFilter::getNumberOfSections(void) const
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Class not initialized");
        return -1;
    }
    return pSOS_->getNumberOfSections();
}

bool SOSFilter::isInitialized(void) const
{
    return pSOS_->isInitialized();
}
