#include <cstdio>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/enums.hpp"
#include "private/throw.hpp"
#include "rtseis/log.h"
#include "rtseis/utilities/filterImplementations/multiRateFIRFilter.hpp"

using namespace RTSeis::Utilities::FilterImplementations;

template<class T>
class MultiRateFIRFilter<T>::MultiRateFIRImpl
{
public:
    /// Default constructor.
    MultiRateFIRImpl() = default;
    /// Copy constructor.
    MultiRateFIRImpl(const MultiRateFIRImpl &firmr)
    {
        *this = firmr;
    }
    /// (Deep) copy operator.
    MultiRateFIRImpl& operator=(const MultiRateFIRImpl &firmr)
    {
        if (&firmr == this){return *this;}
        if (!firmr.linit_){return *this;}
        int ierr = initialize(firmr.upFactor_, firmr.downFactor_,
                              firmr.tapsLen_, firmr.tapsRef_,
                              firmr.mode_,
                              firmr.precision_,
                              firmr.chunkSize_);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to inititalize");
            clear();
            return *this;
        }
        ippsCopy_64f(firmr.zi_, zi_, mbDly_);
        if (firmr.precision_ == RTSeis::Precision::DOUBLE)
        {
            ippsCopy_64f(firmr.work64_, work64_, downFactor_);
            ippsCopy_64f(firmr.pSrc64_, pSrc64_, nwork_);
            ippsCopy_64f(firmr.pDst64_, pDst64_, nwork_);
            ippsCopy_64f(firmr.pDlySrc64_, pDlySrc64_, mbDly_);
            ippsCopy_64f(firmr.pDlyDst64_, pDlyDst64_, mbDly_);
            if (specSize_ > 0)
            {
                Ipp8u *pIn = reinterpret_cast<Ipp8u *> (firmr.pSpec64_);
                Ipp8u *pOut = reinterpret_cast<Ipp8u *> (pSpec64_);
                ippsCopy_8u(pIn, pOut, specSize_);
            }
        }
        else
        {
            ippsCopy_32f(firmr.work32_, work32_, downFactor_);
            ippsCopy_32f(firmr.pSrc32_, pSrc32_, nwork_);
            ippsCopy_32f(firmr.pDst32_, pDst32_, nwork_);
            ippsCopy_32f(firmr.pDlySrc32_, pDlySrc32_, mbDly_);
            ippsCopy_32f(firmr.pDlyDst32_, pDlyDst32_, mbDly_);
            if (specSize_ > 0)
            {
                Ipp8u *pIn = reinterpret_cast<Ipp8u *> (firmr.pSpec32_);
                Ipp8u *pOut = reinterpret_cast<Ipp8u *> (pSpec32_);
                ippsCopy_8u(pIn, pOut, specSize_);
            }
        }
        if (bufferSize_ > 0){ippsCopy_8u(firmr.pBuf_, pBuf_, bufferSize_);}
        chunkSize_ = firmr.chunkSize_;
        nExcess_ = firmr.nExcess_;
        nbDly_ = firmr.nbDly_;
        mbDly_ = firmr.mbDly_;
        upPhase_ = firmr.upPhase_;
        downPhase_ = firmr.downPhase_;
        upFactor_ = firmr.upFactor_;
        downPhaseIC_ = firmr.downPhaseIC_;
        upPhaseIC_ = firmr.upPhaseIC_;
        downFactor_ = firmr.downFactor_;
        tapsLen_ = firmr.tapsLen_;
        order_ = firmr.order_;
        nwork_ = firmr.nwork_;
        specSize_ = firmr.specSize_;
        bufferSize_ = firmr.bufferSize_;
        mode_ = firmr.mode_;
        precision_ = firmr.precision_;
        linit_ = firmr.linit_; 
        return *this;
    }
    /// Destructor.
    ~MultiRateFIRImpl()
    {
        clear();
    }
    //========================================================================//
    /// Clears the filter/releases the memory.
    void clear() noexcept
    {
        if (pSpec64_   != nullptr){ippsFree(pSpec64_);}
        if (pTaps64_   != nullptr){ippsFree(pTaps64_);}
        if (work64_    != nullptr){ippsFree(work64_);}
        if (pSrc64_    != nullptr){ippsFree(pSrc64_);}
        if (pDst64_    != nullptr){ippsFree(pDst64_);}
        if (pDlySrc64_ != nullptr){ippsFree(pDlySrc64_);}
        if (pDlyDst64_ != nullptr){ippsFree(pDlyDst64_);} 
        if (pSpec32_   != nullptr){ippsFree(pSpec32_);}
        if (pTaps32_   != nullptr){ippsFree(pTaps32_);}
        if (work32_    != nullptr){ippsFree(work32_);}
        if (pSrc32_    != nullptr){ippsFree(pSrc32_);}
        if (pDst32_    != nullptr){ippsFree(pDst32_);}
        if (pDlySrc32_ != nullptr){ippsFree(pDlySrc32_);}
        if (pDlyDst32_ != nullptr){ippsFree(pDlyDst32_);}
        if (pBuf_      != nullptr){ippsFree(pBuf_);}
        if (tapsRef_   != nullptr){ippsFree(tapsRef_);}
        if (zi_        != nullptr){ippsFree(zi_);}
        pSpec64_ = nullptr;
        pTaps64_ = nullptr;
        work64_ = nullptr;
        pSrc64_ = nullptr;
        pDst64_ = nullptr;
        pDlySrc64_ = nullptr;
        pDlyDst64_ = nullptr;
        pSpec32_ = nullptr;
        pTaps32_ = nullptr;
        work32_ = nullptr;
        pSrc32_ = nullptr;
        pDst32_ = nullptr;
        pDlySrc32_ = nullptr;
        pDlyDst32_ = nullptr;
        pBuf_ = nullptr;
        tapsRef_ = nullptr; 
        zi_ = nullptr;
        chunkSize_ = defaultChunkSize_;
        nExcess_ = 0;
        nbDly_ = 0;
        mbDly_ = 0;
        upPhase_ = 0;
        downPhase_ = 0;
        upFactor_ = 0;
        downPhaseIC_ = 0;
        upPhaseIC_ = 0;
        downFactor_ = 0;  
        tapsLen_ = 0;
        order_ = 0;
        nwork_ = 0;
        specSize_ = 0;
        bufferSize_ = 0;
        mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        precision_ = RTSeis::Precision::DOUBLE;
        linit_ = false;
    }
    //========================================================================//
    /// Gets the length of the intiial conditions
    [[nodiscard]] int getInitialConditionLength() const
    {
        return nbDly_;
    }
    /// Sets the initial conditions
    int setInitialConditions(const int nz, const double zi[])
    {
        resetInitialConditions();
        int nzRef = getInitialConditionLength();
        if (nzRef != nz){RTSEIS_WARNMSG("%s", "Shouldn't happen");}
        if (nzRef > 0){ippsCopy_64f(zi, zi_, nzRef);}
        ippsCopy_64f(zi, zi_, nzRef);
        if (precision_ == RTSeis::Precision::DOUBLE)
        {
            ippsCopy_64f(zi_, pDlySrc64_, nzRef); 
        }
        else
        {
            ippsConvert_64f32f(zi_, pDlySrc32_, nzRef);
        }
        return 0;
    }
    /// Resets the initial conditions
    int resetInitialConditions()
    {
        upPhase_ = 0;
        downPhase_ = 0;
        nExcess_ = 0; 
        if (precision_ == RTSeis::Precision::DOUBLE)
        {
            ippsZero_64f(work64_, downFactor_);
            ippsZero_64f(pDlySrc64_, mbDly_);
            if (nbDly_ > 0)
            {
                ippsCopy_64f(zi_, pDlySrc64_, nbDly_);
            }
            IppStatus status = ippsFIRMRInit_64f(pTaps64_, tapsLen_,
                                                 upFactor_, upPhase_,
                                                 downFactor_, downPhase_,
                                                 pSpec64_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error reinitializing state structure");
                return -1;
            }
        }
        else
        {
            ippsZero_32f(work32_, downFactor_);
            ippsZero_32f(pDlySrc32_, mbDly_);
            if (nbDly_ > 0)
            {   
                ippsConvert_64f32f(zi_, pDlySrc32_, nbDly_);
            }   
            IppStatus status = ippsFIRMRInit_32f(pTaps32_, tapsLen_,
                                                 upFactor_, upPhase_,
                                                 downFactor_, downPhase_,
                                                 pSpec32_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Error reinitializing state structure");
                return -1;
            }
        }
        return 0;
    }
    /// Initializes the filter.
    int initialize(const int upFactor, const int downFactor,
                   const int nb, const double b[],
                   const RTSeis::ProcessingMode mode,
                   const RTSeis::Precision precision,
                   const int chunkSize = 1024)
    {
        clear();
        // Initialize up/down phases are zero - set in initial conditions?
        upPhaseIC_ = 0;
        downPhaseIC_ = 0;
        upPhase_ = 0;
        downPhase_ = 0;
        chunkSize_ = chunkSize;
        // Copy filter basics
        upFactor_ = upFactor;
        downFactor_ = downFactor;
        order_ = nb - 1;
        tapsLen_ = nb;
        // Copy filter
        tapsRef_ = ippsMalloc_64f(tapsLen_);
        ippsCopy_64f(b, tapsRef_, tapsLen_);
        // Copy basics and set sizes
        mode_ = mode;
        precision_ = precision; 
        order_ = nb - 1;
        nbDly_ = (tapsLen_ + upFactor_ - 1)/upFactor_;
        mbDly_ = std::max(128, nbDly_);
        // Reserve space for copy of initial conditions
        zi_ = ippsMalloc_64f(mbDly_);
        ippsZero_64f(zi_, mbDly_);
        // Set workspace size
        nwork_ = chunkSize_*std::max(upFactor_, downFactor_);
        if (precision == RTSeis::Precision::DOUBLE)
        {
            // Get state size
            IppStatus status = ippsFIRMRGetSize(tapsLen_, upFactor_,
                                                downFactor_, ipp64f,
                                                &specSize_, &bufferSize_);
            if (status != ippStsNoErr)
            {   
                RTSEIS_ERRMSG("%s", "Failed to get state size");
                clear();
                return -1; 
            }
            // Workspace
            work64_ = ippsMalloc_64f(downFactor_);
            ippsZero_64f(work64_, downFactor_);
            pSrc64_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(pSrc64_, nwork_);
            pDst64_ = ippsMalloc_64f(nwork_);
            ippsZero_64f(pDst64_, nwork_);
            // Filter state 
            pSpec64_ = reinterpret_cast<IppsFIRSpec_64f *>
                           (ippsMalloc_8u(specSize_));
            // Filter taps
            pTaps64_ = ippsMalloc_64f(tapsLen_);
            if (upFactor_ > 1)
            {
                // By injecting upFactor_ times the points into the signal
                // to filter we have to gain the FIR filter upFactor_ times.
                // This is not done automatically by Matlab.
                auto gain = static_cast<double> (upFactor_);
                ippsMulC_64f(b, gain, pTaps64_, tapsLen_); 
            }
            else
            {
                ippsCopy_64f(b, pTaps64_, tapsLen_);
            }
            // Delay lines
            pDlySrc64_ = ippsMalloc_64f(mbDly_);
            ippsZero_64f(pDlySrc64_, mbDly_);
            pDlyDst64_ = ippsMalloc_64f(mbDly_);
            ippsZero_64f(pDlyDst64_, mbDly_);
            // Workspace
            pBuf_ = ippsMalloc_8u(bufferSize_);
            ippsZero_8u(pBuf_, bufferSize_);
            // Initialize
            status = ippsFIRMRInit_64f(pTaps64_, tapsLen_,
                                       upFactor_, upPhase_,
                                       downFactor_, downPhase_,
                                       pSpec64_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to initialize filter state");
                clear();
                return -1;
            }
        }
        else
        {
            // Get state size
            IppStatus status = ippsFIRMRGetSize(tapsLen_, upFactor_,
                                                downFactor_, ipp32f,
                                                &specSize_, &bufferSize_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to get state size");
                clear();
                return -1; 
            }
            // Workspace
            work32_ = ippsMalloc_32f(downFactor_);
            ippsZero_32f(work32_, downFactor_);
            pSrc32_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(pSrc32_, nwork_);
            pDst32_ = ippsMalloc_32f(nwork_);
            ippsZero_32f(pDst32_, nwork_);
            // Filter state 
            pSpec32_ = reinterpret_cast<IppsFIRSpec_32f *>
                       (ippsMalloc_8u(specSize_));
            // Filter taps
            pTaps32_ = ippsMalloc_32f(tapsLen_);
            ippsConvert_64f32f(b, pTaps32_, tapsLen_);
            if (upFactor_ > 1)
            {
                // By injecting upFactor_ times the points into the signal
                // to filter we have to gain the FIR filter upFactor_ times.
                // This is not done automatically by Matlab.
                auto gain = static_cast<float> (upFactor_);
                ippsMulC_32f_I(gain, pTaps32_, tapsLen_); 
            }
            // Delay lines
            pDlySrc32_ = ippsMalloc_32f(mbDly_);
            ippsZero_32f(pDlySrc32_, mbDly_);
            pDlyDst32_ = ippsMalloc_32f(mbDly_);
            ippsZero_32f(pDlyDst32_, mbDly_);
            // Workspace
            pBuf_ = ippsMalloc_8u(bufferSize_);
            ippsZero_8u(pBuf_, bufferSize_);
            // Initialize
            status = ippsFIRMRInit_32f(pTaps32_, tapsLen_,
                                       upFactor_, upPhase_,
                                       downFactor_, downPhase_,
                                       pSpec32_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("%s", "Failed to initialize filter state");
                clear();
                return -1;
            }
        }
        mode_ = mode;
        precision_ = precision;
        linit_ = true;
        return 0;
    }

    /// Applies the filter
    int apply(const int n, const double x[], const int ny,
              int *len, double y[])
    {
        if (n <= 0){return 0;} // Nothing to do
        if (precision_ == RTSeis::Precision::FLOAT)
        {
            Ipp32f *x32 = ippsMalloc_32f(n);
            Ipp32f *y32 = ippsMalloc_32f(ny);
            ippsConvert_64f32f(x, x32, n);
            int ierr = apply(n, x32, ny, len, y32);
            ippsFree(x32);
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Failed to apply filter");
                ippsFree(y32);
                return -1;
            }
            ippsConvert_32f64f(y32, y, *len);
            ippsFree(y32);
            return 0;
        }
        // Compute output length
        int upFactor   = upFactor_;
        int downPhase  = downPhase_;
        int downFactor = downFactor_;
        int downPhaseNew = (downFactor + downPhase - n%downFactor)%downFactor;
        *len = (upFactor*n + downFactor - 1 - downPhase)/downFactor;
        // Get pointers
        Ipp64f *pSrc = nullptr;
        Ipp64f *pDst = nullptr;
        bool lfree = false;
        if (n + downFactor <= chunkSize_)
        {
            pSrc = pSrc64_;
            pDst = pDst64_;
        }   
        else
        {
            pSrc = ippsMalloc_64f((n + downFactor)*downFactor_);
            pDst = ippsMalloc_64f((n + downFactor)*upFactor_);
            lfree = true;
        }
        ippsZero_64f(pSrc, (n + downFactor)*downFactor_);
        ippsZero_64f(pDst, (n + downFactor)*upFactor_);
        Ipp64f *dlysrc = pDlySrc64_;
        Ipp64f *dlydst = pDlyDst64_;
        int nUse = n;
        int nprev = 0;
        int nexcess = 0;
        if (mode_ == RTSeis::ProcessingMode::POST_PROCESSING)
        {
            downPhaseNew = 0; // Downsampling phase shouldn't change
            nUse = n;
            ippsZero_64f(dlysrc, mbDly_);
            *len = (upFactor*n + downFactor - 1 - downPhase)/downFactor;
            ippsCopy_64f(x, pSrc, nUse);
        }
        else
        {
            // Get the straggler elements from the previous application
            nprev = nExcess_;
            if (nprev > 0){ippsCopy_64f(work64_, pSrc, nprev);}
            int nNew = nprev + n; // Number of new samples
            // It's possible that there's nothing to do
            if (nNew/downFactor < 1)
            {
                *len = 0;
                ippsCopy_64f(x, &pSrc[nprev], n);
                ippsCopy_64f(pSrc, work64_, nprev + n);
                nExcess_ = nprev + n;
                if (lfree)
                {
                    ippsFree(pSrc);
                    ippsFree(pDst);
                }
                return 0;
            }
            nUse = nNew - nNew%downFactor; // Number of samples to process
            nexcess = (nprev + n) - nUse; // Number of excess samples
            int nCopy = n - nexcess; // Number of samples to copy
            nexcess = nNew - nUse;
            if (nCopy > 0){ippsCopy_64f(x, &pSrc[nprev], nCopy);}
            if (nUse%downFactor != 0)
            {
                RTSEIS_ERRMSG("mod(nUse,downFactor) = mod(%d,%d) = %d != 0",
                              nUse, downFactor, nUse%downFactor);
            }
        }
        *len = 0;
        if (nUse/downFactor > 0)
        {
            *len = (upFactor*nUse + downFactor - 1 - downPhase)/downFactor;
            if (*len > ny)
            {
                RTSEIS_ERRMSG("ny=%d must be at least %d\n", ny, *len);
                return -2;
            }
            // Apply it
            IppStatus status = ippsFIRMR_64f(pSrc, pDst,
                                             nUse/downFactor, pSpec64_,
                                             dlysrc, dlydst, pBuf_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("Error in FIRMR with %d/%d=%d samples",
                              nUse, downFactor, nUse/downFactor);
                return -1;
            }
            ippsCopy_64f(pDst, y, *len);
            if (mode_ == RTSeis::ProcessingMode::REAL_TIME)
            {
                ippsCopy_64f(dlydst, dlysrc, nbDly_);
                ippsZero_64f(work64_, downFactor);
                if (nexcess > 0)
                {
                    ippsCopy_64f(&x[n-nexcess], work64_, nexcess);
                }
                nExcess_ = nexcess;
            }
        }
        // Nothing is coming back
        else
        {
            *len = 0;
            if (mode_ == RTSeis::ProcessingMode::REAL_TIME)
            {
                RTSEIS_WARNMSG("%s", "Shouldn't be here\n");
                ippsCopy_64f(pSrc, work64_, nprev);
                ippsCopy_64f(&x[n-nexcess], &work64_[nprev], nexcess);
                nExcess_ = nexcess;
            }
        }
        if (lfree)
        {
            ippsFree(pSrc);
            ippsFree(pDst);
        }
        downPhase_ = downPhaseNew;
        return 0;
    }
 
    /// Applies the filter
    int apply(const int n, const float x[], const int ny,
              int *len, float y[])
    {
        if (n <= 0){return 0;} // Nothing to do
        if (precision_ == RTSeis::Precision::DOUBLE)
        {
            Ipp64f *x64 = ippsMalloc_64f(n);
            Ipp64f *y64 = ippsMalloc_64f(ny);
            ippsConvert_32f64f(x, x64, n); 
            int ierr = apply(n, x64, ny, len, y64);
            ippsFree(x64);
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("%s", "Failed to apply filter");
                ippsFree(y64);
                return -1; 
            }
            ippsConvert_64f32f(y64, y, *len);
            ippsFree(y64);
            return 0;
        }
        // Compute output length
        int upFactor   = upFactor_;
        int downPhase  = downPhase_;
        int downFactor = downFactor_;
        int downPhaseNew = (downFactor + downPhase - n%downFactor)%downFactor;
        *len = (upFactor*n + downFactor - 1 - downPhase)/downFactor;
        // Get pointers
        Ipp32f *pSrc = nullptr;
        Ipp32f *pDst = nullptr;
        bool lfree = false;
        if (n + downFactor <= chunkSize_)
        {
            pSrc = pSrc32_;
            pDst = pDst32_;
        }   
        else
        {
            pSrc = ippsMalloc_32f((n + downFactor)*downFactor_);
            pDst = ippsMalloc_32f((n + downFactor)*upFactor_);
            lfree = true;
        }
        ippsZero_32f(pSrc, (n + downFactor)*downFactor_);
        ippsZero_32f(pDst, (n + downFactor)*upFactor_);
        Ipp32f *dlysrc = pDlySrc32_;
        Ipp32f *dlydst = pDlyDst32_;
        int nUse = n;
        int nprev = 0;
        int nexcess = 0;
        if (mode_ == RTSeis::ProcessingMode::POST_PROCESSING)
        {
            nUse = n;
            ippsZero_32f(dlysrc, mbDly_);
            *len = (upFactor*n + downFactor - 1 - downPhase)/downFactor;
            ippsCopy_32f(x, pSrc, nUse);
        }
        else
        {
            // Get the straggler elements from the previous application
            nprev = nExcess_;
            if (nprev > 0){ippsCopy_32f(work32_, pSrc, nprev);}
            int nNew = nprev + n; // Number of new samples
            // It's possible that there's nothing to do
            if (nNew/downFactor < 1)
            {
                *len = 0;
                ippsCopy_32f(x, &pSrc[nprev], n);
                ippsCopy_32f(pSrc, work32_, nprev + n);
                nExcess_ = nprev + n;
                if (lfree)
                {
                    ippsFree(pSrc);
                    ippsFree(pDst);
                }
                return 0;
            }
            nUse = nNew - nNew%downFactor; // Number of samples to process
            nexcess = (nprev + n) - nUse; // Number of excess samples
            int nCopy = n - nexcess; // Number of samples to copy
            nexcess = nNew - nUse;
            if (nCopy > 0){ippsCopy_32f(x, &pSrc[nprev], nCopy);}
            if (nUse%downFactor != 0)
            {
                RTSEIS_ERRMSG("mod(nUse,downFactor) = mod(%d,%d) = %d != 0",
                              nUse, downFactor, nUse%downFactor);
            }
        }
        *len = 0;
        if (nUse/downFactor > 0)
        {
            *len = (upFactor*nUse + downFactor - 1 - downPhase)/downFactor;
            if (*len > ny)
            {
                RTSEIS_ERRMSG("ny=%d must be at least %d\n", ny, *len);
                return -2;
            }
            // Apply it
            IppStatus status = ippsFIRMR_32f(pSrc, pDst,
                                             nUse/downFactor, pSpec32_,
                                             dlysrc, dlydst, pBuf_);
            if (status != ippStsNoErr)
            {
                RTSEIS_ERRMSG("Error in FIRMR with %d/%d=%d samples",
                              nUse, downFactor, nUse/downFactor);
                return -1;
            }
            ippsCopy_32f(pDst, y, *len);
            if (mode_ == RTSeis::ProcessingMode::REAL_TIME)
            {
                ippsCopy_32f(dlydst, dlysrc, nbDly_);
                ippsZero_32f(work32_, downFactor);
                if (nexcess > 0)
                {
                    ippsCopy_32f(&x[n-nexcess], work32_, nexcess);
                }
                nExcess_ = nexcess;
            }
        }
        // Nothing is coming back
        else
        {
            *len = 0;
            if (mode_ == RTSeis::ProcessingMode::REAL_TIME)
            {
                RTSEIS_WARNMSG("%s", "Shouldn't be here\n");
                ippsCopy_32f(pSrc, work32_, nprev);
                ippsCopy_32f(&x[n-nexcess], &work32_[nprev], nexcess);
                nExcess_ = nexcess;
            }
        }
        if (lfree)
        {
            ippsFree(pSrc);
            ippsFree(pDst);
        }
        downPhase_ = downPhaseNew;
        return 0;
    }

    /// Determines if the filter is initialized.
    bool isInitialized(void) const {return linit_;}          

    /// Estimates space
    int estimateSpace(const int n) const
    {
        int len = (upFactor_*n + downFactor_ - 1 - downPhase_)/downFactor_;
        return len;
    }

private:
    /// Default chunk-size for processing real-time blocks
    const int defaultChunkSize_ = 1024;
    /// The filter state.  This has dimension [specSize_].
    IppsFIRSpec_64f *pSpec64_ = nullptr;
    /// The taps.  This has dimension [tapsLen_].
    Ipp64f *pTaps64_ = nullptr;
    /// Helps with real-time filtering.  This has dimension [downFactor_].
    Ipp64f *work64_ = nullptr; 
    /// Helps with real-time filtering.  This has dimension [nwork_].
    Ipp64f *pSrc64_ = nullptr;
    /// Helps with real-time filtering.  This has dimension [nwork_].
    Ipp64f *pDst64_ = nullptr;
    /// The input delay line.  This has dimension [mbDly_].
    Ipp64f *pDlySrc64_ = nullptr;
    /// The output delay line.  This has dimension [mbDly_].
    Ipp64f *pDlyDst64_ = nullptr; 
    /// The filter state.  This has dimension [specSize_].
    IppsFIRSpec_32f *pSpec32_ = nullptr;
    /// The taps.  This has dimension [tapsLen_].
    Ipp32f *pTaps32_ = nullptr;
    /// Helps with real-time filtering.  This has dimension [downFactor_].
    Ipp32f *work32_ = nullptr;
    /// Helps with real-time filtering.  This has dimension [nwork_].
    Ipp32f *pSrc32_ = nullptr;
    /// Helps with real-time filtering.  This has dimension [nwork_].
    Ipp32f *pDst32_ = nullptr;
    /// The input delay line.  This has dimension [mbDly_].
    Ipp32f *pDlySrc32_ = nullptr;
    /// The output delay line.  This has dimension [mbDly_].
    Ipp32f *pDlyDst32_ = nullptr;
    /// The workspace buffer.  This has dimension [bufferSize_].
    Ipp8u *pBuf_ = nullptr; 
    /// A copy of the initial conditions.  This has dimension [nbDly_].
    double *zi_ = nullptr; 
    /// A copy of the FIR filter coefficients.
    double *tapsRef_ = nullptr;
    /// Length of the chunks for real-time processing.
    int chunkSize_ = defaultChunkSize_;
    /// Left-over signals from previous real-time iteration.
    int nExcess_ = 0; 
    /// Length of the delay lines.
    int nbDly_ = 0;
    /// Max space reserved to delay lines.
    int mbDly_ = 0;
    /// Phase of the upsampler.
    int upPhase_ = 0;
    /// Phase of downsampler.
    int downPhase_ = 0;
    /// Initial conditions for phase of upsampler.
    int upPhaseIC_ = 0;
    /// Initial conditions for phase of downsampler. 
    int downPhaseIC_ = 0;
    /// Upsampling factor.  Upsampling factor - 1 zeros will be stuffed
    /// between samples.
    int upFactor_ = 0;
    /// Downsampling factor.  Downsampling factor - 1 points will be 
    /// skipped while sampling.
    int downFactor_ = 0; 
    /// Number of filter taps.
    int tapsLen_ = 0;
    /// Filter order.
    int order_ = 0;
    /// The length of the workspace array.
    int nwork_ = 0;
    /// The length of the filter state. 
    int specSize_ = 0;
    /// The length of pBuf_.
    int bufferSize_ = 0;
    /// By default the module does post-procesing.
    RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
    /// The default module implementation.
    RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
    /// Flag indicating this is initialized.
    bool linit_ = false;
};

/// C'tor
template<class T>
MultiRateFIRFilter<T>::MultiRateFIRFilter() :
    pFIR_(std::make_unique<MultiRateFIRImpl> ())
{
}

/// Copy c'tor
template<class T>
MultiRateFIRFilter<T>::MultiRateFIRFilter(const MultiRateFIRFilter &firmr)
{
    *this = firmr;
}

/// Move c'tor
template<class T>
MultiRateFIRFilter<T>::MultiRateFIRFilter(MultiRateFIRFilter &&firmr) noexcept
{
    *this = std::move(firmr);
}

/// Copy assignment
template<class T>
MultiRateFIRFilter<T>&
MultiRateFIRFilter<T>::operator=(const MultiRateFIRFilter &firmr)
{
    if (&firmr == this){return *this;}
    if (pFIR_){pFIR_->clear();}
    pFIR_ = std::make_unique<MultiRateFIRImpl> (*firmr.pFIR_);
    return *this;
}

/// Move assignment
template<class T>
MultiRateFIRFilter<T>&
MultiRateFIRFilter<T>::operator=(MultiRateFIRFilter<T> &&firmr) noexcept
{
    if (&firmr == this){return *this;}
    pFIR_ = std::move(firmr.pFIR_);
    return *this;
}

/// Initialize
template<>
void MultiRateFIRFilter<double>::initialize(
    const int upFactor, const int downFactor,
    const int nb, const double b[],
    const int chunkSize,
    const RTSeis::ProcessingMode mode)
{
    constexpr RTSeis::Precision precision = RTSeis::Precision::DOUBLE;
    clear();
    if (upFactor < 1 || downFactor < 1 || nb < 1 ||
        b == nullptr || chunkSize < downFactor)
    {
        if (upFactor < 1)
        {
            RTSEIS_THROW_RTE("Upsampling factor=%d must be positive", upFactor);
        }
        if (downFactor < 1)
        {
            RTSEIS_THROW_RTE("Downsampling factor=%d must be positive",
                             downFactor);
        }
        if (nb < 1){RTSEIS_THROW_RTE("No filter taps; nb=%d", nb);}
        if (downFactor < 1){RTSEIS_THROW_RTE("%s", "b is NULL");}
        if (chunkSize < 1)
        {
            RTSEIS_ERRMSG("Chunksize=%d must be positive", chunkSize);
        }
    }
#ifdef DEBUG
    int ierr = pFIR_->initialize(upFactor, downFactor,
                                 nb, b, mode, precision, chunkSize);
    assert(ierr == 0);
#else
    pFIR_->initialize(upFactor, downFactor,
                      nb, b, mode, precision, chunkSize);
#endif
}

template<>
void MultiRateFIRFilter<float>::initialize(
    const int upFactor, const int downFactor,
    const int nb, const double b[],
    const int chunkSize,
    const RTSeis::ProcessingMode mode)
{
    constexpr RTSeis::Precision precision = RTSeis::Precision::FLOAT;
    clear();
    if (upFactor < 1 || downFactor < 1 || nb < 1 ||
        b == nullptr || chunkSize < downFactor)
    {
        if (upFactor < 1)
        {
            RTSEIS_THROW_RTE("Upsampling factor=%d must be positive", upFactor);
        }
        if (downFactor < 1)
        {
            RTSEIS_THROW_RTE("Downsampling factor=%d must be positive",
                             downFactor);
        }
        if (nb < 1){RTSEIS_THROW_RTE("No filter taps; nb=%d", nb);}
        if (downFactor < 1){RTSEIS_THROW_RTE("%s", "b is NULL");}
        if (chunkSize < 1)
        {
            RTSEIS_ERRMSG("Chunksize=%d must be positive", chunkSize);
        }
    }
#ifdef DEBUG
    int ierr = pFIR_->initialize(upFactor, downFactor,
                                 nb, b, mode, precision, chunkSize);
    assert(ierr == 0);
#else
    pFIR_->initialize(upFactor, downFactor,
                      nb, b, mode, precision, chunkSize);
#endif
}

/// Destructor
template<class T>
MultiRateFIRFilter<T>::~MultiRateFIRFilter() = default;

template<class T>
void MultiRateFIRFilter<T>::clear() noexcept
{
    pFIR_->clear();
}

template<class T>
void MultiRateFIRFilter<T>::initialize(
    const int upFactor, const int downFactor,
    const int nb, const double b[],
    const RTSeis::ProcessingMode mode)
{
    initialize(upFactor, downFactor, nb, b, 1024, mode);
}

template<class T>
int MultiRateFIRFilter<T>::estimateSpace(const int n) const
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Module is not initialized");
    }
    if (n <= 0){return 0;} // No points
    int len = pFIR_->estimateSpace(n);
    return len;
}

template<class T>
int MultiRateFIRFilter<T>::getInitialConditionLength() const
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    return pFIR_->getInitialConditionLength();
}

template<class T>
void MultiRateFIRFilter<T>::setInitialConditions(
    const int nz, const double zi[])
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    int nzRef = getInitialConditionLength();
    if (nz != nzRef || zi == nullptr)
    {
        if (nz != nzRef){RTSEIS_THROW_IA("nz=%d should equal %d", nz, nzRef);}
        RTSEIS_THROW_IA("%s", "zi is NULL");
    }
    pFIR_->setInitialConditions(nz, zi);
}

template<class T>
void MultiRateFIRFilter<T>::apply(const int n, const T x[],
                                  const int nywork, int *ny, T *yIn[])
{
    *ny = 0;
    if (n <= 0){return;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Module is not initialized");
    }
    if (x == nullptr)
    {
        RTSEIS_THROW_IA("%s", "x is NULL");
    }
    int nworkEst = pFIR_->estimateSpace(n);
    if (nywork < nworkEst)
    {
        RTSEIS_THROW_IA("Resize nywork =  %d to %d", nywork, nworkEst);
    }
    T *y = *yIn;
    if (y == nullptr)
    {
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
#ifdef DEBUG
    int ierr = pFIR_->apply(n, x, nywork, ny, y);
    assert(ierr == 0);
#else
    pFIR_->apply(n, x, nywork, ny, y);
#endif
}

/*
int MultiRateFIRFilter::apply(const int n, const float x[],
                              const int nywork, int *ny, float y[])
{
    *ny = 0;
    if (n <= 0){return 0;} // Nothing to do
    if (x == nullptr)
    {
        RTSEIS_ERRMSG("%s", "x is NULL");
        return -1; 
    }
    if (!pFIR_->isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Module is not initialized");
        return -1; 
    }
    int nworkEst = pFIR_->estimateSpace(n);
    if (nywork < nworkEst)
    {
        RTSEIS_WARNMSG("May have insufficient space %d %d", nywork, nworkEst);
    }
    if (y == nullptr)
    {
        RTSEIS_ERRMSG("%s", "y is NULL");
        return -1; 
    }
    int ierr = pFIR_->apply(n, x, nywork, ny, y); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        *ny = 0;
        return -1;
    }
    return 0;
}
*/

/// Reset initial conditions
template<class T>
void MultiRateFIRFilter<T>::resetInitialConditions()
{
    if (!pFIR_->isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Module is not initialized");
    }
#ifdef DEBUG
    int ierr = pFIR_->resetInitialConditions();
    assert(ierr == 0);
#else
    pFIR_->resetInitialConditions();
#endif
}

/// Initialized?
template<class T>
bool MultiRateFIRFilter<T>::isInitialized() const noexcept
{
    return pFIR_->isInitialized();
}

/// Template instantiation
template class RTSeis::Utilities::FilterImplementations::MultiRateFIRFilter<double>;
template class RTSeis::Utilities::FilterImplementations::MultiRateFIRFilter<float>;
