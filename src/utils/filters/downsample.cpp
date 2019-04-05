#include <cstdio>
#include <cstdlib>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#define IPPS_CORE_SRC 1
#include "rtseis/utilities/filterImplementations/downsample.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterImplementations;

class Downsample::DownsampleImpl
{
public:
    /// Default constructor
    DownsampleImpl(void)
    {
        return;
    }
    /// Copy constructor
    DownsampleImpl(const DownsampleImpl &downsample)
    {
        *this = downsample;
        return;
    }
    /// (Deep) copy operator
    DownsampleImpl& operator=(const DownsampleImpl &downsample)
    {
        if (&downsample == this){return *this;}
        phase0_ = downsample.phase0_;
        downFactor_ = downsample.downFactor_;
        phase_ = downsample.phase_;
        mode_ = downsample.mode_;
        precision_ = downsample.precision_;
        linit_ = downsample.linit_; 
        return *this;
    }
    /// Default destructor
    ~DownsampleImpl(void) = default;
    /// Clears memory and resets the module
    void clear(void)
    {
        phase0_ = 0;
        downFactor_ = 0;
        phase_ = 0;
        mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        precision_ = RTSeis::Precision::DOUBLE;
        linit_ = false;
        return;
    }
    //--------------------------------------------------------------------//
    /// Initialize the filter
    int initialize(const int downFactor,
                   const RTSeis::ProcessingMode mode,
                   const RTSeis::Precision precision)
    {
        clear();
        downFactor_ = downFactor;
        phase0_ = 0;
        phase_ = 0;
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
    /// Estimates the space required to hold the downsampled signal
    int estimateSpace(const int n) const
    {
        int phase = 0;
        if (mode_ ==  RTSeis::ProcessingMode::REAL_TIME){phase = phase_;}
        int pDstLen = (n + downFactor_ - 1 - phase)/downFactor_;
        return pDstLen;  
    }
    /// Gets the downsampling factor
    int getDownsampleFactor(void) const
    {
        return downFactor_;
    }
    /// Sets the initial conditions
    int setInitialConditions(const int phase)
    {
        resetInitialConditions();
        if (phase < 0 || phase > downFactor_ - 1) 
        {
            return -1;
        }
        phase0_ = phase;
        phase_  = phase;
        return 0;
    }
    /// Resets the initial conditions
    int resetInitialConditions(void)
    {
        phase_ = phase0_;
        return 0;
    }
    /// Apply the downsampler
    int apply(const int nx, const double x[],
              int *nyDown, double y[])
    {
        *nyDown = 0;
        if (nx == 0){return 0;} 
        // There's really nothing to do
        if (downFactor_ == 1)
        {
            ippsCopy_64f(x, y, nx);
            *nyDown = nx; 
            return 0;
        }
        // Apply downsampler
        int pDstLen = estimateSpace(nx);
        int phase = 0;
        int pEst = pDstLen;
        if (mode_ == RTSeis::ProcessingMode::REAL_TIME){phase = phase_;}
        IppStatus status = ippsSampleDown_64f(x, nx, y, &pEst,
                                              downFactor_, &phase);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to downsample signal");
            return -1; 
        }
        // TODO - this is hack'ish.  I think I need to try cacheing vectors
        // that are too small to downsample then hit them on the next
        // iteration. 
        // Weird unhandled exception
        if (pEst < 0 || pEst >= nx){pEst = pDstLen;}
        *nyDown = pEst;
        if (mode_ == RTSeis::ProcessingMode::REAL_TIME){phase_ = phase;}
        return 0;
    }
    /// Apply the downsampler
    int apply(const int nx, const float x[],
              int *nyDown, float y[])
    {
        *nyDown = 0;
        if (nx == 0){return 0;}
        // There's really nothing to do
        if (downFactor_ == 1)
        {
            ippsCopy_32f(x, y, nx);
            *nyDown = nx;
            return 0;
        }
        // Apply downsampler
        int pDstLen = estimateSpace(nx);
        int phase = 0;
        int pEst = pDstLen;
        if (mode_ == RTSeis::ProcessingMode::REAL_TIME){phase = phase_;}
        IppStatus status = ippsSampleDown_32f(x, nx, y, &pEst,
                                              downFactor_, &phase);
        if (status != ippStsNoErr)
        {
            RTSEIS_ERRMSG("%s", "Failed to downsample signal");
           return -1;
        }
        // Weird unhandled exception
        if (pEst < 0 || pEst >= nx){pEst = pDstLen;}
        *nyDown = pEst;
        if (mode_ == RTSeis::ProcessingMode::REAL_TIME){phase_ = phase;}
        return 0;
    }
private:
    /// Initial conditions for phase.
    int phase0_ = 0;
    /// Downsampling factor.
    int downFactor_ = 0;
    /// The phase.
    int phase_ = 0;
    /// By default the module does post-procesing.
    RTSeis::ProcessingMode mode_ = RTSeis::ProcessingMode::POST_PROCESSING;
    /// The default module implementation.
    RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
    /// Flag indicating the module is initialized.
    bool linit_ = false; 
};

//============================================================================//

Downsample::Downsample(void) :
    pDownsample_(std::make_unique<DownsampleImpl>())
{
    return;
}

Downsample::~Downsample(void)
{
    clear();
    return;
}

Downsample::Downsample(const Downsample &downsample)
{
    *this = downsample;
    return;
}

Downsample::Downsample(Downsample &&downsample)
{
    *this = std::move(downsample); 
    return;
}

Downsample& Downsample::operator=(const Downsample &downsample)
{
    if (&downsample == this){return *this;}
    pDownsample_ = std::make_unique<DownsampleImpl> (*downsample.pDownsample_);
    //pDownsample_ = std::unique_ptr<DownsampleImpl>
    //               (new DownsampleImpl(*downsample.pDownsample_));
    return *this;
}

Downsample& Downsample::operator=(Downsample &&downsample)
{
    if (&downsample == this){return *this;}
    pDownsample_ = std::move(downsample.pDownsample_);
    return *this; 
}

int Downsample::initialize(const int downFactor,
                           const RTSeis::ProcessingMode mode,
                           const RTSeis::Precision precision)
{
    clear();
    if (downFactor < 1)
    {
        RTSEIS_ERRMSG("Downsampling factor=%d must be positive", downFactor);
        return -1;
    }
    pDownsample_->initialize(downFactor, mode, precision);
    return 0;
}

void Downsample::clear(void)
{
    if (pDownsample_){pDownsample_->clear();}
    return;
}

int Downsample::setInitialConditions(const int phase)
{
    if (!isInitialized())
    {   
        RTSEIS_ERRMSG("%s", "Downsampler not initialized");
        return -1;
    }
    int downFactor = pDownsample_->getDownsampleFactor();
    if (phase < 0 || phase > downFactor - 1)
    {
        RTSEIS_ERRMSG("phase=%d must be in range[0,%d]", phase, downFactor-1);
        return -1; 
    }
    pDownsample_->setInitialConditions(phase);
    return 0;
}

int Downsample::resetInitialConditions(void)
{
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Downsampler not initialized");
        return -1;
    }
    pDownsample_->resetInitialConditions();
    return 0;
}

int Downsample::apply(const int nx, const double x[], const int ny,
                      int *nyDown, double y[])
{
    *nyDown = 0;
    if (nx <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "ippsDS structure not intitialized");
        return -1;
    }
    int pDstLen = pDownsample_->estimateSpace(nx); 
    if (ny < pDstLen || pDstLen < 0)
    {
        if (pDstLen < 1)
        {
            RTSEIS_ERRMSG("%s", "Space estimate error");
            return -1;
        }
        RTSEIS_ERRMSG("ny=%d must be at least length=%d", ny, pDstLen);
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1;
    }
    int ierr = pDownsample_->apply(nx, x, nyDown, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply downsampler");
        return -1;
    }
    return 0;
}

int Downsample::apply(const int nx, const float x[], const int ny, 
                          int *nyDown, float y[])
{
    *nyDown = 0;
    if (nx <= 0){return 0;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_ERRMSG("%s", "ippsDS structure not intitialized");
        return -1;
    }
    int pDstLen = estimateSpace(nx); 
    if (ny < pDstLen || pDstLen < 0)
    {
        if (pDstLen < 1)
        {
            RTSEIS_ERRMSG("%s", "Space estimate error");
            return -1; 
        }
        RTSEIS_ERRMSG("ny=%d must be at least length=%d", ny, pDstLen);
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "Error x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "Error y is NULL");}
        return -1; 
    }
    int ierr = pDownsample_->apply(nx, x, nyDown, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply downsampler");
        return -1;
    }
    return 0;
}

bool Downsample::isInitialized(void) const
{
    return pDownsample_->isInitialized();
}

int Downsample::estimateSpace(const int n) const
{
    if (!isInitialized() || n <= 0)
    {
        if (!isInitialized()){RTSEIS_ERRMSG("%s", "Class not initialized");}
        if (n < 0){RTSEIS_ERRMSG("n=%d cannot be negative", n);}
        return -1;
    }
    return pDownsample_->estimateSpace(n);
}
 
int Downsample::getDownsampleFactor(void) const
{
    return pDownsample_->getDownsampleFactor();
}
