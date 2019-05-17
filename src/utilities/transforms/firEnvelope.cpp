#include <cstdio>
#include <cstdlib>
#include <ipps.h>
#include "rtseis/utilities/transforms/firEnvelope.hpp"
#include "rtseis/utilities/design/fir.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#include "rtseis/private/throw.hpp"

using namespace RTSeis::Utilities;
using namespace RTSeis::Utilities::Transforms;

class FIREnvelope::FIREnvelopeImpl
{
public:
    FIREnvelopeImpl() = default;
    ~FIREnvelopeImpl() = default;
    FIREnvelopeImpl(const FIREnvelopeImpl &firEnvelope)
    {
        *this = firEnvelope;
    }
    FIREnvelopeImpl& operator=(const FIREnvelopeImpl &firEnvelope)
    {
        if (&firEnvelope == this){return *this;}
        mRealFIRFilter = firEnvelope.mRealFIRFilter;
        mImagFIRFilter = firEnvelope.mImagFIRFilter;
        mNumberOfTaps = firEnvelope.mNumberOfTaps;
        mZeroPhase = firEnvelope.mZeroPhase;
        mType3 = firEnvelope.mType3;
        mHaveInitialCondition = firEnvelope.mHaveInitialCondition;
        mInitialized = firEnvelope.mInitialized;
        mMode = firEnvelope.mMode;
        return *this; 
    }

    FilterImplementations::FIRFilter mRealFIRFilter;
    FilterImplementations::FIRFilter mImagFIRFilter;
    int mNumberOfTaps = 0;
    bool mZeroPhase = true;
    bool mType3 = false;
    bool mHaveInitialCondition = false;
    bool mInitialized = false;
    ProcessingMode mMode = ProcessingMode::POST_PROCESSING;
};

/// Constructor
FIREnvelope::FIREnvelope() :
    pImpl(std::make_unique<FIREnvelopeImpl> ())
{
}

/// Copy constructor
FIREnvelope::FIREnvelope(const FIREnvelope &firEnvelope)
{
    *this = firEnvelope;
}

/// Move constructor
FIREnvelope::FIREnvelope(FIREnvelope &&firEnvelope) noexcept
{
    *this = std::move(firEnvelope);
}

/// Copy assignment
FIREnvelope& FIREnvelope::operator=(const FIREnvelope &firEnvelope)
{
    if (&firEnvelope == this){return *this;}
    pImpl = std::make_unique<FIREnvelopeImpl> (*firEnvelope.pImpl);
    return *this; 
}

/// Move assignment
FIREnvelope& FIREnvelope::operator=(FIREnvelope &&firEnvelope) noexcept
{
    if (&firEnvelope == this){return *this;}
    pImpl = std::move(firEnvelope.pImpl);
    return *this;
}

/// Destructor
FIREnvelope::~FIREnvelope() = default;

/// Clear the filter
void FIREnvelope::clear() noexcept
{
    pImpl->mRealFIRFilter.clear();
    pImpl->mImagFIRFilter.clear();
    pImpl->mNumberOfTaps = 0;
    pImpl->mZeroPhase = true;
    pImpl->mType3 = false;
    pImpl->mHaveInitialCondition = false;
    pImpl->mInitialized = false;
    pImpl->mMode = ProcessingMode::POST_PROCESSING;
}

/// Check if initialized
bool FIREnvelope::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Initialize
void FIREnvelope::initialize(const int ntaps,
                             const RTSeis::ProcessingMode mode,
                             const RTSeis::Precision precision)
{
    clear();
    if (ntaps < 1)
    {
        RTSEIS_THROW_IA("ntaps = %d must be positive", ntaps);
    }
    pImpl->mMode = mode;
    pImpl->mType3 = false;
    if (ntaps%2 == 1){pImpl->mType3 = true;}
    pImpl->mNumberOfTaps = ntaps;
    std::pair<FilterRepresentations::FIR, FilterRepresentations::FIR> zfir;
    constexpr double beta = 8;
    constexpr auto direct = FilterImplementations::FIRImplementation::DIRECT;
    // Create an FIR hilbert transform
    try
    {
        auto zfir = FilterDesign::FIR::HilbertTransformer(ntaps - 1, beta);
        std::vector<double> rfir = zfir.first.getFilterTaps();
        pImpl->mRealFIRFilter.initialize(rfir.size(), rfir.data(),
                                        mode, precision,
                                        direct);
        auto cfir = zfir.second.getFilterTaps();
        pImpl->mImagFIRFilter.initialize(cfir.size(), cfir.data(),
                                        mode, precision,
                                        direct);
    }
    catch (const std::exception &e) 
    {
        clear();
        RTSEIS_THROW_RTE("%s; Failed to initialize Hilbert transformer",
                         e.what());
    }
    pImpl->mInitialized = true;
}

/// Get the initial condition length
int FIREnvelope::getInitialConditionLength() const
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Envelope class not initialized");
    }
    return pImpl->mImagFIRFilter.getInitialConditionLength();
}

/// Sets the initial conditions
void FIREnvelope::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Envelope class not initialized");
    }
    int nzRef = getInitialConditionLength();
    if (nz != nzRef)
    {
        RTSEIS_THROW_IA("nz = %d must equal %d", nz, nzRef);
    }
    pImpl->mRealFIRFilter.setInitialConditions(nz, zi);
    pImpl->mImagFIRFilter.setInitialConditions(nz, zi);
    pImpl->mHaveInitialCondition = true;
}

/// Resets the initial conditions
void FIREnvelope::resetInitialConditions()
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Envelope class not initialized");
    }
    pImpl->mRealFIRFilter.resetInitialConditions();
    pImpl->mImagFIRFilter.resetInitialConditions();
}

void FIREnvelope::transform(const int n, const double x[], double y[])
{
    if (n < 1){return;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Failed to initialize");
    }
/*
    if (x == nullptr || y)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL";}
        RTSEIS_THROW_IA("%s", "y is NULL"); 
    }
    // Post-processing removes the phase shift
    if (pImpl->mMode = RTSeis::ProcessingMode::POST_PROCESSING)
    {
        // Get the mean

        // Copy and remove the mean 
        int nhalf = pImpl->mNumberOfTaps/2;
        double *xpad  = ippsMalloc_64f(n + nhalf);
        ippsSubC_64f(x, mean, xpad, n);
        ippsZero_64f(&xpad[n], nhalf);
        double *ypadr = ippsMalloc_64f(n + nhalf);
        double *ypadi = ippsMalloc_64f(n + nhalf);
        pImpl->mImagFIRFilter.apply(n+nhalf, xpad, ypadi); 
        ippsFree(xpad);
        // Take the absolute value
        ippsMag_64f(x, &ypadi[nhalf], y, n);
        ippsFree(ypadr);
        ippsFree(ypadi);
      
    }
    else
    {
        // Apply the filter and stitch them together
        for (auto i=0; i<n; i=i+pImpl->chunkSize)
        {
            auto nploc = std::min(n-i, pImpl->chunkSize);
            pImpl->mRealFIRFilter.apply(nploc, &x[i], pImpl->yrWork);
            pImpl->mImagFIRFilter.apply(nploc, &x[i], pImpl->yiWork);
            ippsMag_64f(pImpl->yrWork, pImpl->yiWork, &y[i], nploc);
        }
    }
*/
}

