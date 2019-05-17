#include <cstdio>
#include <cstdlib>
#include <array>
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
        mMean = firEnvelope.mMean;
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
    double mMean = 0;
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
    pImpl->mMean = 0;
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
    pImpl->mMean = 0;
    if (n < 1){return;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Failed to initialize");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    // Post-processing removes the phase shift
    if (pImpl->mMode == RTSeis::ProcessingMode::POST_PROCESSING)
    {
        // Compute the mean
        double pMean;
        ippsMean_64f(x, n, &pMean);
        pImpl->mMean = pMean;
        // Remove the mean and pad out the signal
        // N.B. The group delay is actually + 1 but C wants to shift relative to
        // a base address so we subtract the one.  Hence, n/2 instead of n/2+1.
        int groupDelay = pImpl->mNumberOfTaps/2;
        int npad = n + groupDelay;
        double *xPad = ippsMalloc_64f(npad);
        ippsSubC_64f(x, pMean, xPad, n);
        ippsZero_64f(&xPad[n], groupDelay); // Post-pad with zeros
        // Now apply the filter and compute absolute value - Type III
        if (pImpl->mType3)
        {
            double *yPadr = xPad;
            double *yPadi = ippsMalloc_64f(npad);
            pImpl->mImagFIRFilter.apply(npad, xPad, yPadi); 
            ippsMagnitude_64f(yPadr, &yPadi[groupDelay], y, n);
            ippsFree(yPadi);
        }
        else
        {
            double *yPadr = ippsMalloc_64f(npad);
            pImpl->mRealFIRFilter.apply(npad, xPad, yPadr);
            double *yPadi = ippsMalloc_64f(npad);
            pImpl->mImagFIRFilter.apply(npad, xPad, yPadi);
            ippsMagnitude_64f(&yPadr[groupDelay], &yPadi[groupDelay], y, n);
            ippsFree(yPadr);
            ippsFree(yPadi);
        }
        ippsFree(xPad);
        // Reconstitute the mean
        ippsAddC_64f_I(pMean, y, n);
    }
    else
    {
        // TODO - add a sparse FIR filter for the type III case
        constexpr int chunkSize = 1024;
        std::array<double, chunkSize> yrTemp;
        std::array<double, chunkSize> yiTemp;
        for (auto ic=0; ic<n; ic=ic+chunkSize)
        {
            auto npfilt = std::min(n - ic, chunkSize);
            pImpl->mRealFIRFilter.apply(npfilt, &x[ic], yrTemp.data());
            pImpl->mImagFIRFilter.apply(npfilt, &x[ic], yiTemp.data());
            ippsMagnitude_64f(yrTemp.data(), yiTemp.data(), &y[ic], npfilt);
        }
    }
}

