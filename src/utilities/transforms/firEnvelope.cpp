#include <cstdio>
#include <cstdlib>
#include "rtseis/private/throw.hpp"

using namespace RTSeis::Utilities;
using namespace RTSeis::Utilities::Transforms;

class FIREnvelope::FIREnvelopeImpl
{
public:
    void initialize(const int ntaps,
                    const RTSeis::ProcessingMode mode,
                    const RTSeis::Precision precision)
    {
        mNumberOfTaps = ntaps;
        std::pair<FilterRepresentations::FIR, FilterRepresentations::FIR> zfir;
        constexpr double beta = 8;
        zfir = HilbertTransformer(pImpl->mNumberOfTaps - 1, beta);
        std::vector<double> realfir = zfir.first().getFilterTaps();
        std::vector<double> imagfir = zfir.second().getFilterTaps();
        realFIRFilter.initialize(realfir.size(), realfir.data(),
                                 mode, precision,
                                 FIRImplementation::DIRECT);
        imagFIRFilter.initialize(imagfir.size(), imagfir.data(),
                                 mode, precision,
                                 FIRImplementation::DIRECT);

    }

    FilterImplementations::FIRFilter realFIRFilter;
    FilterImplementations::FIRFilter imagFIRFilter;
    
    int mNumberOfTaps = 0;
    bool mZeroPhase = true;
    bool mType3 = false;
    bool mInitialized = false;
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
FIREnvelope::FIREnvelope(FIREnvelope &&firEnvelope)
{
    *this = firEnvelope;
}

/// Copy assignment
FIREnvelope& FIREnvelope::operator=(const FIREnvelope &firEnvelope)
{
    if (&firEnvelope == this){return *this;}
   
}

/// Move assignment
FIREnvelope& FIREnvelope::operator=(FIREnvelope &&firEnvelope)
{
    if (&firEnvelope == this){return *this;}
    pImpl = std::move(firEnvelope.pImpl);
}

/// Destructor
FIREnvelope::~FIREnvelope() = default;

/// Clear the filter
void FIREnvelope::clear() noexcept
{
    pImpl->realFIRFilter.clear();
    pImpl->imagFIRFilter.clear();
    pImpl->mNumberOfTaps = 0;
    pImpl->mZeroPhase = true;
    pImpl->isType3 = false;
    pImpl->mInitialized = false;
}

/// Check if initialized
void FIREnvelope::isInitialized() const noexcept
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
    pImpl->mType3 = false;
    if (ntaps%2 == 1){pImpl->mType3 = true;}
    pImpl->mNumberOfTaps = ntaps;
    std::pair<FilterRepresentations::FIR, FilterRepresentations::FIR> zfir;
    constexpr double beta = 8;
    // Create an FIR hilbert transform
    try
    {
        zfir = HilbertTransformer(pImpl->mNumberOfTaps - 1, beta);
        std::vector<double> rfir = zfir.first().getFilterTaps();
        pImpl->realFIRFilter.initialize(rfir.size(), rfir.data(),
                                        mode, precision,
                                        FIRImplementation::DIRECT);
        std::vector<double> cfir = zfir.second().getFilterTaps();
        pImpl->imagFIRFilter.initialize(cfir.size(), cfir.data(),
                                        mode, precision,
                                        FIRImplementation::DIRECT);
    }
    else
    {
        clear();
        RTSEIS_THROW_RTE("%s", "Failed to initialize Hilbert transformer");
    }
    pImpl->mInitialzed = true;
}

void FIREnvelope::apply(const int n, const double x[], double y[])
{
    if (n < 1){return;} // Nothing to do
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Failed to initialize");
    }
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
        pImpl->imagFIRFilter.apply(n+nhalf, xpad, ypadi); 
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
            pImpl->realFIRFilter.apply(nploc, &x[i], pImpl->yrWork);
            pImpl->imagFIRFilter.apply(nploc, &x[i], pImpl->yiWork);
            ippsMag_64f(pImpl->yrWork, pImpl->yiWork, &y[i], nploc);
        }
    }
}

