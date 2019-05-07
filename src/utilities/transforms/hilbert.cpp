#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <valarray>
#ifdef __INTEL_COMPILER
#include <pstl/algorithm>
#include <pstl/execution>
#endif
#define RTSEIS_LOGGING 1
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/transforms/enums.hpp"
#include "rtseis/utilities/transforms/hilbert.hpp"
#include "rtseis/utilities/transforms/dftRealToComplex.hpp"
#include "rtseis/utilities/transforms/dft.hpp"


using namespace RTSeis::Utilities::Transforms;

class Hilbert::HilbertImpl
{
public:
    /// Constructor
    HilbertImpl() = default;
    /// Destructor
    ~HilbertImpl() = default;
    /// Copy constructor
    HilbertImpl(const HilbertImpl &hilbert)
    {
        *this = hilbert;
    }
    /// Copy assignment operator
    HilbertImpl& operator=(const HilbertImpl &hilbert)
    {
        if (&hilbert == this){return *this;}
        mInitializedFlag = hilbert.mInitializedFlag; 
        if (!mInitializedFlag){return *this;}
        mDFTR2C = hilbert.mDFTR2C;
        mDFT = hilbert.mDFT;
        mTransformLength = hilbert.mTransformLength;
        mPrecision = hilbert.mPrecision;
        return *this;
    }
    /// Clears the memory
    void clear() noexcept
    {
        mDFTR2C.clear();
        mDFT.clear();
        mTransformLength = 0;
        mInitializedFlag = false;
        mPrecision = RTSeis::Precision::DOUBLE;
    }
    /// Initialize the DFT's
    void initialize(const int n, RTSeis::Precision precision)
    {
        mTransformLength = n;
        mPrecision = precision; 
        mDFTR2C.initialize(n, FourierTransformImplementation::DFT, mPrecision);
        mDFT.initialize(n, FourierTransformImplementation::DFT, mPrecision);
        mInitializedFlag = true;
    }
    /// Apply the Hilbert transform
    void transform(const double xr[], std::complex<double> h[])
    {
        // Transform the input signal
        constexpr std::complex<double> zero = std::complex<double> (0, 0);
        std::valarray<std::complex<double>> x(zero, mTransformLength);
        auto xptr = std::begin(x); 
        mDFTR2C.forwardTransform(mTransformLength, xr,
                                 mTransformLength, xptr);
        // Compute the analytic signal see: 
        // https://en.wikipedia.org/wiki/Analytic_signal
        int nw = (mTransformLength + 1)/2;
        if (mTransformLength%2 == 0){nw = mTransformLength/2 - 1;}
        // (i)   S_a(f) = S(f) for f = 0
        // (ii)  S_a(f) = 2*S(f) for f > 0
        // (iii) S_a(f) = 0 for f < 0
        int i1 = 1;
        int i2 = i1 + nw;
#ifdef __INTEL_COMPILER
        std::for_each(std::execution::unseq,
                      std::begin(x)+i1, std::begin(x)+i2,
                      [](auto &xi){ xi *= 2.0;});
#else
        std::for_each(std::begin(x)+i1, std::begin(x)+i2,
                      [](auto &xi){ xi *= 2.0;});
#endif
        // Inverse transform to obtain the time domain signal
        xptr = std::begin(x);
        mDFT.inverseTransform(mTransformLength, xptr, mTransformLength, h);
    }
    /// Apply the Hilbert transform
    void transform(const float xr[], std::complex<float> h[])
    {
        // Transform the input signal
        constexpr std::complex<float> zero = std::complex<float> (0, 0); 
        std::valarray<std::complex<float>> x(zero, mTransformLength);
        auto xptr = std::begin(x); 
        mDFTR2C.forwardTransform(mTransformLength, xr, 
                                 mTransformLength, xptr);
        // Compute the analytic signal see: 
        // https://en.wikipedia.org/wiki/Analytic_signal
        int nw = (mTransformLength + 1)/2;
        if (mTransformLength%2 == 0){nw = mTransformLength/2 - 1;} 
        // (i)   S_a(f) = S(f) for f = 0
        // (ii)  S_a(f) = 2*S(f) for f > 0
        // (iii) S_a(f) = 0 for f < 0
        int i1 = 1;
        int i2 = i1 + nw; 
#ifdef __INTEL_COMPILER
        std::for_each(std::execution::unseq,
                      std::begin(x)+i1, std::begin(x)+i2,
                      [](auto &xi){ xi *= 2.0f;});
#else
        std::for_each(std::begin(x)+i1, std::begin(x)+i2,
                      [](auto &xi){ xi *= 2.0f;});
#endif
        // Inverse transform to obtain the time domain signal
        xptr = std::begin(x);
        mDFT.inverseTransform(mTransformLength, xptr, mTransformLength, h);
    }
    /// Variables
    DFTRealToComplex mDFTR2C;
    DFT mDFT;
    int mTransformLength = 0;
    RTSeis::Precision mPrecision = Precision::DOUBLE;
    bool mInitializedFlag = false;
};

Hilbert::Hilbert() :
    pImpl(std::make_unique<HilbertImpl> ())
{
}

Hilbert::Hilbert(const Hilbert &hilbert)
{
    *this = hilbert;
}

/*
Hilbert::Hilbert(Hilbert &&hilbert) noexcept
{
    *this = std::move(hilbert);
}
*/

Hilbert& Hilbert::operator=(const Hilbert &hilbert)
{
    if (&hilbert == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<HilbertImpl> (*hilbert.pImpl);
    return *this;
}

/*
Hilbert& Hilbert::operator=(Hilbert &&hilbert) noexcept
{
    if (&hilbert == this){return *this;}
    pImpl = std::move(hilbert.pImpl);
    return *this;
}
*/

Hilbert::~Hilbert() = default;

void Hilbert::clear(void) noexcept
{
    pImpl->clear();
}

bool Hilbert::isInitialized() noexcept
{
    return pImpl->mInitializedFlag;
}

int Hilbert::getTransformLength()
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    return pImpl->mTransformLength;
}

void Hilbert::initialize(const int n,
                         const RTSeis::Precision precision)
{
    clear();
    if (n < 1)
    {
        RTSEIS_THROW_IA("n = %d must be positive", n);
    }
    pImpl->initialize(n, precision);
}

void Hilbert::transform(const int n, const double x[], 
                        std::complex<double> h[])
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (n != pImpl->mTransformLength)
    {
        RTSEIS_THROW_IA("n = %d must equal %d", n, pImpl->mTransformLength);
    }
    if (pImpl->mPrecision != RTSeis::Precision::DOUBLE)
    {
        RTSEIS_THROW_RTE("%s", "Precision switch not yet implemented");
    }
    pImpl->transform(x, h);
    return;
}


