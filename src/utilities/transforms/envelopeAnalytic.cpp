#include <cstdio>
#include <cstdlib>
#include <complex>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/transforms/envelopeAnalytic.hpp"
#include "rtseis/utilities/transforms/hilbert.hpp"

using namespace RTSeis::Utilities::Transforms;

class EnvelopeAnalytic::EnvelopeAnalyticImpl
{
public: 
    /// Constructor
    EnvelopeAnalyticImpl() = default;
    /// Destructor
    ~EnvelopeAnalyticImpl()
    {
        clear();
    }
    /// Copy constructor
    EnvelopeAnalyticImpl(const EnvelopeAnalyticImpl &envelope)
    {
        *this = envelope;
    }
    /// Clears the memory
    void clear() noexcept
    {
       mHilbert.clear();
       if (mHilb64fc){ippsFree(mHilb64fc);}
       if (mHilb32fc){ippsFree(mHilb32fc);}
       mHilb64fc = nullptr;
       mHilb32fc = nullptr;
       mEnvelopeLength = 0;
       mMean = 0;
       mInitialized = false;
       mPrecision = RTSeis::Precision::DOUBLE;
    }
    /// Copy assignment operator
    EnvelopeAnalyticImpl& operator=(const EnvelopeAnalyticImpl &envelope)
    {
       if (&envelope == this){return *this;}
       clear();
       mHilbert = envelope.mHilbert;
       mEnvelopeLength = envelope.mEnvelopeLength;
       mMean = envelope.mMean;
       mInitialized = envelope.mInitialized;
       mPrecision = envelope.mPrecision;
       if (mEnvelopeLength > 0)
       {
           if (mPrecision == RTSeis::Precision::DOUBLE)
           {
               mHilb64fc = ippsMalloc_64fc(mEnvelopeLength);
               ippsCopy_64fc(envelope.mHilb64fc, mHilb64fc, mEnvelopeLength);
           }
           else
           {
               mHilb32fc = ippsMalloc_32fc(mEnvelopeLength);
               ippsCopy_32fc(envelope.mHilb32fc, mHilb32fc, mEnvelopeLength);
           }
       }
       return *this;
    }
 
//private::
    Transforms::Hilbert mHilbert;
    Ipp64fc *mHilb64fc = nullptr;
    Ipp32fc *mHilb32fc = nullptr;
    double mMean = 0;
    int mEnvelopeLength = 0;
    bool mInitialized = false;
    RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
};

/// Constructor
EnvelopeAnalytic::EnvelopeAnalytic() :
    pImpl(std::make_unique<EnvelopeAnalyticImpl> ())
{
}

/// Copy constructor
EnvelopeAnalytic::EnvelopeAnalytic(const EnvelopeAnalytic &envelope)
{
    *this = envelope;
} 

/// Move constructor
EnvelopeAnalytic::EnvelopeAnalytic(EnvelopeAnalytic &&envelope) noexcept
{
    *this = std::move(envelope);
}

/// Copy assignment operator
EnvelopeAnalytic& EnvelopeAnalytic::operator=(const EnvelopeAnalytic &envelope)
{
    if (&envelope == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<EnvelopeAnalyticImpl> (*envelope.pImpl);
    return *this;
}

/// Move assignment oeprator
EnvelopeAnalytic&
EnvelopeAnalytic::operator=(EnvelopeAnalytic &&envelope) noexcept
{
    if (&envelope == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(envelope.pImpl);
    return *this;
}
/// Destructor
EnvelopeAnalytic::~EnvelopeAnalytic() = default;

/// Clears the module
void EnvelopeAnalytic::clear() noexcept
{
    pImpl->clear();
}

/// Checks if the module is initialized
bool EnvelopeAnalytic::isInitialized() noexcept
{
    return pImpl->mInitialized;
}

/// Gets the envelope length
int EnvelopeAnalytic::getTransformLength()
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    return pImpl->mEnvelopeLength;
}

/// Initializes the class
void EnvelopeAnalytic::initialize(const int n,
                                  const RTSeis::Precision precision)
{
    clear();
    if (n < 1){RTSEIS_THROW_IA("n = %d must be positive", n);}
    if (precision != RTSeis::Precision::DOUBLE)
    {
        RTSEIS_THROW_IA("%s", "Only double precision implemented at moment");
    }
    pImpl->mHilbert.initialize(n, precision); 
    pImpl->mInitialized = true;
}

/// Compute the upper and lower envelope
void EnvelopeAnalytic::transform(const int n, const double x[],
                                 double yupper[], double ylower[])
{
    if (pImpl){pImpl->mMean = 0;}
    if (ylower == nullptr)
    {
        RTSEIS_THROW_IA("%s", "ylower is NULL");
    }
    // Compute the upper envelope
    transform(n, x, yupper); // will throw
    // Compute the lower envelope from the upper envelope.
    // |Hilbert| = yUpper - mean
    // yLower =-|Hilbert| + mean
    // yLower =-(yUpper - mean) + mean =-yUpper + mean
    ippsSubCRev_64f(yupper, pImpl->mMean, ylower, n);
}

/// Compute the envelope
void EnvelopeAnalytic::transform(const int n, const double x[], double y[])
{
    pImpl->mMean = 0;
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    if (n != pImpl->mEnvelopeLength)
    {
        RTSEIS_THROW_IA("n = %d must equal %d", n, pImpl->mEnvelopeLength);
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    // Remove the mean
    double pMean;
    ippsMean_64f(x, n, &pMean);
    ippsSubC_64f(x, pMean, y, n); 
    pImpl->mMean = pMean;
    // Hilbert transform the signal (yields the analytic signal)
    auto yhilb = reinterpret_cast<std::complex<double> *> (pImpl->mHilb64fc);
    pImpl->mHilbert.transform(n, y, yhilb);
    // Take the absolute value to obtain the envelope
    ippsMagnitude_64fc(pImpl->mHilb64fc, y, n); 
    // Restore the mean
    ippsAddC_64f_I(pMean, y, n);
}

/// Compute the envelope
/*
void EnvelopeAnalytic::transform(const int n, const float x[], float y[])
{
    pImpl->mMean = 0;
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    if (n != pImpl->mEnvelopeLength)
    {
        RTSEIS_THROW_IA("n = %d must equal %d", n, pImpl->mEnvelopeLength);
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y is NULL");
    }
    // Remove the mean
    float pMean;
    ippsMean_32f(x, n, &pMean, ippAlgHintAccurate);
    ippsSubC_32f(x, pMean, y, n); 
    pImpl->mMean = static_cast<double> (pMean);
    // Hilbert transform the signal (yields the analytic signal)
    auto yhilb = reinterpret_cast<std::complex<float> *> (pImpl->mHilb32fc);
    pImpl->mHilbert.transform(n, y, yhilb);
    // Take the absolute value to obtain the envelope
    ippsMagnitude_32fc(pImpl->mHilb32fc, y, n); 
    // Restore the mean
    ippsAddC_32f_I(pMean, y, n);
}
*/
