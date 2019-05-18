#include <cstdio>
#include <cstdlib>
#include <complex>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/transforms/envelope.hpp"
#include "rtseis/utilities/transforms/hilbert.hpp"

using namespace RTSeis::Utilities::Transforms;

class Envelope::EnvelopeImpl
{
public: 
    /// Constructor
    EnvelopeImpl() = default;
    /// Destructor
    ~EnvelopeImpl()
    {
        clear();
    }
    /// Copy constructor
    EnvelopeImpl(const EnvelopeImpl &envelope)
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
    EnvelopeImpl& operator=(const EnvelopeImpl &envelope)
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
Envelope::Envelope() :
    pImpl(std::make_unique<EnvelopeImpl> ())
{
}

/// Copy constructor
Envelope::Envelope(const Envelope &envelope)
{
    *this = envelope;
} 

/// Move constructor
Envelope::Envelope(Envelope &&envelope) noexcept
{
    *this = std::move(envelope);
}

/// Copy assignment operator
Envelope& Envelope::operator=(const Envelope &envelope)
{
    if (&envelope == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<EnvelopeImpl> (*envelope.pImpl);
    return *this;
}

/// Move assignment oeprator
Envelope&
Envelope::operator=(Envelope &&envelope) noexcept
{
    if (&envelope == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(envelope.pImpl);
    return *this;
}
/// Destructor
Envelope::~Envelope() = default;

/// Clears the module
void Envelope::clear() noexcept
{
    pImpl->clear();
}

/// Checks if the module is initialized
bool Envelope::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Gets the envelope length
int Envelope::getTransformLength() const
{
    if (!isInitialized())
    {
        RTSEIS_THROW_RTE("%s", "Class not initialized");
    }
    return pImpl->mEnvelopeLength;
}

/// Initializes the class
void Envelope::initialize(const int n,
                          const RTSeis::Precision precision)
{
    clear();
    if (n < 1){RTSEIS_THROW_IA("n = %d must be positive", n);}
    if (precision != RTSeis::Precision::DOUBLE)
    {
        RTSEIS_THROW_IA("%s", "Only double precision implemented at moment");
    }
    if (n > 1)
    {
        pImpl->mHilbert.initialize(n, precision); 
    }
    if (precision == RTSeis::Precision::DOUBLE)
    {
        pImpl->mHilb64fc = ippsMalloc_64fc(n);
    }
    else
    {
        pImpl->mHilb32fc = ippsMalloc_32fc(n);
    }
    pImpl->mEnvelopeLength = n;
    pImpl->mPrecision = precision;
    pImpl->mInitialized = true;
}

/// Compute the upper and lower envelope
void Envelope::transform(const int n, const double x[],
                         double yupper[], double ylower[])
{
    if (pImpl){pImpl->mMean = 0;}
    if (ylower == nullptr)
    {
        RTSEIS_THROW_IA("%s", "ylower is NULL");
    }
    // Compute the upper envelope
    transform(n, x, yupper); // will throw
    // Special case
    if (n == 1)
    {
        ylower[0] = yupper[0];
        return;
    } 
    // Compute the lower envelope from the upper envelope.
    //  (1) yupper = |Hilbert| + mean
    //  (2) ylower = |Hilbert| - mean
    // Expressing |Hilbert| as a function of yupper in Eqn 2
    //   |Hilbert| = yupper - mean
    // And substituting into the second equation
    //   ylower = yupper - mean - mean 
    //          = yupper - 2*mean
    ippsSubCRev_64f(yupper, 2*pImpl->mMean, ylower, n);
}

/// Compute the envelope
void Envelope::transform(const int n, const double x[], double y[])
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
    // Special case
    if (n == 1)
    {
        pImpl->mMean = x[0];
        y[0] = x[0];
        return;
    }
    // Remove the mean
    double pMean;
    ippsMean_64f(x, n, &pMean);
    ippsSubC_64f(x, pMean, y, n); 
    pImpl->mMean = pMean;
    // Hilbert transform the signal (yields the analytic signal)
    auto yhilb = reinterpret_cast<std::complex<double> *> (pImpl->mHilb64fc);
if (yhilb == nullptr){printf("problem\n");}
    pImpl->mHilbert.transform(n, y, yhilb);
    // Take the absolute value to obtain the envelope
    ippsMagnitude_64fc(pImpl->mHilb64fc, y, n); 
    // Restore the mean
    ippsAddC_64f_I(pMean, y, n);
}

/// Compute the envelope
/*
void Envelope::transform(const int n, const float x[], float y[])
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
