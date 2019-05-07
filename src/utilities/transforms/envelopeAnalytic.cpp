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
    ~EnvelopeAnalyticImpl() = default;
    /// Copy constructor
    EnvelopeAnalyticImpl(const EnvelopeAnalyticImpl &envelope)
    {
       *this = envelope;
    }
    /// Clears the memory
    void clear() noexcept
    {
       mHilbert.clear();
       mInitialized = false;
    }
    /// Copy assignment operator
    EnvelopeAnalyticImpl& operator=(const EnvelopeAnalyticImpl &envelope)
    {
       if (&envelope == this){return *this;}
       mHilbert = envelope.mHilbert;
       mInitialized = envelope.mInitialized;
       return *this;
    }
 
//private::
    Transforms::Hilbert mHilbert;
    bool mInitialized = false;
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

/// Initializes the class
void EnvelopeAnalytic::initialize(const int n,
                                  const RTSeis::Precision precision)
{
    clear();
    if (n < 1){RTSEIS_THROW_IA("n = %d must be positive", n);}
    pImpl->mHilbert.initialize(n, precision); 
    pImpl->mInitialized = true;
}
