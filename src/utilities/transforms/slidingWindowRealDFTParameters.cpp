#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <ipps.h>
#include "rtseis/utilities/windowFunctions.hpp"
#include "rtseis/utilities/transforms/slidingWindowRealDFTParameters.hpp"
#include "rtseis/utilities/transforms/enums.hpp"

using namespace RTSeis::Utilities::Transforms;

class SlidingWindowRealDFTParameters::SlidingWindowRealDFTParametersImpl
{
public:
    /// Holds the window function
    std::vector<double> mWindow;
    /// The number of samples in the signal
    int mSamples = 0;
    /// The length of the DFT
    int mDFTLength = 0;
    /// The number of samples in the overlap.
    int mSamplesInOverlap = 0;
    /// Describes the window function
    SlidingWindowType mWindowType = SlidingWindowType::BOXCAR;
    /// Describes the detrend strategy
    SlidingWindowDetrendType mDetrendType = SlidingWindowDetrendType::REMOVE_NONE;
};

/// Constructor
SlidingWindowRealDFTParameters::SlidingWindowRealDFTParameters() :
    pImpl(std::make_unique<SlidingWindowRealDFTParametersImpl> ())
{
}

/// Copy constructor
SlidingWindowRealDFTParameters::SlidingWindowRealDFTParameters(
    const SlidingWindowRealDFTParameters &parameters)
{
    *this = parameters;
} 

/// Move constructor
SlidingWindowRealDFTParameters::SlidingWindowRealDFTParameters(
    SlidingWindowRealDFTParameters &&parameters) noexcept
{
    *this = std::move(parameters);
}

/// Copy operator
SlidingWindowRealDFTParameters& 
SlidingWindowRealDFTParameters::operator=(
    const SlidingWindowRealDFTParameters &parameters)
{
    if (&parameters == this){return *this;}
    pImpl = std::make_unique<SlidingWindowRealDFTParametersImpl> (*parameters.pImpl);
    return *this;
}

/// Move operator
SlidingWindowRealDFTParameters&
SlidingWindowRealDFTParameters::operator=(
    SlidingWindowRealDFTParameters &&parameters) noexcept
{
    if (&parameters == this){return *this;}
    pImpl = std::move(parameters.pImpl);
    return *this;
}

/// Destructor
SlidingWindowRealDFTParameters::~SlidingWindowRealDFTParameters() = default;

/// Clears the memory and resets the class
void SlidingWindowRealDFTParameters::clear() noexcept
{
    pImpl->mWindow.clear();
    pImpl->mSamples = 0;
    pImpl->mDFTLength = 0;
    pImpl->mSamplesInOverlap = 0;
    pImpl->mWindowType = SlidingWindowType::BOXCAR;
    pImpl->mDetrendType = SlidingWindowDetrendType::REMOVE_NONE;
}

/// Set number of samples
void SlidingWindowRealDFTParameters::setNumberOfSamples(const int nSamples)
{
    pImpl->mSamples = 0;
    if (nSamples < 1)
    {
        throw std::invalid_argument("Number of samples = " 
                                  + std::to_string(nSamples)
                                  + " must be positive");
    }
    pImpl->mSamples = nSamples; 
}

/// Get number of samples
int SlidingWindowRealDFTParameters::getNumberOfSamples() const
{
    if (pImpl->mSamples < 1)
    {
        throw std::runtime_error("Number of samples not yet set");
    }
    return pImpl->mSamples;
}

/// Sets the window
void SlidingWindowRealDFTParameters::setWindow(
    const int windowLength,
    const SlidingWindowType windowType)
{
    // Error checks
    pImpl->mDFTLength = 0;
    pImpl->mSamplesInOverlap = 0;
    pImpl->mWindow.clear();
    if (windowLength < 1)
    {
        throw std::invalid_argument("Window length must be positive");
    }
    if (windowType == SlidingWindowType::CUSTOM)
    {
        throw std::invalid_argument("Predefined window type cannot be custom");
    } 
    // Resize 
    pImpl->mDFTLength = windowLength; 
    pImpl->mWindowType = windowType;
    pImpl->mWindow.resize(windowLength);
    // Create window
    double *windowPtr = pImpl->mWindow.data();
    if (windowType == SlidingWindowType::HAMMING)
    {
        Utilities::WindowFunctions::hamming(windowLength, &windowPtr);
    }
    else if (windowType == SlidingWindowType::HANN)
    {
        Utilities::WindowFunctions::hann(windowLength, &windowPtr);
    }
    else if (windowType == SlidingWindowType::BLACKMAN)
    {
        Utilities::WindowFunctions::blackman(windowLength, &windowPtr);
    }
    else if (windowType == SlidingWindowType::BARTLETT)
    {
        Utilities::WindowFunctions::bartlett(windowLength, &windowPtr);
    }
    else if (windowType == SlidingWindowType::BOXCAR)
    {
        ippsSet_64f(1.0, windowPtr, windowLength);
    }
#ifndef NDEBUG
    else
    {
        std::cerr << "Impossible location" << std::endl;
        assert(false);
    }
#endif
}

void SlidingWindowRealDFTParameters::setWindow(
    const int windowLength, const double window[])
{
    // Error checks
    pImpl->mDFTLength = 0;
    pImpl->mSamplesInOverlap = 0;
    pImpl->mWindow.clear();
    if (windowLength < 1)
    {
        throw std::invalid_argument("Window length must be positive");
    }
    if (window == nullptr)
    {
        throw std::invalid_argument("window cannot be NULL");
    }
    // Resize and copy
    pImpl->mDFTLength = windowLength;
    pImpl->mWindow.resize(windowLength);
    pImpl->mWindowType = SlidingWindowType::CUSTOM;
    ippsCopy_64f(window, pImpl->mWindow.data(), windowLength);
}

std::vector<double> SlidingWindowRealDFTParameters::getWindow() const
{
    if (pImpl->mWindow.empty()){throw std::runtime_error("Window not yet set");}
    return pImpl->mWindow;
}

int SlidingWindowRealDFTParameters::getWindowLength() const
{
    if (pImpl->mWindow.empty()){throw std::runtime_error("Window not yet set");}
    return static_cast<int> (pImpl->mWindow.size());
}

SlidingWindowType SlidingWindowRealDFTParameters::getWindowType() const
{
    if (pImpl->mWindow.empty()){throw std::runtime_error("Window not yet set");}
    return pImpl->mWindowType;
}

/// Sets the DFT length
void SlidingWindowRealDFTParameters::setDFTLength(const int dftLength)
{
    int windowLength = getWindowLength(); // Will throw if window not set
    pImpl->mDFTLength = windowLength;
    if (dftLength < windowLength)
    {
        throw std::invalid_argument("dftLength = " + std::to_string(dftLength)
                                  + " cannot be less than window length = "
                                  + std::to_string(windowLength));
    }
    pImpl->mDFTLength = dftLength;
} 

int SlidingWindowRealDFTParameters::getDFTLength() const
{
    if (pImpl->mWindow.empty()){throw std::runtime_error("Window not yet set");}
    return pImpl->mDFTLength;
}

/// Sets the number of samples in the ovlerap
void SlidingWindowRealDFTParameters::setNumberOfSamplesInOverlap(
    const int nSamplesInOverlap)
{
    pImpl->mSamplesInOverlap = 0;
    int windowLength = getWindowLength(); // Will throw
    if (nSamplesInOverlap < 0 || nSamplesInOverlap >= windowLength)
    {
        throw std::invalid_argument("nSamplesInOverlap = "
                                  + std::to_string(nSamplesInOverlap)
                                  + " must be in range [0,"
                                  + std::to_string(windowLength - 1) + "]");
    }
    pImpl->mSamplesInOverlap = nSamplesInOverlap;
}

int SlidingWindowRealDFTParameters::getNumberOfSamplesInOverlap() const noexcept
{
    return pImpl->mSamplesInOverlap;
}

/// Defines the detrend type
void SlidingWindowRealDFTParameters::setDetrendType(
    const SlidingWindowDetrendType detrendType) noexcept
{
    pImpl->mDetrendType = detrendType;
}

SlidingWindowDetrendType
SlidingWindowRealDFTParameters::getDetrendType() const noexcept
{
    return pImpl->mDetrendType;
}

/// Gets/sets the precision
/*
void SlidingWindowRealDFTParameters::setPrecision(
    const RTSeis::Precision precision) noexcept 
{
    pImpl->mPrecision = precision;
}

RTSeis::Precision SlidingWindowRealDFTParameters::getPrecision() const noexcept
{
    return pImpl->mPrecision;
}
*/

/// Check if class is usable
bool SlidingWindowRealDFTParameters::isValid() const noexcept
{
    if (pImpl->mSamples < 1){return false;}
    if (pImpl->mWindow.empty()){return false;}
    return true;
}

/// Check for equality
bool RTSeis::Utilities::Transforms::operator==(
    const SlidingWindowRealDFTParameters &lhs,
    const SlidingWindowRealDFTParameters &rhs)
{
    if (lhs.isValid() != rhs.isValid()){return false;}
    // Do the no except checks
    if (lhs.getDetrendType() != rhs.getDetrendType()){return false;}
    if (lhs.getNumberOfSamplesInOverlap() != rhs.getNumberOfSamplesInOverlap())
    {
        return false;
    }
    // Do the checks that would throw exceptions
    if (lhs.isValid() == false && rhs.isValid() == false){return true;}
    // This is valid so we can safely call these
    if (lhs.getNumberOfSamples() != rhs.getNumberOfSamples()){return false;}
    auto v1 = lhs.getWindow();
    auto v2 = rhs.getWindow();
    if (v1.size() != v2.size()){return false;}
    double error = 0;
    ippsNormDiff_Inf_64f(v1.data(), v2.data(), v1.size(), &error);
    if (error > 0){return false;}
    return true;
}

bool RTSeis::Utilities::Transforms::operator!=(
    const SlidingWindowRealDFTParameters &lhs,
    const SlidingWindowRealDFTParameters &rhs)
{
    return !(lhs == rhs);
}
