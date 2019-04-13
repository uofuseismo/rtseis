#ifndef RTSEIS_UTILITIES_TRANSFORMS_ENUMS_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_ENUMS_HPP 1

namespace RTSeis
{
namespace Utilities
{
namespace Transforms 
{
/*!
 * @brief Defines the Fourier transform implementation.
 */
enum class FourierTransformImplementation
{
    DFT, /*!< Perform a Discrete Fourier Transform computation. */
    FFT  /*!< Force an Fast Fouerier Transform computation.  The 
              implementation will have to zero-pad the signal so that
              its length is a power of 2. */
};

}; // End Transforms
}; // End Utilities
}; // End RTSeis

#endif
