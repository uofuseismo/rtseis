#ifndef RTSEIS_UTILS_DESIGN_ENUMS_HPP
#define RTSEIS_UTILS_DESIGN_ENUMS_HPP 1

namespace RTSeis
{
namespace Utilities
{
namespace FilterDesign
{

/*! 
 * @brief defines the windows available for FIR design.
 * @ingroup rtseis_utils_design_fir
 */
enum class FIRWindow
{
    HAMMING,     /*!< Hamming window. */
    BARTLETT,    /*!< Bartlett (triangle) window. */
    HANN,        /*!< Hann window. */
    BLACKMAN_OPT /*!< Optimal Blackman window. */
};
/*! 
 * @brief Defines the analog prototype for the IIR filter design.
 * @ingroup rtseis_utils_design_iir
 */
enum class IIRPrototype
{   
    BUTTERWORTH, /*!< Butterworth filter design. */
    BESSEL,      /*!< Bessel filter design. */
    CHEBYSHEV1,  /*!< Chebyshev I filter design. */
    CHEBYSHEV2   /*!< Chebyshev II filter design. */
};  
/*!
 * @brief Differentiates the IIR filter design as being for analog or
 *        digital filtering.
 * @ingroup rtseis_utils_design_iir
 */
enum class IIRFilterDomain
{
    DIGITAL, /*!< Digital filter. */
    ANALOG   /*!< Analog filter. */
};
/*! 
 * @brief Defines the bandtype for the IIR filter design.
 * @ingroup rtseis_utils_design_iir
 */
enum class Bandtype
{   
    LOWPASS,  /*!< Lowpass filter. */
    HIGHPASS, /*!< Highpass filter. */
    BANDPASS, /*!< Bandpass filter. */
    BANDSTOP  /*!< Bandstop filter. */
};  
/*! 
 * @brief Pairing strategy when converting a ZPK filter to an SOS filter.
 * @ingroup rtseis_utils_design_iir
 */
enum class SOSPairing
{
    NEAREST = 0,  /*!< This attempts to minimize the peak gain. */
    KEEP_ODD = 1  /*!< This attempts to minimize the peak gain
                       subject to the constraint that odd-order
                       systems should retain one section as first order. */
};

}; // End Filter design
}; // End Utilities
}; // End RTSeis

#endif
