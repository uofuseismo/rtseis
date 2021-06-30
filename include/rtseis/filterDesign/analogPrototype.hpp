#ifndef RTSEIS_FILTERDESIGN_IIR_AP_HPP
#define RTSEIS_FILTERDESIGN_IIR_AP_HPP 1
#include <complex>
#include <vector>
namespace RTSeis::FilterRepresentations
{
class ZPK;
}
/// @defgroup rtseis_utils_design_iir_ap Analog Prototype
/// @brief Utility functions for analog prototype IIR filter design.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_filterdesign_iir
namespace RTSeis::FilterDesign::IIR::AnalogPrototype
{
/// @brief Computes an n'th order analog prototype Bessel filter design.
/// @param[in] n   The filter order which is the number of poles.
///                This must be positive.
/// @result The corresonding Bessel filter in zero, pole, gain format.
/// @throws std::invalid_argument if the order is invalid.
/// @ingroup rtseis_filterdesign_iir_ap
[[nodiscard]] FilterRepresentations::ZPK bessel(int n);
/// @brief Computes an n'th order analog prototype Butterworth filter design.
/// @param[in] n   The filter order which is the number of poles.
///                This must be positive.
/// @result The corresonding Butterworth filter in zero, pole, gain format.
/// @throws std::invalid_argument if the order is invalid.
/// @ingroup rtseis_filterdesign_iir_ap
[[nodiscard]] FilterRepresentations::ZPK butter(int n);
/// @brief Computes an n'th order Chebyshev I filter with \f$ r_p \f$
///        decibels of ripple in the passband.  This is based on the
///        SciPy implementation.
/// @param[in] n   The order of the filter.  This must be positive.
/// @param[in] rp  Controls the size of the ripples in the passband.
///                This is measured in decibels and must be positive.
/// @result The corresonding Chebyshev I filter in zero, pole, gain format.
/// @throws std::invalid_argument if the order or rp is invalid.
/// @ingroup rtseis_filterdesign_iir_ap
[[nodiscard]] FilterRepresentations::ZPK cheb1ap(int n, double rp);
/// @brief Computes an n'th order Chebyshev II filter with \f$ r_s \f$
///        decibels of ripple in the stopband.  This is based on the
///        SciPy implementation.
/// @param[in] n   The order of the filter.  This must be positive.
/// @param[in] rs  Controls the size of the ripples in the stopband.
///                This is measured in decibels and must be positive.
/// @result The corresonding Chebyshev II filter in zero, pole, gain format.
/// @throws std::invalid_argument if the order or rs is invalid.
/// @ingroup rtseis_filterdesign_iir_ap
[[nodiscard]] FilterRepresentations::ZPK cheb2ap(int n, double rs);
} 
#endif
