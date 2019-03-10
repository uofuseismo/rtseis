#ifndef RTSEIS_UTILS_DESIGN_FIR_HPP
#define RTSEIS_UTILS_DESIGN_FIR_HPP 1
#include <utility>
#include "rtseis/utilities/design/enums.hpp"

namespace RTSeis
{
namespace Utilities
{
namespace FilterRepresentations
{
class FIR;
}
/*!
 * @defgroup rtseis_utils_filterDesign Filter Design
 * @brief Utility functions for FIR and IIR filter design.
 * @ingroup rtseis_utils
 */
namespace FilterDesign
{
/*!
 * @defgroup rtseis_utils_design_fir FIR Design
 * @brief Utility functions for window-based FIR filter design.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filterDesign
 */
namespace FIR
{

/*!
 * @brief Designs an FIR lowpass filter using the window method.
 * @param[in] order   Order of filter.  The number of taps is order + 1.
 *                    This must be at least 4.
 * @param[in] r       The normalized cutoff frequency where 1 is
 *                    the Nyquist.
 * @param[out] fir    The FIR lowpass filter corresponding to the
 *                    design parameters.
 * @param[in] window  FIR window design.  The default is a Hamming window.
 * @throws std::invalid_argument if any arguments are incorrect.
 * @ingroup rtseis_utils_design_fir
 */
void FIR1Lowpass(const int order, const double r,
                 FilterRepresentations::FIR &fir,
                 const FIRWindow window = FIRWindow::HAMMING);
/*!
 * @brief Designs an FIR highpass filter using the window method.
 * @param[in] order   Order of filter.  The number of taps is order + 1.
 *                    This must be at least 4.
 * @param[in] r       The normalized cutoff frequency where 1 is
 *                    the Nyquist.
 * @param[out] fir    The highpass filter corresponding to the 
 *                    design parameters.
 * @param[in] window  FIR window design.  The default is a Hamming window.
 * @throws std::invalid_argument if any arguments are incorrect.
 * @ingroup rtseis_utils_design_fir
 */
void FIR1Highpass(const int order, const double r,
                  FilterRepresentations::FIR &fir,
                  const FIRWindow window = FIRWindow::HAMMING);
/*!
 * @brief Designs an FIR bandpass filter using the window method.
 * @param[in] order   Order of filter.  The number of taps is order + 1.
 *                    This must be at least 4.
 * @param[in] r       Normalized low cutoff frequency and high cutoff
 *                    frequencies where 1 is the Nyquist.  Here,
 *                    r.first is the low cutoff and r.second is the
 *                    high cutoff.
 * @param[out] fir    The bandpass filter corresponding to the design
 *                    parameters.
 * @param[in] window  FIR window design.  The default is a Hamming window.
 * @throws std::invalid_argument if any arguments are incorrect.
 * @ingroup rtseis_utils_design_fir
 */
void FIR1Bandpass(const int order, const std::pair<double,double> &r,
                  FilterRepresentations::FIR &fir,
                  const FIRWindow window = FIRWindow::HAMMING);
/*! 
 * @brief Designs an FIR bandstop (notch) filter using the window method.
 * @param[in] order   Order of filter.  The number of taps is order + 1.
 *                    This must be at least 4.
 * @param[in] r       Normalized low cutoff frequency and high cutoff
 *                    frequencies where 1 is the Nyquist.  Here,
 *                    r.first is the low cutoff and r.second is the
 *                    high cutoff.
 * @param[out] fir    The bnadstop filter corresponding to the design
 *                    parameters.
 * @param[in] window  FIR window design.  The default is a Hamming window.
 * @throws std::invalid_argument if any arguments are incorrect.
 * @ingroup rtseis_utils_design_fir
 */
void FIR1Bandstop(const int order, const std::pair<double,double> &r,
                  FilterRepresentations::FIR &fir,
                  const FIRWindow window = FIRWindow::HAMMING);

}; /* End FIR */

}; /* End FilterDesign */
}; /* End Utils */

}; /* End RTseis */

#endif
