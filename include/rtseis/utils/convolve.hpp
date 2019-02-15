#ifndef RTSEIS_UTILS_POLYNOMIAL_HPP
#define RTSEIS_UTILS_POLYNOMIAL_HPP 1
#include "rtseis/config.h"
#include <complex>
#include <vector>

namespace RTSeis
{
namespace Utilities
{
namespace Math 
{
namespace Convolve
{
    enum class Implementation
    {
        AUTO = 0,   /*!< Let IPP decide. */
        DIRECT = 1, /*!< Time domain implementation. */
        FFT = 2     /*!< FFT implementation. */
    };
    enum class Mode
    {
        FULL = 0,   /*!< A full discrete convolution or correlation of
                         inputs which will have length \f$ m + n - 1 \f$.
                         Because the signals do not overlap completely
                         at the convolution edges boundary effects can be
                         seen. */
        VALID = 1,  /*!< The output consists only of those elements that
                         do not rely on zero-padding.  The return 
                         convolution or correlation will have length
                         \f$ \max(m,n) - \min(m,n) + 1 \f$.  This will
                         only be computed where the input signals completely
                         overlap so that there will not be edge effects. */
        SAME  = 2,  /*!< The output is the same size as the first input
                         and centered with respect to the FULL output.
                         The resulting convolution or correlation will
                         have length \f$ \max(m, n) \f$. */
    };
         
    int convolve(const std::vector<double> a,
                 const std::vector<double> b,
                 std::vector<double> &c,
                 const Mode mode = Mode::FULL,
                 const Implementation implementation = Implementation::AUTO);
    int correlate(const std::vector<double> a,
                  const std::vector<double> b,
                  std::vector<double> &c, 
                  const Mode mode = Mode::FULL,
                  const Implementation implementation = Implementation::AUTO);
    int autocorrelate(const std::vector<double> a,
                      std::vector<double> &c, 
                      const Mode mode = Convolve::Mode::FULL,
                      const Implementation implementation = Implementation::AUTO);
}; /* End Convolve. */
}; /* End Math. */
}; /* End Utils */
}; /* End RTSeis */
#endif
