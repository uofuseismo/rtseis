#ifndef RTSEIS_UTILS_TRANSFORMS_UTILITIES_HPP
#define RTSEIS_UTILS_TRANSFORMS_UTILITIES_HPP 1
#include <cmath>
#include <complex>

namespace RTSeis
{
namespace Utilities
{
namespace Transforms 
{
/*!
 * @defgroup rtseis_utils_transforms_utils Utility Functions
 * @brief Utility routines accompanying the DFT.
 * @ingroup rtseis_utils_transforms
 */
namespace DFTUtilities
{
/*! 
 * @brief Unwraps the phase, p, by changing the absolute jumps
 *        greater than \f$ \pi \f$ to their \f$ 2 \pi \f$ 
 *        complement.  If tolerance is set, then the jump is
 *        \f$ 2 \cdot tol \f$.
 * @param[in] n    The number of points in the phase signals.
 * @param[in] p    The array of phase angles in radians to unwrap.
 *                 This has dimension [n].
 * @param[in] q    The unwrapped phase angles in radians.  This
 *                 has dimension [n].
 * @param[in] tol  The jump tolerance specified in radians which
 *                 must be positive.  The default is \f$ \pi \f$.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_transforms_utils
 */
int unwrap(const int n, const double p[], double q[], 
           const double tol = M_PI);

/*!
 * @brief Computes the phase angle, i.e., the angle between the
 *        imaginary and real parts of z 
 *        \f[ \phi
 *          = \arctan \left ( \frac{ \Im \{z \} }
 *                                 { \Re \{z \} } \right ) \f].
 * @param[in] n     Number of points in array z.
 * @param[in] z     The array of complex numbers of which to compute
 *                  the angle.  This has dimension [n].
 * @param[out] phi  The phase angle between the real and imaginary 
 *                  parts of z.  This has dimension [n].  By
 *                  default this will be in radians but it can
 *                  optionally be given in degrees.
 * @param[in] lwantDeg  If true then phi is given in degrees.  
 *                      If false then phi is given in radians.
 *                      This is the default.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_transforms_utils
 */
int phase(const int n, const std::complex<double> z[], double phi[],
          const bool lwantDeg = false);

/*!
 * @brief Finds the next number, n2, such that n2 is a power of 2 and
 *        n2 is greater than or equal to n.
 * @param[in] n  Non-negative number of which to find the next power
 *               of 2.
 * @result On successful exit this is a number that is a power of 2
 *         and is greater than or equal to n.
 * @ingroup rtseis_utils_transforms_utils
 */
int nextPow2(const int n);

}; // End DFT utilities
}; // End Transforms
}; // End Utilities
}; // End RTSeis
#endif
