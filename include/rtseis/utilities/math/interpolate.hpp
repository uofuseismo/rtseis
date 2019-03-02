#ifndef RTSEIS_UTILITIES_INTERPOLATE_HPP
#define RTSEIS_UTILITIES_INTERPOLATE_HPP 1
#include <vector>

namespace RTSeis
{
namespace Utilities
{
namespace Math
{
/*!
 * @defgroup rtseis_utils_math_interp Interpolation
 * @brief Utilities for interpolation.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_math
 */
namespace Interpolate
{

/*!
 * @brief Methodologies for 1D interpolation.
 */
enum Method
{
    NEAREST,          /*!< Nearest neighbor interpolation.  At least
                           1 data point must be defined. */
    LINEAR,           /*!< Linear interpolation.  At least 2 points
                           must be defined. */
    CSPLINE_NATURAL,  /*!< Natural cubic spline.  At least 3 points
                           must be defined. */
    CSPLINE_PERIODIC, /*!< Cubic spline with periodic boundary conditions.
                           At least 3 points must be specified with the
                           further constraint that the first and last
                           specified function values are equal. */
    AKIMA,            /*!< Akima spline.  At least 3 points must be
                           defined. */
    AKIMA_PERIODIC    /*!< Akima spline with periodic boundary conditions.
                           At least 3 points must be specified with the
                           further constraint that the first and last
                           specified function values are equal. */   
};

/*!
 * @brief Interpolates a function $\f v = f(x) \f$ at interpolatoin points
 *        \f$ x_q \f$. 
 * @param[in] x       The x ordinates at which the function, \f$ v \f$, has
 *                    been evaluated.
 * @param[in] v       The function values \f$ v = f(x) \f$.  The size of v
 *                    must match the size of x.
 * @param[in] xq      The points at which to interpolate \f$ v \f$.
 * @param[out] vq     The interpolated function values at \f$ x_q \f$.
 * @param[in] method  The interpolation strategy.
 * @param[in] lxSorted   This is a data hint that indicates that the x ordinates
 *                       are in increasing monotonic order.
 * @param[in] lxqSorted  This is a data hint that indicates that the xq 
 *                       ordinates at which to interpolate are in increasing
 *                       monotonic order.
 * @throws std::invalid_argument for incorrect values which typically pertain
 *         to the requisite number of data and function points.
 */
void interp1d(const std::vector<double> &x,
              const std::vector<double> &v,
              const std::vector<double> &xq,
              std::vector<double> &vq,
              const Interpolate::Method method = Interpolate::Method::LINEAR,
              const bool lxSorted=false,
              const bool lxqSorted=false);

/*!
 * @brief FFT interpolation of a signal x.
 * @param[in] x   The signal to interpolate.
 * @param[in] npnew  This is the desired number of points in the interpolated
 *                   signal, yint.  This must be positive will likely be greater
 *                   than the length of x.  Otherwise, one should be sure to
 *                   lowpass filter x prior to calling this function to avoid
 *                   aliasing.
 * @param[out] yint  The interpolate version of x.
 * @throws std::invalid_argument if npnew is invalid or x is empty.
 */
void interpft(const std::vector<double> &x, const int npnew,
              std::vector<double> &yint);


}; // End Interpoalte
}; // End Math
}; // End Utilities
}; // End RTSeis
#endif
