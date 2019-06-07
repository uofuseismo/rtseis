#ifndef RTSEIS_UTILITIES_INTERPOLATE_HPP
#define RTSEIS_UTILITIES_INTERPOLATE_HPP 1
#include <vector>
#include <memory>

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

class Interp1D
{
/*! 
 * @brief Methodologies for 1D interpolation.
 */
public:
    enum Method
    {   
        NEAREST,          /*!< Nearest neighbor interpolation.  At least
                               1 data point must be defined. */
        LINEAR,           /*!< Linear interpolation.  At least 2 points
                               must be defined. */
        CSPLINE_NATURAL,  /*!< Natural cubic spline.  At least 3 points
                               must be defined. */
        CSPLINE_PERIODIC, /*!< Cubic spline with periodic boundary
                               conditoins.  At least 3 points must be
                               specified with the further constraint that
                               the first and last specified function values
                               are equal. */
        AKIMA,            /*!< Akima spline.  At least 3 points must be
                               defined. */
        AKIMA_PERIODIC  /*!< Akima spline with periodic boundary conditions.
                             At least 3 points must be specified with the
                             further constraint that the first and last
                             specified function values are equal. */
    };
public:
    /*!
     * @brief Default constructor.
     */
    Interp1D(void);
    /*!
     * @brief Default destructor.
     */
    ~Interp1D(void);
    /*!
     * @brief Creates an interpolator given the coordinates (x,v(x)).
     * @param[in] x   The x ordinates at which the function, \f$ v \f$, has
     *                been evaluated.  It is required that x be sorted in
     *                increasing order.
     * @param[in] v   The function values \f$ v = f(x) \f$.  The size of v
     *                must match the size of x.  For periodic boundary
     *                conditions v[0] must equal v[v.size()-1].
     * @param[in] method  The interpolation method.
     * @throws std::invalid_argument for invalid inputs.
     */
    void initialize(const std::vector<double> &x,
                    const std::vector<double> &v,
                    const Method method);
    /*!
     * @brief Creates an interpolator on the uniform grid (x,v(x)).
     *        whose domain is \f$ x \in [x.first(), x.second()] \f$.
     * @param[in] npts  Number of points in partition.
     * @param[in] x     The end points of the domain where x.first must
     *                  be less than x.second.
     * @param[in] v     The function values \f$ v = f(x) \f$.  The size of
     *                  v must match the size of npts.  For periodic
     *                  boundary conditions v[0] must equal v[npts-1].
     * @param[in] method  The interpolation method.
     * @throws std::invalid_argument for invalid inputs.
     */
    void initialize(const int npts,
                    const std::pair<double,double> x,
                    const std::vector<double> &v, 
                    const Interp1D::Method method);
    /*!
     * @brief Interpolates points \f$ \hat{y} = v_q(x_q) \f$.
     * @param[in] xq   The ordinates at which to interpoalte \f$ v_q \f$.
     * @param[out] vq  The interpolated function values at \f$ x_q \f$.
     *                 Note that vq.size() must equal xq.size().
     * @throws std::invalid_argument if xq.size() != vq.size() or the
     *         interpolator has yet to be initialized. 
     */
    void apply(const std::vector<double> &xq,
               std::vector<double> &vq);
    /*!
     * @brief Interpolates points \f$ \hat{y} = v_q(x_q) \f$.
     * @param[in] xq   The ordinates at which to interpoalte \f$ v_q \f$.
     * @param[out] vq  The interpolated function values at \f$ x_q \f$.
     *                 Note that vq.size() must equal xq.size().
     * @param[in] lxqSorted If true, then this indicates that the
     *            interpolation points, \f$ x_q \f$, are sorted in
     *            non-decreasing order.
     * @throws std::invalid_argument if xq.size() != vq.size() or the
     *         interpolator has yet to be initialized. 
     */
    void apply(const std::vector<double> &xq,
               std::vector<double> &vq,
               bool lxqSorted);
private:
    class Interp1DImpl;
    std::unique_ptr<Interp1DImpl> pImpl;
}; // End Interp1D


/*!
 * @brief FFT interpolation of a signal x.
 * @param[in] x   The signal to interpolate.
 * @param[in] npnew  This is the desired number of points in the
 *                   interpolated signal, yint.  This must be positive
 *                   will likely be greater than the length of x.
 *                   Otherwise, one should be sure to lowpass filter x
 *                   prior to calling this function to avoid aliasing.
 * @result The interpolated variant of x.  This will have 
 *         dimension [npnew].
 * @throws std::invalid_argument if npnew is invalid or x is empty.
 * @throws std::runtime_error if an internal error has occurred. 
 */
std::vector<double> interpft(const std::vector<double> &x, const int npnew);


}; // End Interpoalte
}; // End Math
}; // End Utilities
}; // End RTSeis
#endif
