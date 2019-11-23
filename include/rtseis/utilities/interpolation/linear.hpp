#ifndef RTSEIS_UTILITIES_INTERPOLATION_LINEAR_HPP
#define RTSEIS_UTILITIES_INTERPOLATION_LINEAR_HPP 1
#include <memory>

namespace RTSeis::Utilities::Interpolation
{
/*!
 * @class Linear linear.hpp "rtseis/utilities/interpolation/linear.hpp"
 * @brief A class for performing linear interpolation.
 * @author Ben Baker (University of Utah)
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_math_interpolation
 */
class Linear
{
public:
    /*! @name Constructors
     * @{
     */
    /*! 
     * @brief Default constructor.
     */
    Linear();
    /*! 
     * @brief Move constructor.
     * @param[in,out] linear  The linear interpolation class from which to
     *                        initialize this class.  On exit, linear's behavior
     *                        is undefined.
     */
    Linear(CubicSpline &&spline) noexcept;
    /*! @} */

    /*!
     * @brief Move assignment operator.
     * @param[in,out] linear  The linear interoplation class to move.
     *                        On exit, linear's behavior is undefined.
     * @result Linear's memory onto this.
     */
    Linear& operator=(Linear &&linear) noexcept;

    /*! @name Destructor
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~Linear();
    /*!
     * @brief Releases all memory on the module and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the linear interpolation over the closed interval 
     *        [xInterval.first, xInterval.second] with regularly spaced
     *        function values \f$ y(x) \f$.
     * @param[in] npts       The number of data points in y.  This must be at
     *                       least 2.
     * @param[in] xInterval  The closed interval which begins at xInterval.first
     *                       and ends at xInterval.second.  xInterval.first must
     *                       be less than xInterval.second.
     * @param[in] y          The function values.  This is an array of dimension
     *                       [npts].
     * @throws std::invalid_argument if any of the arguments are invalid.
     */
    void initialize(int npts,
                    const std::pair<double, double> xInterval,
                    const double y[]);
    /*!
     * @brief Initializes the linear interpolation for data points defined
     *        at \f$ y(x) \f$.
     * @param[in] npts     The number of data points in x and y.  This must be
     *                     at least 2.
     * @param[in] x        The ordinates at which y is evaluated.  This is an
     *                     array of dimension [npts] with the further caveat
     *                     that
     *                     \f$ x_i < x_{i+1} \forall i=0,\cdots,n_{pts}-2 \f$.
     * @param[in] y        The function values.  This is an array of dimension
     *                     [npts].
     * @throws std::invalid_argument if any of the arguments are invalid.
     */
    void initialize(int npts,
                    const double x[],
                    const double y[]);

    /*! 
     * @brief Determines if the cubic spline has been inititalized or not.
     * @retval True indicates that the spline was inititalized.
     */
    bool isInitialized() const noexcept;
    /*! 
     * @brief Gets the minimum x ordinate that can be interpolated.
     * @retval The minimum x value that can be interpolated.
     * @throws std::runtime_error if the class is not initialized.
     * @sa isInitialized()
     */
    double getMinimumX() const;
    /*! 
     * @brief Gets the maximum x ordinate that can be interpolated.
     * @retval The maximum x value that can be interpolated.
     * @throws std::runtime_error if the class is not initialized.
     * @sa isInitialized()
     */
    double getMaximumX() const;

    /*!
     * @brief Interpolates function values \f$ y_q = f(x_q) \f$ where \f$ f \f$
     *        is the linear interpolating function.
     * @param[in] nq   The number of points at which to interpolate.
     * @param[in] xq   The abscissas at which to interpolate yq.  This is
     *                 an array of dimension [nq].  Additionally, each
     *                 xq must be in the range 
     *                 [\c getMinimumX(), \c getMaximumX()].
     * @param[out] yq  The interpolated values at \f$ x_q \f$.  This is an
     *                 array of dimension [nq].
     * @throws std::runtime_error if the class was not initialized.
     * @throws std::invalid_argument if xq or yq is NULL or any xq is
     *         out of the interpolation range.
     * @sa isInitialized(), getMinimumX(), getMaximum()
     */
    void interpolate(int nq, const double xq[], double *yq[]) const;
    /*!
     * @brief Interpolates the function faluves \f$ y_q = f(x_q) \f$ where
     *        \f$ f \f$ is the linear interpolation function.  This is for
     *        evenly spaced points on the interval [xq.first, xq.second].
     * @param[in] nq   The number of points at which to interpolate.
     * @param[in] xq   The interval over which to interpolate.  Note, that
     *                 xq.first must be less than xq.second.
     * @param[out] yq  The interpolated values at the evenly spaced points on
     *                 the interval given by \f$ x_q \f$.  This is an array
     *                 of dimension [nq].  Note, \f$ y_q[0] = f(x_q.first) \f$
     *                 and \f$ y_q[nq-1] = f(x_q.second \f$.
     * @throws std::runtime_error if the class was not initialized.
     * @throws std::invalid_argument if xq.first is less than \c getMinimumX()
     *         or xq.second is greater than \c getMaximumX()
     *         or xq.first >= xq.second.
     * @sa isInitialized(), getMinimumX(), getMaximum()
     */
    void interpolate(int nq,
                     const std::pair<double, double> xq,
                     double *yq[]) const;
private:
    class LinearImpl;
    std::unique_ptr<LinearImpl> pImpl;
    Linear& operator=(const Linear &linear) = delete;
};
}
#endif
