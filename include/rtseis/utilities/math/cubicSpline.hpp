#ifndef RTSEIS_UTILITIES_MATH_CUBICSPLINE_HPP__
#define RTSEIS_UTILITIES_MATH_CUBICSPLINE_HPP__ 1
#include <memory>

namespace RTSeis::Utilities::Math::Interpolation
{
/*!
 * @brief Defines the types of boundary conditions available for use in the
 *        cubic spline interpolation.
 * @author Ben Baker (University of Utah)
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_math_interpolation
 */
enum class CubicSplineBoundaryConditionType
{
    NOT_A_KNOT, /*!< The first and second segment at a curve end are the
                     same.  This is useful when there is no information on
                     the boundary conditions. */
    NATURAL,    /*!< The second derivative at the curve ends are set to 0, i.e.,
                     \f$ y''[0] = y''[ny-1] \f$. */
    CLAMPED,    /*!< The first derivative at the curve ends are set to 0, i.e.,
                     \f$  y'[0] = y'[ny-1] \f$. */
    PERIODIC    /*!< The interpolated function is assumed to be periodic.
                     In this case, the first y value must equal the last
                     y value.  The resulting boundary condition will be
                     \f$  y'[0] = y'[ny-1] \f$ and
                     \f$ y''[0] = y''[ny-1] \f$. */
};

/*!
 * @brief A class for performing cubic spline interoplation.
 * @author Ben Baker (University of Utah)
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_math_interpolation
 */
class CubicSpline
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    CubicSpline();
    /*!
     * @brief Move constructor.
     * @param[in,out] spline  The cubic spline from which to initialize this
     *                        class.  On exit, spline's behavior is undefined.
     */
    CubicSpline(CubicSpline &&spline) noexcept;
    /*! @} */

    /*!
     * @brief Move assignment operator.
     * @param[in,out] spline  The cubic spline to move. 
     *                        On exit, spline's behavior is undefined.
     * @result Spline's memory moved onto this.
     */
    CubicSpline& operator=(CubicSpline &&spline) noexcept;

    /*! @name Destructor
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~CubicSpline();
    /*!
     * @brief Releases all memory on the module and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the cubic spline over the closed interval 
     *        [xInterval.first, xInterval.second] with function values
     *        \f$ y(x) \f$.
     * @note This will assume the y values are distributed uniformly over
     *       interval.
     * @param[in] npts       The number of data points in y.  This must be at
     *                       least 4.
     * @param[in] xInterval  The closed interval which begins at xPair.first
     *                       and ends at xPair.second.  xPair.first must be
     *                       less than xPair.second.
     * @param[in] y          The function values.  This is an array of dimension
     *                       [npts].  If the boundary conditions are periodic
     *                       then y[0] must equal y[npts-1].
     * @param[in] bcType     Defines the boundary conditions and by extension
     *                       the type of cubic spline.
     * @throws std::invalid_argument if any of the arguments are invalid.
     */ 
    void initialize(const int npts,
                    const std::pair<double, double> xInterval,
                    const double y[],
                    const CubicSplineBoundaryConditionType boundaryCondition);
    /*!
     * @brief Initializes the cubic spline for data points defined
     *        at \f$ y(x) \f$.
     * @param[in] npts     The number of data points in x and y.  This must be
     *                     at least 4.
     * @param[in] x        The ordinates at which y is evaluated.  This is an
     *                     array of dimension [npts] with the further caveat
     *                     that
     *                     \f$ x_i < x_{i+1} \forall i=0,\cdots,n_{pts}-2 \f$.
     * @param[in] y        The function values.  This is an array of dimension
     *                     [npts].  If the boundary conditions are periodic
     *                     then y[0] must equal y[npts-1]..
     * @param[in] bcType   Defines the boundary conditions and by extension
     *                     the type of cubic spline.
     * @throws std::invalid_argument if any of the arguments are invalid.
     */
    void initialize(const int npts,
                    const double x[],
                    const double y[],
                    const CubicSplineBoundaryConditionType boundaryCondition);
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
     *        is the cubic spline.
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
     * @sa isInitialized(), getMinimumX(), getMaximum().
     */
    void interpolate(const int nq, const double xq[], double *yq[]) const;
    /*!
     * @brief Integrates the cubic spline over the interval.
     * @param[in] interval  Defines the integration interval which goes
     * @throws std::runtime_error if the class was not intitialized.
     * @throws std::invalid_argument if the interval exceeds the interval
     *         on which the spline was constructed.
     * @sa isInitialized(), getMinimumX(), getMaximum().
     */
    double integrate(const std::pair<double, double> &interval) const;
private:
    class CubicSplineImpl;
    std::unique_ptr<CubicSplineImpl> pImpl;
    CubicSpline& operator=(const CubicSpline &spline) = delete;
}; 

}
#endif
