#ifndef RTSEIS_UTILITIES_INTERPOLATION_WEIGHTEDAVERAGESLOPES_HPP
#define RTSEIS_UTILITIES_INTERPOLATION_WEIGHTEDAVERAGESLOPES_HPP 1
#include <memory>
namespace RTSeis::Utilities::Interpolation
{
/*!
 * @class WeightedAverageSlopes weightedAverageSlopes.hpp "rtseis/utilities/interpolation/weightedAverageSlopes.hpp" 
 * @brief Performs the weighted average slopes interpolation described in
 */
template<class T>
class WeightedAverageSlopes
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    WeightedAverageSlopes();
    /*!
     * @brief Copy constructor.
     * @param[in] slopes  The weighted average slopes from which to initialize
     *                    this class. 
     */  
    WeightedAverageSlopes(const WeightedAverageSlopes &slopes);
    /*!
     * @brief Move constructor.
     * @param[in,out] slopes  The weightd average slopes from which to
     *                        initialize this class.  On exit, slopes's
     *                        behavior will be undefined.
     */
    WeightedAverageSlopes(WeightedAverageSlopes &&slopes);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] slopes  The slopes class to copy.
     * @result A deep copy of the slopes class.
     */
    WeightedAverageSlopes& operator=(const WeightedAverageSlopes &slopes);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] slopes  The slopes class whose memory will be moved to
     *                        this.  On exit, slopes's behavior is undefined.
     * @result The memory from slopes moved to this. 
     */
    WeightedAverageSlopes& operator=(WeightedAverageSlopes &&slopes) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~WeightedAverageSlopes();
    /*!
     * @brief Clears all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the weighted average over the closed interval
     *        [xInterval.first, xInterval.second] with regularly spaced
     *        function values \f$ y(x) \f$.
     * @param[in] npts       The number of data points in y.  This must be at
     *                       least 4.
     * @param[in] xInterval  The closed interval which begins at xInterval.first
     *                       and ends at xInterval.second.  xInterval.first must
     *                       be less than xInterval.second.
     * @param[in] y          The function values.  This is an array of dimension
     *                       [npts].
     * @throws std::invalid_argument if any of the arguments are invalid.
     */
    void initialize(int npts,
                    const std::pair<T, T> xInterval,
                    const T y[]);
    /*!
     * @brief Determines if the class has been initialized yet.
     * @result True indicates that the class was inititalized.
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
     * @brief Interolates function values \f$ y_q = f(x_q) \f$ where \f$ f \f$
     *        is the weighted average slopes interpolator.
     * @param[in] nq   The number of points at which to interpolate.
     * @param[in] xq   The abscissas at which to interpolate yq.  This is
     *                 an array of dimension [nq].  Additionally, each 
     *                 xq must be in the range 
     *                 [\c getMinimumX(), \c getMaximumX()].
     * @param[out] yq  The interpolates values at \f$ x_q \f$.  This is an
     *                 array of dimension [nq].
     * @throws std::runtime_error if the class was not initialized.
     * @throws std::invalid_argument if xq or yq is NULL or any xq is
     *         out of the interpolation range.
     * @sa isInitialized(), getMinimumX(), getMaximum().
     */
    void interpolate(int nq, const T xq[], T *yq[]) const; 
private:
    class WeightedAverageSlopesImpl;
    std::unique_ptr<WeightedAverageSlopesImpl> pImpl;
};
}
#endif
