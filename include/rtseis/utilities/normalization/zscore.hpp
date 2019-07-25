#ifndef RTSEIS_UTILITIES_NORMALIZATION_ZSCORE_HPP
#define RTSEIS_UTILITIES_NORMALIZATION_ZSCORE_HPP 1
#include <memory>
namespace RTSeis::Utilities::Normalization
{
/*!
 * @class ZScore zscore.hpp "include/rtseis/utilities/normalization/zscore.hpp"
 * @brief Standardizes data by applying the transform
 *        \f[
 *            z = \frac{x - \mu}{\sigma}
 *        \f]
 *        where \f$ \mu$ \f$ is the mean and \f$ \sigma \f$ the standard
 *        deviation.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class ZScore
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    ZScore();
    /*!
     * @brief Copy constructor
     * @param[in] zscore  The z-score normalization class from to initialize
     *                    this class.
     */
    ZScore(const ZScore &zscore);
    /*!
     * @brief Move constructor
     * @param[in,out] zscore  The z-score normalization class whose memory will
     *                        be moved to this.  On exit zscore's behavior is
     *                        undefined.
     */
    ZScore(ZScore &&zscore) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] zscore   The zscore class to copy.
     * @result A deep copy of zscore.
     */ 
    ZScore& operator=(const ZScore &zscore);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] zscore  The zscore class to move.  On exit, zscore's
     *                        behavior will be undefined.
     * @result zscore's memory moved to this.
     */
    ZScore& operator=(ZScore &&zscore) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~ZScore();
    /*!
     * @brief Resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the class.
     * @param[in] mean     The mean of the data.
     * @param[in] stdDev   The (sample) standard deviation of the data.
     * @throws std::invalid_argument if the standard deviation is not positive.
     */
    void initialize(const double mean, const double stdDev);
    /*!
     * @brief Initializes the class by computing the mean and sample standard
     *        deviation.
     * @param[in] nx   The number of samples in x.  This must be at least 2.
     * @param[in] x    The array of which to compute the mean and standard
     *                 deviation.  This is an array of dimension [nx] and
     *                 cannot be comprised of uniform values.
     * @throws std::invalid_argument if nx is too small, x is NULL, or x is
     *         comprised of all identical values.
     */
    void initialize(const int nx, const double x[]);
    /*!
     * @brief Determines if the class is initialized.
     * @retval True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Applies the z-score normalization.
     * @param[in] nx     The number of points in the time series.
     * @param[in] x      The input time series to standardize.  This is an
     *                   array of dimension [nx].
     * @param[out] y     The standardized variant of x.  This is an array of
     *                   dimension [nx].
     * @throws std::invalid_argument if nx is positive and x or y is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void apply(const int npts, const double x[], double *y[]);
    /*!@copydoc apply */
    void apply(const int npts, const float x[], float *y[]);
private:
    class ZScoreImpl;
    std::unique_ptr<ZScoreImpl> pImpl; 
};
/*!
 * @brief Computes the mean and standard deviation of a vector.
 * @param[in] nx   The number of samples in x.
 * @param[in] x    The values of which to compute the mean and standard
 *                 deviation.  This is an array of dimension [nx].
 * @result A pair where pair.first contains the mean and pair.second
 *         contains the standard deviation.
 */
std::pair<double,double>
computeMeanAndStanardDeviation(const int nx, const double x[]);
}
#endif
