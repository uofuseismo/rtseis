#ifndef RTSEIS_UTILITIES_NORMALIZATION_MINMAX_HPP
#define RTSEIS_UTILITIES_NORMALIZATION_MINMAX_HPP 1
#include <memory>
namespace RTSeis::Utilities::Normalization
{
/*!
 * @class MinMax minMax.hpp "rtseis/utilities/normalization/minMax.hpp"
 * @brief Performs min/max normalization 
 *        \f[
 *            y = \frac{y - min}{max - min}
 *        \f]
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class MinMax
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    MinMax();
    /*!
     * @brief Copy constructor.
     * @param[in] minMax  The min-max class from which to initialize this class.
     */
    MinMax(const MinMax &minMax);
    /*!
     * @brief Move constructor.
     * @param[in,out] minMax  The min-max class from which to initialize this
     *                        class.  On exit, minMax's behavior is undefined.
     */
    MinMax(MinMax &&minMax) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment oeprator
     * @param[in] minMax  The min-max class to copy.
     * @result A deep copy of the minMax class.
     */
    MinMax& operator=(const MinMax &minMax);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] minMax  The min-max class to move to this.
     *                        On exit minMax's behavior is undefined.
     * @result The memory of minMax moved to this.
     */
    MinMax& operator=(MinMax &&minMax) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Move assignment operator.
     *
     */
    ~MinMax();
    /*!
     * @brief Resets the class and releases memory.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the min-max rescaling.
     * @param[in] dataRange    Defines the data range where dataRange.first is
     *                         minimum and dataRange.second is the maximum.
     * @param[in] targetRange  Defines the target data range where 
     *                         targetRange.first is the minimum and
     *                         targetRange.second is the maximum.
     *                         By default this will rescale to [0,1].
     * @throws std::invalid_argument if dataRange.first == dataRange.second.
     */
    void initialize(std::pair<double, double> dataRange,
                    std::pair<double, double> targetRange = std::make_pair<double, double> (0.0, 1.0));
    /*!
     * @brief Initializes the min-max rescaling.  This will compute the data
     *        range from the given array.  
     * @param[in] npts         The number of data points in x.
     * @param[in] x            The data from which to extract the min and max.
     *                         This is an array of dimension [npts].
     * @param[in] targetRange  Defines the target data range where
     *                         targetRange.first is the minimum and
     *                         targetRange.second is the maximum.
     *                         By default this will rescale to [0,1].
     * @throws std::invalid_argument if npts is not positive, x is NULL, or all 
     *         of the samples in x have the value.
     */
    void initialize(const int npts, const double x[],
                    std::pair<double, double> targetRange = std::make_pair<double, double> (0.0, 1.0));
    /*!
     * @brief Determines if the class is initialized.
     * @retval True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;

    /*!
     * @brief Applies the min-max normalization.
     * @param[in] npts   The number of data points in x and y.
     * @param[in] x      The data to normalize.  This is an array of 
     *                   dimension [npts].
     * @param[out] y     The min-max normalized data.  This is an array of
     *                   dimension [npts].
     * @throws std::invalid_argument if npts is positive and x or y is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * \c \isInitialized()
     */
    void apply(const int npts, const double x[], double *y[]); 
    /*! @copydoc apply */
    void apply(const int npts, const float x[], float *y[]); 
private:
    class MinMaxImpl;
    std::unique_ptr<MinMaxImpl> pImpl;
};
}
#endif
