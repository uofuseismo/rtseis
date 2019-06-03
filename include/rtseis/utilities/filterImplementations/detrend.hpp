#ifndef RTSEIS_UTILITIES_FILTER_DETREND_HPP
#define RTSEIS_UTILITIES_FILTER_DETREND_HPP 1
#include <memory>
#include "rtseis/utilities/filterImplementations/enums.hpp"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{
namespace FilterImplementations
{

/*!
 * @class Detrend detrend.hpp "include/rtseis/utilities/filterImplementations/detrend.hpp"
 * @brief Removes the mean or trend from a signal.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
*/
class Detrend
{
public:
    /*! @name Constructors
     * @{
     */ 
    /*!
     * @brief Default constructor.
     */
    Detrend();
    /*!
     * @brief Copy constructor.
     * @param[in] detrend  Detrend class from which to initialize this class.
     */
    Detrend(const Detrend &detrend);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] detrend  The detrend class to copy.
     * @result A deep copy of the detrend class.
     */
    Detrend& operator=(const Detrend &detrend);
    /*! @} */
     

    /*! @name Destructors
     * @{
     */
    ~Detrend();
    /*!
     * @brief Resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the class.
     * @param[in] type       Defines the trend removal strategy.
     * @param[in] precision  Defines the precision of the underlying
     *                       calculation.
     */
    void initialize(const DetrendType type,
                    const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Determines whether or not the class is inititalized.
     * @result True indicates that the class is inititalized.
     */
    bool isInitialized() const noexcept;

    /*!
     * @brief Detrends the data.
     * @param[in] nx  The number of samples in the signal.
     * @param[in] x   The signal to demean or detrend.  This is an array of
     *                dimension [nx].
     * @param[out] y  The demeaned or detrended data.  This is an arra of
     *                dimension [nx].
     * @throws std::invalid_argument if nx is positive and x or y is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @note If the detrend type is linear and nx is 1 then only the mean will
     *       be removed.
     */
    void apply(const int nx, const double x[], double *y[]);
    /*! @copydoc apply */
    void apply(const int nx, const float x[], float *y[]);
private:
    class DetrendImpl;
    std::unique_ptr<DetrendImpl> pImpl;
};

/*!
 * @brief Function that removes the mean from the data.
 * @param[in] nx      The number of samples in the signal.
 * @param[in] x       The signal from which to remove the mean.  This is
 *                    an array of dimension [nx].
 * @param[out] y      The demeaned version of x.  This is an array of
 *                    dimension [nx].
 * @param[out] mean   The mean of the dataset.
 * @throws std::invalid_argument if nx is positive and x or y is NULL.
 */
void removeMean(const int nx, const double x[], double *y[], double *mean);
/*! @copydoc removeMean */
void removeMean(const int nx, const float x[], float *y[], float *mean);

/*!
 * @brief Function that removes the linear trend from the data.
 * @param[in] nx          The number of samples in the signal.
 * @param[in] x           The signal from which to remove the trend.  This is
 *                        an array of dimension [nx].
 * @param[out] y          The detrended version of x.  This is an array of
 *                        dimension [nx].
 * @param[out] intercept  The intercept of the best-fitting line.
 * @param[out] slope      The slope of the best fitting line.
 * @throws std::invalid_argument if nx is positive and x or y is NULL.
 * @note If nx is less than 2 then the mean will be removed and the slope
 *       set to 0.
 */
void removeTrend(const int nx, const double x[], double *y[],
                 double *intercept, double *slope);
/*! @copydoc removeTrend */
void removeTrend(const int nx, const float x[], float *y[],
                 float *intercept, float *slope);

} // End filter implementations
} // End utilities
} // End rtseis
#endif

