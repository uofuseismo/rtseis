#ifndef RTSEIS_UTILS_FILTER_MEDIAN_HPP
#define RTSEIS_UTILS_FILTER_MEDIAN_HPP 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{
namespace FilterImplementations
{
/*!
 * @class MedianFilter medianFilter.hpp "include/rtseis/utilities/filterImplementations/medianFilter.hpp"
 * @brief This is the core implementation for median filtering.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
class MedianFilter
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Copy constructor.
     */
    MedianFilter(void);
    /*!
     * @brief Copy constructor.
     * @param[in] median  Median class from which to initialize.
     */
    MedianFilter(const MedianFilter &median);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy operator.
     * @param[in] median  The median class to copy.
     * @result A deep copy of the median filter class.
     */
    MedianFilter& operator=(const MedianFilter &median);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*! 
     * @brief Destructor.
     */
    ~MedianFilter(void);
    /* @} */

    /*!
     * @brief Initializes the median filter.
     * @param[in] n   The window size of the median filter.  This must
     *                be a positive and odd number.  If n is not odd
     *                then it's length will be increased by 1.
     * @param[in] mode  The processing mode.  By default this
     *                  is for post-processing.
     * @param[in] precision   The precision of the filter.  By default
     *                        this is double precision.
     * @result 0 indicates success.
     */
    int initialize(const int n,
                   const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                   const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    bool isInitialized(void) const;
    /*!
     * @brief Utility routine to determine the initial condition length.
     * @retval A non-negative number is the length of the initial
     *         condition array.
     * @retval -1 Indicates failure.
     */
    int getInitialConditionLength(void) const;
    /*!
     * @brief Returns the group delay of the filter.  Note, that this
     *        shift is required to get a correspondence to Matlab.
     * @result The group delay.
     */
    int getGroupDelay(void) const;
    /*!
     * @brief Sets the initial conditions for the filter.  This should
     *        be called prior to filter application as it will reset
     *        the filter.
     * @param[in] nz   The median filter initial conditions.  This
     *                 should be equal to getInitialConditionLength().
     * @param[in] zi   The initial conditions.  This has dimension [nz].
     * @result 0 indicates success.
     */
    int setInitialConditions(const int nz, const double zi[]);
    /*!
     * @{
     * @brief Appplies the median filter to the array x.
     * @param[in] n   Number of points in x.
     * @param[in] x   The signal to filter.  This has dimension [n].
     * @param[out] y  The filtered signal.  This has dimension [n].
     * @result 0 indicates success.
     */
    int apply(const int n, const double x[], double y[]);
    int apply(const int n, const float x[], float y[]);
    /*! @} */
    /*!
     * @brief Resets the initial conditions on the source delay line to
     *        the default initial conditions or the initial conditions
     *        set when MedianFilter::setInitialConditions() was called.
     * @result 0 indicates success.
     */ 
    int resetInitialConditions(void);
    /*! 
     * @brief Clears the module and resets all parameters.
     */
    void clear(void);
 private:
    class MedianFilterImpl;
    std::unique_ptr<MedianFilterImpl> pMedian_;
}; // End medianFilter 
}; // End filter implementations
}; // End utilties
}; // End RTSeis
#endif
