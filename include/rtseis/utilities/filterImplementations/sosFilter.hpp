#ifndef RTSEIS_UTILITIES_FILTER_SOS_HPP
#define RTSEIS_UTILITIES_FILTER_SOS_HPP 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{
namespace FilterImplementations
{

/*!
 * @class SOSFilter sosFilter.hpp "include/rtseis/utilities/filterImplementations/sosFilter.hpp"
 * @brief This is the core implementation for second order section (biquad)
 *        infinite impulse response filtering.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
*/
class SOSFilter
{
public:
    /*!
     * @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    SOSFilter();
    /*!
     * @brief A copy constructor.
     * @param[in] sos  The SOS class from which to initialize.
     */
    SOSFilter(const SOSFilter &sos);
    /*!
     * @brief Copy operator.
     * @param[in] sos  The class to copy.
     * @result A deep copy of the input SOS class.
     */
    SOSFilter& operator=(const SOSFilter &sos);
    /*! @} */

    /*!
     * @brief Default destructor.
     */
    ~SOSFilter();

    /*!
     * @brief Initializes the second order section filter.
     * @param[in] ns    The number of second order sections.
     * @param[in] bs    Numerator coefficients.  This is an array of
     *                  dimension [3 x ns] with leading dimension 3.
     *                  There is a further requirement that b[3*is]
     *                  for \f$ i_s=0,1,\cdots,n_s-1 \f$ not be zero.
     * @param[in] as    Denominator coefficients.  This is an array of
     *                  dimension [3 x ns] with leading dimension 3. 
     *                  There is a further requirement that a[3*is]
     *                  for \f$ i_s=0,1,\cdots,n_s-1 \f$ not be zero.
     * @param[in] mode  The processing mode.  By default this
     *                  is for post-processing.
     * @param[in] precision   The precision of the filter.  By default
     *                        this is double precision.
     * @result 0 indicates success.
     * @throws std::invalid_argument if ns, bs, or as is invalid.
     */
    void initialize(const int ns,
                    const double bs[],
                    const double as[],
                    const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                    const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Returns the length of the initial conditions.
     * @result The length of the initial condtions array.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getInitialConditionLength() const;
    /*!
     * @brief Sets the initial conditions for the filter.  This should
     *        be called prior to filter application as it will reset
     *        the filter.
     * @param[in] nz   The second order section filter initial
     *                 conditions.  This should be equal to
     *                 getInitialConditionLength().
     * @param[in] zi   The initial conditions.  This has dimension [nz].
     * @throws std::invalid_argument if nz is invalid or if nz is positive
     *         and zi is NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void setInitialConditions(const int nz, const double zi[]);

    /*! @name Filter Application
     * @{
     */
    /*!
     * @brief Applies the second order section filter to the data.
     * @param[in] n   Number of points in signals.
     * @param[in] x   The signal to filter.  This has dimension [n].
     * @param[out] y  The filtered signal.  This has dimension [n].
     * @throws std::invalid_argument if n is positive and x or y is NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void apply(const int n, const double x[], double *y[]);
    /*! @copydoc apply */
    void apply(const int n, const float x[], float *y[]);
    /*! @} */

    /*!
     * @brief Resets the initial conditions on the source delay line
     *        to the default initial conditions or the initial
     *        conditions set when SOSFilter::setInitialConditions()
     *        was called.
     * @throws std::runtime_error if the class is not initialized.
     */
    void resetInitialConditions();
    /*! 
     * @brief Clears the module and resets all parameters.
     */

    void clear() noexcept;
    /*!
     * @brief Gets the number of second order sections in the filter.
     * @result The number of second order sections.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getNumberOfSections() const;

private:
    class SOSFilterImpl;
    std::unique_ptr<SOSFilterImpl> pSOS_;
}; // sos
} // filterimplementations
} // utilties
} // rtseis

#endif
