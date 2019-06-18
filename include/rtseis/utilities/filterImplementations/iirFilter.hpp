#ifndef RTSEIS_UTILS_FILTER_IIR_HPP
#define RTSEIS_UTILS_FILTER_IIR_HPP 1
#include <memory>
#include "rtseis/enums.h"
#include "rtseis/utilities/filterImplementations/enums.hpp"

namespace RTSeis
{
namespace Utilities
{
namespace FilterImplementations
{

/*!
 * @class IIRFilter iirFilter.hpp "include/rtseis/utilities/filterImplementations/iirFilter.hpp"
 * @brief This is the core implementation for IIR filtering.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
class IIRFilter
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    IIRFilter();
    /*!
     * @brief Copy constructor.
     * @param[in] iir  IIR filter class from which to initialize this class.
     */
    IIRFilter(const IIRFilter &iir);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*! 
     * @brief Copy assignent operator.
     * @param[in] iir  IIR filter class to copy.
     * @result A deep copy of the IIR filter class.
     */
    IIRFilter &operator=(const IIRFilter &iir);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~IIRFilter();
    /*! @} */

    /*!
     * @brief Initializes the IIR filter.
     * @param[in] nb     Number of numerator coefficients.
     * @param[in] b      The numerator coefficients.  This has
     *                   dimension [nb].
     * @param[in] na     Number of denominator coefficients.
     * @param[in] a      The denominator coefficients.  This has
     *                   dimension [na].
     * @param[in] mode   Defines whether the filter is for real-time
     *                   or post-processing.
     * @param[in] precision  Defines the precision of the filter
     *                        implementation.
     * @param[in] implementation  Defines the algorithmic filter
     *                            implementation.
     * @throws std::invalid_argument if any arguments are invalid.
     */
    void initialize(const int nb, const double b[],
                    const int na, const double a[],
                    const RTSeis::ProcessingMode mode =  RTSeis::ProcessingMode::POST_PROCESSING,
                    const RTSeis::Precision precision = RTSeis::Precision::DOUBLE,
                    const IIRDFImplementation implementation = IIRDFImplementation::DF2_FAST);
    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Gets the length of the initial conditions array.
     * @result The length of the initial conditions array.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getInitialConditionLength() const;
    /*!
     * @brief Sets the initial conditions for the filter.  This should
     *        be called prior to filter application as it will reset
     *        the filter.
     * @param[in] nz   The IIR filter initial condition length.
     *                 This should be equal to
     *                 getInitialConditionLength().
     * @param[in] zi   The initial conditions.  This has dimension [nz].
     * @throws std::invalid_argument if nz is invalid or nz is positive
     *         zi is NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void setInitialConditions(const int nz, const double zi[]);
    /*!
     * @brief Applies the IIR filter.  Note, the class must be
     *        initialized prior to using this function.
     * @param[in] n   The number of points in the signal.
     * @param[in] x   The input signal to filter.  This has dimension [n].
     * @param[out] y  The filtered signal.  This has dimension [n].
     * @throws std;:invalid_argument if n is positive and x or y is NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void apply(const int n, const double x[], double *y[]);
    /*! @copydoc apply() */
    void apply(const int n, const float x[], float *y[]);
    /*!
     * @brief Resets the initial conditions to those set in
     *        \c setInitialConditions().
     * @throws std::runtime_error if the class is not initialized.
     */
    void resetInitialConditions();
    /*!
     * @brief Clears memory and resets the filter.  After applying
     *        this function the filter must be re-initialized prior
     *        to being applied to the data.
     */
    void clear() noexcept;
private:
    class IIRFilterImpl;
    std::unique_ptr<IIRFilterImpl> pIIR_;
}; // IIR
}; // filterimplementations
}; // utilties
}; // rtseis
#endif
