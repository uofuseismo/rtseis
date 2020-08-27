#ifndef RTSEIS_UTILS_FILTERIMPLEMENTATIONS_IIRIIR_HPP
#define RTSEIS_UTILS_FILTERIMPLEMENTATIONS_IIRIIR_HPP 1
#include <memory>

namespace RTSeis::Utilities::FilterImplementations
{
/*!
 * @class IIRIIRFilter iiriirFilter.hpp "include/rtseis/utilities/filterImplementations/iiriirFilter.hpp"
 * @brief Implements a zero-phase IIR filter.  This is for
 *        post-processing only.
 * @ingroup rtseis_utils_filters
 * @copyright Ben Baker distributed under the MIT license.
 */
template<class T = double>
class IIRIIRFilter
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    IIRIIRFilter();
    /*! 
     * @brief A copy constructor.
     * @param[in] iiriir  The IIRIIRFilter class to initialize from.
     */
    IIRIIRFilter(const IIRIIRFilter &iiriir);
    /*!
     * @brief Move constructor.
     * @param[in,out] iiriir   The IIRIIR filter class from which to initialize
     *                         this class.  On exit, iiriir's behavior
     *                         is undefined.
     */
    IIRIIRFilter(IIRIIRFilter &&iiriir) noexcept;
    /*! 
     * @brief Copy assignment operator.
     * @param[in] iiriir  The IIRIIRFilter class to copy.
     * @result A deep copy of the iiriir filter class.
     */
    IIRIIRFilter &operator=(const IIRIIRFilter &iiriir);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] iiriir   The IIRIIIR filter class to move to this.
     *                         On exit, iiriir's behavior is undefined.
     * @return The memory from iiriir moved to this.
     */
    IIRIIRFilter &operator=(IIRIIRFilter &&iiriir) noexcept;
    /*! @} */

    /*!
     * @brief Default destructor.
     */
    ~IIRIIRFilter();
    /*!
     * @brief Initializes the zero-phase IIR filter.
     * @param[in] nb   The number of numerator coefficients.  This must
     *                 be positive.
     * @param[in] b    The numerator (feedforward) coefficients.  This
     *                 is an array of dimension [nb].
     * @param[in] na   The number of denominator coefficients.  This
     *                 must be positive.
     * @param[in] a    The denominator (feedback) coefficients.  This is
     *                 an array of dimension [na].
     * @param[in] precision  The precision of the filter.  By default
     *                       this is double precision.
     * @throws std::invalid_argument if any of the arguments are invalid.
     */
    void initialize(int nb, const double b[],
                    int na, const double a[]);
    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    [[nodiscard]] bool isInitialized() const noexcept;
    /*!
     * @brief Gets the length of the initial conditions array.
     * @result The length of the initial conditions array.
     * @throws std::runtime_error if the class is not initialized.
     */
    [[nodiscard]] int getInitialConditionLength() const;
    /*!
     * @brief Sets the initial conditions.
     * @param[in] nz   The length of the initial conditions.  This
     *                 should equal getInitialConditionLength().
     * @param[in] zi   The initial conditions to set.  This has
     *                 has dimension [nwork_]. 
     * @throws std::invalid_argument if nz is invalid or nz is positive
     *         and zi is NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void setInitialConditions(int nz, const double zi[]);
    /*!
     * @brief Applies the zero-phase IIR filter to the data.  Note,
     *        the class must be initialized prior to using this function.
     * @param[in] n   Number of points in signal.
     * @param[in] x   The signal to filter.  This has dimension [n].
     * @param[out] y  The zero-phase IIR filtered signal.  This has
     *                dimension [n].
     * @throws std::invalid_argument if n is positive and x or y is NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void apply(int n, const T x[], T *y[]);
    /*!
     * @brief Resets the initial conditions to those set in
     *        setInitialConditions().  Note, this will not do anything
     *        as the final conditions are never extracted from the
     *        IIR filter.
     * @throws std::runtime_error if the class is not initialized.
     */
    void resetInitialConditions(); 
    /*!
     * @brief Clears memory and resets the filter.  After applying
     *        this function the filter must be re-initialized prior
     *        to being applied to the data.
     */
    void clear() noexcept;
    /*!
     * @brief Utilility function to get the filter order.
     * @result The filter order.
     * @throws std::runtime_error if the class is not initialized.
     */
    [[nodiscard]] int getFilterOrder() const;
private:
    class IIRIIRImpl;
    std::unique_ptr<IIRIIRImpl> pIIRIIR_;
}; // iiriirfilter
} // rtseis
#endif
