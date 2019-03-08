#ifndef RTSEIS_UTILS_FILTER_IIRIIR_HPP
#define RTSEIS_UTILS_FILTER_IIRIIR_HPP 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{
namespace FilterImplementations
{
/*!
 * @class IIRIIRFilter iiriirFilter.hpp "include/rtseis/utilities/filterImplementations/iiriirFilter.hpp"
 * @brief Implements a zero-phase IIR filter.  This is for
 *        post-processing only.
 * @ingroup rtseis_utils_filters
 * @copyright Ben Baker distributed under the MIT license.
 */
class IIRIIRFilter
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.  Note, this class is not yet 
     *        initialized and cannot be used.
     */
    IIRIIRFilter(void);
    /*! 
     * @brief A copy constructor.
     * @param[in] iiriir  The IIRIIRFilter class to initialize from.
     */
    IIRIIRFilter(const IIRIIRFilter &iiriir);
    /*! 
     * @brief A copy operator.
     * @param[in] iiriir  The IIRIIRFilter class to copy.
     * @result A deep copy of the iiriir filter class.
     */
    IIRIIRFilter &operator=(const IIRIIRFilter &iiriir);
    /*! @} */

    /*!
     * @brief Default destructor.
     */
    ~IIRIIRFilter(void);
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
     * @result 0 indicates success.
     */
    int initialize(const int nb, const double b[],
                   const int na, const double a[],
                   const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    bool isInitialized(void) const;
    /*!
     * @brief Gets the length of the initial conditions array.
     * @result On successful exit this is the length of the initial
     *         conditions array.  On failure this is negative.
     */
    int getInitialConditionLength(void) const;
    /*!
     * @brief Sets the initial conditions.
     * @param[in] nz   The length of the initial conditions.  This
     *                 should equal getInitialConditionLength().
     * @param[in] zi   The initial conditions to set.  This has
     *                 has dimension [nwork_]. 
     * @result 0 indicates success.
     */
    int setInitialConditions(const int nz, const double zi[]);
    /*!
     * @brief Applies the zero-phase IIR filter to the data.  Note,
     *        the class must be initialized prior to using this function.
     * @param[in] n   Number of points in signal.
     * @param[in] x   The signal to filter.  This has dimension [n].
     * @param[out] y  The zero-phase IIR filtered signal.  This has
     *                dimension [n].
     * @result 0 indicates success.
     */
    int apply(const int n, const double x[], double y[]);
    /*!
     * @brief Applies the zero-phase IIR filter to the data.  Note,
     *        the class must be initialized prior to using this function.
     * @param[in] n   Number of points in signal.
     * @param[in] x   The signal to filter.  This has dimension [n].
     * @param[out] y  The zero-phase IIR filtered signal.  This has
     *                dimension [n].
     * @result 0 indicates success.
     */
    int apply(const int n, const float x[], float y[]); 
    /*!
     * @brief Resets the initial conditions to those set in
     *        setInitialConditions().  Note, this will not do anything
     *        as the final conditions are never extracted from the
     *        IIR filter.
     * @result 0 indicates success.
     */
    int resetInitialConditions(void); 
    /*!
     * @brief Clears memory and resets the filter.  After applying
     *        this function the filter must be re-initialized prior
     *        to being applied to the data.
     */
    void clear(void);
    /*!
     * @brief Utilility function to get the filter order.
     * @result The filter order.
     */
    int getFilterOrder(void) const;
private:
    class IIRIIRImpl;
    std::unique_ptr<IIRIIRImpl> pIIRIIR_;
}; // iiriirfilter
}; // filterimplementations
}; // utilities
}; // rtseis
#endif
