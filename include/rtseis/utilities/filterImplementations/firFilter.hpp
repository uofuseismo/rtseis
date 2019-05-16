#ifndef RTSEIS_UTILS_FILTER_FIR_HPP
#define RTSEIS_UTILS_FILTER_FIR_HPP 1
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
 * @class FIRFilter firFilter.hpp "include/rtseis/utilities/filterImplementations/firFilter.hpp"
 * @brief This is the core implementation for FIR filtering.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
class FIRFilter
{
public:
    /*!
     * @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    FIRFilter();
    /*!
     * @brief Copy constructor.
     * @param[in] fir   FIR class from which to initialize.
     */
    FIRFilter(const FIRFilter &fir);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy operator.
     * @param[in] fir   FIR class to copy.
     * @result A deep copy of the FIR class.
     */
    FIRFilter& operator=(const FIRFilter &fir);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~FIRFilter();
    /*! @} */

    /*!
     * @brief Initializes the FIR filter.
     * @param[in] nb    Number of numerator coefficients.
     * @param[in] b     Numerator coefficients.  This is an array of
     *                  dimension [nb].
     * @param[in] mode  The processing mode.  By default this
     *                  is for post-processing.
     * @param[in] precision   The precision of the filter.  By default
     *                        this is double precision.
     * @param[in] implementation  Defines the implementation.
     *                            The default is to use the direct form.
     * @result 0 indicates success.
     */
    int initialize(const int nb, const double b[],
                   const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                   const RTSeis::Precision precision = RTSeis::Precision::DOUBLE,
                   FIRImplementation implementation = FIRImplementation::DIRECT);
    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Utility routine to determine the initial condition length.
     * @result A non-negative number is the length of the initial
     *         condition array.
     */
    int getInitialConditionLength() const;
    /*!
     * @brief Sets the initial conditions for the filter.  This should
     *        be called prior to filter application as it will reset
     *        the filter.
     * @param[in] nz   The FIR filter initial condition length.
     *                 This should be equal to
     *                 getInitialConditionLength().
     * @param[in] zi   The initial conditions.  This has dimension [nz].
     * @result 0 indicates success.
     */
    int setInitialConditions(const int nz, const double zi[]);
    /*!
     * @brief Gets a copy of the initial conditions.
     * @param[in] nz   The FIR filter initial condition length.
     *                 This should be equal to
     *                 getInitialConditionLength().
     * @param[out] zi  The initial conditions.  This has dimension [nz].
     * @result 0 indicate success.
     */
    int getInitialConditions(const int nz, double zi[]) const;
    /*!
     * @{
     * @brief Applies the FIR filter to the data.
     * @param[in] n   Number of points in signals.
     * @param[in] x   Signal to filter.  This has dimension [n].
     * @param[out] y  The filtered signal.  This has dimension [n].
     */
    int apply(const int n, const double x[], double y[]);
    int apply(const int n, const float x[], float y[]);
    /*! @} */
    /*!
     * @brief Resets the initial conditions on the source delay line to
     *        the default initial conditions or the initial conditions
     *        set when FIRFilter::setInitialConditions() was called.
     * @result 0 indicates success.
     */
    int resetInitialConditions();
    /*!
     * @brief Clears the module and resets all parameters.
     */
    void clear() noexcept;
private:
    class FIRImpl;
    std::unique_ptr<FIRImpl> pFIR_;
}; // FIR
}; // filterimplementations
}; // utilties
}; // rtseis
#endif
