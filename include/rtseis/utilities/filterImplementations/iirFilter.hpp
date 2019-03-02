#ifndef RTSEIS_UTILS_FILTER_IIR_HPP
#define RTSEIS_UTILS_FILTER_IIR_HPP 1
#include <memory>
#include "rtseis/enums.h"

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
        /*!
         * @brief Defines the direct-form implementation.
         */
        enum Implementation
        {
            DF2_FAST, /*!< A fast IPP implementation that appears 
                           to only retain 3 digits of accuracy when
                            compared to Matlab. */
            DF2_SLOW  /*!< A slow implementation that is consistent
                           with Matlab. */
        };
    public:
        /*!
         * @brief Default destructor.
         */
        IIRFilter(void);
        /*!
         * @brief Copy constructor.
         * @param[in] iir  IIR filter class from which to initialize
         *                 this clsas.
         */
        IIRFilter(const IIRFilter &iir);
        /*! 
         * @brief Copy assignent operator.
         * @param[in] iir  IIR filter class to copy.
         * @result A deep copy of the IIR filter class.
         */
        IIRFilter &operator=(const IIRFilter &iir);
        /*!
         * @brief Default destructor.
         */
        ~IIRFilter(void);
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
         */
        int initialize(const int nb, const double b[],
                       const int na, const double a[],
                       const RTSeis::ProcessingMode mode =  RTSeis::ProcessingMode::POST_PROCESSING,
                       const RTSeis::Precision precision = RTSeis::Precision::DOUBLE,
                       const Implementation implementation = Implementation::DF2_FAST);
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
         * @brief Sets the initial conditions for the filter.  This should
         *        be called prior to filter application as it will reset
         *        the filter.
         * @param[in] nz   The IIR filter initial condition length.
         *                 This should be equal to
         *                 getInitialConditionLength().
         * @param[in] zi   The initial conditions.  This has dimension [nz].
         * @result 0 indicates success.
         */
        int setInitialConditions(const int nz, const double zi[]);
        /*!
         * @brief Applies the IIR filter.  Note, the class must be
         *        initialized prior to using this function.
         * @param[in] n   The number of points in the signal.
         * @param[in] x   The input signal to filter.  This has dimension [n].
         * @param[out] y  The filtered signal.  This has dimension [n].
         * @result 0 indicates success.
         */
        int apply(const int n, const double x[], double y[]);
        /*! @copydoc apply */
        int apply(const int n, const float x[], float y[]);
        /*!
         * @brief Resets the initial conditions to those set in
         *        setInitialConditions().
         * @result 0 indicates success.
         */
        int resetInitialConditions(void);
        /*!
         * @brief Clears memory and resets the filter.  After applying
         *        this function the filter must be re-initialized prior
         *        to being applied to the data.
         */
        void clear(void);
    private:
        class IIRFilterImpl;
        std::unique_ptr<IIRFilterImpl> pIIR_;
}; // IIR
}; // filterimplementations
}; // utilties
}; // rtseis
#endif
