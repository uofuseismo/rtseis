#ifndef RTSEIS_UTILS_FR_SOS_HPP
#define RTSEIS_UTILS_FR_SOS_HPP 1
#include <cstdio>
#include <vector>
#include "rtseis/config.h"

namespace RTSeis
{
namespace Utilities
{
namespace FilterRepresentations
{
/*!
 * @class SOS sos.hpp "include/rtseis/utilities/filterRepresentations/sos.hpp"
 * @brief Represents a filter in terms of a series of cascaded
 *        second-order-sections.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_fr
 */
class SOS
{
    public:
        /*!
         * @brief Default constructor. 
         */
        SOS(void);
        /*!
         * @brief Copy constructor.
         * @param[in] sos  SOS class from which to initialize.
         */
        SOS(const SOS &sos);
        /*!
         * @brief Constructs an SOS class from the given second order sections.
         * @param[in] ns    Number of sections.
         * @param[in] bs    The numerator coefficients to set.  This 
         *                  is an array of dimension [3 x ns] with leading
         *                  dimension 3.  Moreover, the leading coefficient
         *                  of each second-order-section cannot be 0.
         * @param[in] as    The denoiminator coefficients to set.  This is
         *                  is an array of dimension [3 x ns] with leading
         *                  dimension 3.  Moreover, the leading coefficient
         *                  of each second-order-section cannot be 0.
         */
        SOS(const int ns,
            const std::vector<double> &bs,
            const std::vector<double> &as);
        /*!
         * @brief Assignement operator.
         * @param[in] sos  SOS class to copy.
         * @result A deep copy of the input class.
         */
        SOS &operator=(const SOS &sos);
        /*!
         * @brief Equality operator.
         * @param sos  Class to compare to this class.
         * @result True indicates that sos equals this class to within a given
         *         tolerance.
         */
        bool operator==(const SOS &sos) const;
        /*!
         * @brief Inequality operator.
         * @param[in] sos  Class to compare to this class.
         * @result True indicates that sos does not equal this class to within a
         *         given tolerance.
         */
        bool operator!=(const SOS &sos) const;
        /*!
         * @brief Default destructor.
         */
        ~SOS(void);
        /*!
         * @brief Prints the SOS structure.
         * @param[in] fout   File handle to print to.  If fout is NULL then this
         *                   will print to stdout.
         */
        void print(FILE *fout = stdout);
        /*!
         * @brief Clears the class.
         */
        void clear(void);
        /*!
         * @brief Sets the second order sections on the class.
         * @param[in] ns   The number of second order sections.
         * @param[in] bs   The numerator coefficients.  This has dimension
         *                 [3 x ns] with leading order 3.  Futhermore, the
         *                 \f$ b[3 i_s] \f$ cannot be 0 for
         *                 \f$ i=1,2, \cdots n_s \f$.
         * @param[in] as   The denominator coefficients.  This has dimension
         *                 [3 x ns] with leading order 3.  Futhermore,
         *                 the \f$ b[3 i_s] \f$ cannot be 0 for
         *                 \f$ i=1,2, \cdots n_s \f$.
         * @result 0 indicates success.
         */
        int setSecondOrderSections(const int ns,
                                   const std::vector<double> &bs,
                                   const std::vector<double> &as);
        /*!
         * @brief Returns the numerator coefficients of the SOS filter.
         * @brief A vector holding the numerator coefficients.  This has
         *        dimension [3 x ns] with leading dimension 3 where ns is
         *        given by getNumberOfSections().
         */
        std::vector<double> getNumeratorCoefficients(void) const;
        /*!
         * @brief Gets the denominator coefficients of the SOS filter.
         * @brief A vector holding the denominator coefficients.  This has
         *        dimension [3 x ns] with leading dimension 3 where ns is
         *        given by getNumberOfSections().
         */
        std::vector<double> getDenominatorCoefficients(void) const;
        /*!
         * @brief Returns the number of sections.
         * @result The number of second order sections.
         */
        int getNumberOfSections(void) const;
        /*!
         * @brief Sets the tolerance when testing for equality.
         * @param[in] tol   Tolerance.
         */
        void setEqualityTolerance(const double tol = 1.e-12);
    private:
        /*!< The numerator sections. */
        std::vector<double> bs_; 
        /*!< The denominator sections. */
        std::vector<double> as_; 
        /*!< The number of sections. */
        int ns_ = 0; 
        /*!< Default tolerance. */
        const double defaultTol_ = 1.e-12;
        /*!< Tolerance in checking equality. */
        double tol_ = defaultTol_;
}; // End SOS
}; // End FilterRepresentations
}; // End Utilities
}; // End RTSeis

#endif
