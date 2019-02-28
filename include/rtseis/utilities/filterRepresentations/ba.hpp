#ifndef RTSEIS_UTILS_FR_BA_HPP
#define RTSEIS_UTILS_FR_BA_HPP 1
#include <cstdio>
#include <vector>
#include <memory>

namespace RTSeis
{
namespace Utilities
{
/*!
 * @defgroup rtseis_utils_fr Filter Representations
 * @brief Different representations of filters.  
 * @ingroup rtseis_utils
 */
namespace FilterRepresentations
{

/*!
 * @class BA ba.hpp "include/rtseis/utilities/filterRepresentations/ba.hpp"
 * @brief Represents an infinite impulse resonse filter in terms of
 *        numerator and denominator coefficients.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_fr
 */
class BA
{
    public:
        /*!
         * @brief Default constructor.
         */ 
        BA(void);
        /*!
         * @brief Constructs a transfer function class from the given numerator
         *        and denminator coefficients.
         * @param[in] b     The numerator coefficients.
         * @param[in] a     The denomiantor coefficients.
         */
        BA(const std::vector<double> &b, const std::vector<double> &a);
        /*!
         * @brief Copy constructor.
         * @param[in] ba  Class from which to initialized.
         */
        BA(const BA &ba);
        /*!
         * @brief Assignment operator.
         * @param[in] ba   Class to copy.
         * @result A deep copy of the input class.
         */
        BA &operator=(const BA &ba);
        /*!
         * @brief Equality operator.
         * @param[in] ba  Class to compare to this class.
         * @result True indicates that ba equals this class to within a given
         *         tolerance.
         */
        bool operator==(const BA &ba) const;
        /*!
         * @brief Inequality operator.
         * @param[in] ba   Class to compare to this class.
         * @result True indicates that ba does not equal this class to within a
         *         given tolerance. 
         */
        bool operator!=(const BA &ba) const;
        /*!
         * @brief Default destructor.
         */
        ~BA(void);
        /*!
         * @brief Prints the BA structure.
         * @param[in] fout   File handle to print to.  If fout is NULL then this
         *                   will print to stdout.
         */
        void print(FILE *fout = stdout);
        /*!
         * @brief Clears the structure.
         */
        void clear(void);
        /*!
         * @brief Gets the number of numerator coefficients.
         * @result Number of numerator coefficients in transfer function.
         */
        int getNumberOfNumeratorCoefficients(void) const;
        /*!
         * @brief Gets the number of denominator coefficients.
         * @result The number of denominator coefficients in the
         *         transfer function.
         */
        int getNumberOfDenominatorCoefficients(void) const;
        /*!
         * @{
         * @brief Sets the numerator coefficients of the transfer function.
         * @param[in] n    The number of numerator coefficients.
         * @param[in] b    The numerator coefficients to set.  This is an array
         *                 of dimension [n].
         */
        void setNumeratorCoefficients(const size_t n, double b[]);
        /*!
         * @brief Sets the numerator coefficients of the transfer function.
         * @param[in] b  Numerator coefficients to set on the transfer function.
         */
        void setNumeratorCoefficients(const std::vector<double> &b);
        /*! @} */
        /*!
         * @{
         * @brief Sets the denominator coefficients of the transfer function.
         * @param[in] n    The number of numerator coefficients.
         * @param[in] a    The denominator coefficients to set.  This is an
         *                 array of dimension [n].
         */
        void setDenominatorCoefficients(const size_t n, double a[]);
        /*! 
         * @brief Sets the denominator coefficients of the transfer function.
         * @param[in] a  The denominator coefficients to set on the transfer
         *               function.
         */
        void setDenominatorCoefficients(const std::vector<double> &a);
        /*! @} */
        /*!
         * @brief Gets the numerator coefficients.
         * @result The numerator coefficients.
         */
        std::vector<double> getNumeratorCoefficients(void) const;
        /*!
         * @brief Gets the denominator coefficients.
         * @result The denominator coefficients.
         */
        std::vector<double> getDenominatorCoefficients(void) const;
        /*!
         * @brief Sets the equality tolerance.
         */
        void setEqualityTolerance(const double tol = 1.e-12);
    private:
        class BAImpl;
        std::unique_ptr<BAImpl> pImpl_;
}; // End BA
}; // End FilterRepresentations
}; // End Utilities
}; // End RTSeis

#endif
