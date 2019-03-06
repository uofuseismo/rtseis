#ifndef RTSEIS_UTILS_FR_FIR_HPP
#define RTSEIS_UTILS_FR_FIR_HPP 1
#include <vector>
#include <memory>

namespace RTSeis
{
namespace Utilities
{
namespace FilterRepresentations
{
/*!
 * @class FIR fir.hpp "include/rtseis/utilities/filterRepresentations/fir.hpp"
 * @brief Represents a finite impulse response filter in terms of filter taps.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_fr
 */
class FIR
{
    public:
        /*!
         * @brief Default constructor.
         */
        FIR(void);
        /*!
         * @brief Constructs an FIR class from the given filter taps.
         * @param[in] pTaps   Filter taps to set.
         */
        explicit FIR(const std::vector<double> &pTaps);
        /*!
         * @brief Copy constructor.
         * @param[in] fir  Class from which to initialize this class.
         */
        FIR(const FIR &fir);
        /*!
         * @brief Assignment operator.
         * @param[in] fir   Class to copy.
         * @result A deep copy of the input class.
         */
        FIR &operator=(const FIR &fir);
        /*!
         * @brief Equality operator.
         * @param[in] fir  Class to compare to this class.
         * @result True indicates that fir equals this class to within a given
         *         tolerance.
         */
        bool operator==(const FIR &fir) const;
        /*!
         * @brief Inequality operator.
         * @param[in] fir   Class to compare to this class.
         * @result True indicates that fir does not equal this class to within a
         *         given tolerance. 
         */
        bool operator!=(const FIR &fir) const;
        /*!
         * @brief Default destructor.
         */
        ~FIR(void);
        /*!
         * @brief Prints the FIR structure.
         * @param[in] fout   File handle to print to.  If fout is NULL then this
         *                   will print to stdout.
         */
        void print(FILE *fout = stdout);
        /*!
         * @brief Clears the structure.
         */
        void clear(void);
        /*!
         * @brief Gets the number of filter coefficients.
         * @result Number of filter taps.
         */
        int getNumberOfFilterTaps(void) const;
        /*!
         * @{
         * @brief Sets the filter taps.
         * @param[in] n     The number of filter taps.
         * @param[in] taps  The filter coefficients to set.  This is an array
         *                  of dimension [n].
         */
        void setFilterTaps(const size_t n, double taps[]);
        /*!
         * @brief Sets the filter taps.
         * @param[in] taps  The filter coefficients to set.
         */
        void setFilterTaps(const std::vector<double> &taps);
        /*! @} */
        /*!
         * @brief Gets the filter taps.
         * @result The filter coefficients.
         */
        std::vector<double> getFilterTaps(void) const;
        /*!
         * @brief Sets the tolerance when testing for equality.
         * @param[in] tol   Tolerance.
         */
        void setEqualityTolerance(const double tol = 1.e-12);
    private:
        class FIRImpl;
        std::unique_ptr<FIRImpl> pImpl_;
};
};
};
};
#endif
