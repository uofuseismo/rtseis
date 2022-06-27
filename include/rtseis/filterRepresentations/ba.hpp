#ifndef RTSEIS_FILTERREPRESENTATIONS_BA_HPP
#define RTSEIS_FILTERREPRESENTATIONS_BA_HPP 1
#include <vector>
#include <ostream>
#include <cstdio>
#include <ostream>
#include <memory>
namespace RTSeis::FilterRepresentations
{
/// @class BA ba.hpp "include/rtseis/filterRepresentations/ba.hpp"
/// @brief Represents an infinite impulse resonse filter in terms of
///        numerator and denominator coefficients.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_utils_fr
class BA
{
public:
    /// @name Constructors
    /// @{

    /// @brief Default constructor.
    BA();
    /// @brief Constructs a transfer function class from the given numerator
    ///        and denminator coefficients.
    /// @param[in] b     The numerator coefficients.
    /// @param[in] a     The denomiantor coefficients.
    /// @throws std::invalid_argument if b is empty, if a is empty, or
    ///         if a[0] is 0.
    BA(const std::vector<double> &b, const std::vector<double> &a);
    /// @brief Copy constructor.
    /// @param[in] ba  Class from which to initialized.
    BA(const BA &ba);
    /// @brief Move constructor.
    /// @param[in,out] ba   Class to move to this class.  On exit ba's behavior
    ///                     is undefined.
    BA(BA &&ba) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Assignment operator.
    /// @param[in] ba   Class to copy.
    /// @result A deep copy of the input class.
    BA &operator=(const BA &ba);
    /// @brief Move assignment operator.
    /// @param[in,out] ba  BA class to move.  On exit ba's behavior is undefined.
    /// @result The moved version of ba.
    BA &operator=(BA &&ba) noexcept;
    /// @brief Equality operator.
    /// @param[in] ba  Class to compare to this class.
    /// @result True indicates that ba equals this class to within a given
    ///         tolerance.
    [[nodiscard]] bool operator==(const BA &ba) const;
    /// @brief Inequality operator.
    /// @param[in] ba   Class to compare to this class.
    /// @result True indicates that ba does not equal this class to within a
    ///         given tolerance. 
    [[nodiscard]] bool operator!=(const BA &ba) const;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Default destructor.
    ~BA();
    /// @brief Clears the structure.
    void clear() noexcept;
    /// @}

    /// @brief Prints the BA structure.
    /// @param[in] fout   File handle to print to.  If fout is NULL then this
    ///                   will print to stdout.
    void print(FILE *fout = stdout) const noexcept;

    /// @result The number of numerator coefficients of transfer function.
    [[nodiscard]] int getNumberOfNumeratorCoefficients() const noexcept;
    /// @result The number of denominator coefficients of the transfer function.
    [[nodiscard]] int getNumberOfDenominatorCoefficients() const noexcept;
    /// @{
    /// @brief Sets the numerator coefficients of the transfer function.
    /// @param[in] n    The number of numerator coefficients.  This must be
    ///                 positive.
    /// @param[in] b    The numerator coefficients to set.  This is an array
    ///                 of dimension [n].
    /// @throws std::invalid_argument if n is less than 1 or b is NULL.
    void setNumeratorCoefficients(size_t n, const double b[]);
    /// @brief Sets the numerator coefficients of the transfer function.
    /// @param[in] b  Numerator coefficients to set on the transfer function.
    /// @throws std::invalid_argument if b is empty.
    void setNumeratorCoefficients(const std::vector<double> &b);
    /// @}

    /// @{
    /// @brief Sets the denominator coefficients of the transfer function.
    /// @param[in] n    The number of numerator coefficients.  This must be
    ///                 positive.
    /// @param[in] a    The denominator coefficients to set.  This is an
    ///                 array of dimension [n].
    /// @throws std::invalid_argument if a[0] is 0 or n < 1.
    void setDenominatorCoefficients(size_t n, const double a[]);
    /// @brief Sets the denominator coefficients of the transfer function.
    /// @param[in] a  The denominator coefficients to set on the transfer
    ///               function.
    /// @throws std::invalid_argument if a is empty or a[0] is 0.
    void setDenominatorCoefficients(const std::vector<double> &a);
    /// @}

    /// @result The numerator coefficients.
    [[nodiscard]] std::vector<double> getNumeratorCoefficients() const noexcept;
    /// @result The denominator coefficients.
    [[nodiscard]] std::vector<double> getDenominatorCoefficients() const noexcept;
    /// @brief Sets the equality tolerance.
    /// @param[in] tol  The maximum absolute tolerance between coefficients
    ///                 when checking for inequality.  Recommend 1.e-12.
    void setEqualityTolerance(double tol);
private:
    class BAImpl;
    std::unique_ptr<BAImpl> pImpl;
}; // End BA
std::ostream& operator<<(std::ostream &os, const BA &ba);
}
#endif
