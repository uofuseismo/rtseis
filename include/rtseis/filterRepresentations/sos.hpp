#ifndef RTSEIS_FILTERREPRESENTATIONS_SOS_HPP
#define RTSEIS_FILTERREPRESENTATIONS_SOS_HPP 1
#include <vector>
#include <ostream>
#include <memory>
namespace RTSeis::FilterRepresentations
{
/// @class SOS sos.hpp "include/rtseis/filterRepresentations/sos.hpp"
/// @brief Represents an infinite impulse response filter in terms of a
///        series of cascaded second-order-sections.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_utils_fr
class SOS
{
public:
    /// @name Constructor
    /// @{

    /// @brief Default constructor. 
    SOS();
    /// @brief Copy constructor.
    /// @param[in] sos  SOS class from which to initialize.
    SOS(const SOS &sos);
    /// @brief Move constructor.
    /// @param[in,out] sos  SOS class to move to this class.  
    ///                     On exit sos is undefined.
    SOS(SOS &&sos) noexcept;
    /// @brief Constructs an SOS class from the given second order sections.
    /// @param[in] ns    Number of sections.  This must be positive.
    /// @param[in] bs    The numerator coefficients to set.  This 
    ///                  is an array of dimension [3 x ns] with leading
    ///                  dimension 3.  Moreover, the leading coefficient
    ///                  of each second-order-section cannot be 0.
    /// @param[in] as    The denoiminator coefficients to set.  This is
    ///                  is an array of dimension [3 x ns] with leading
    ///                  dimension 3.  Moreover, the leading coefficient
    ///                  of each second-order-section cannot be 0.
    /// @throws std::invalid_argument if the ns is less than 1 or the leading 
    ///         coefficient of as or bs is 0.
    SOS(int ns,
        const std::vector<double> &bs,
        const std::vector<double> &as);
    /// @}

    /// @name Operators
    /// @{

    /// @brief Assignement operator.
    /// @param[in] sos  SOS class to copy.
    /// @result A deep copy of the input class.
    SOS& operator=(const SOS &sos);
    /// @brief Move assignement operator.
    /// @param[in,out] sos  SOS class to move.  On exit sos's behavior is
    ///                     undefined.
    /// @result The moved version of sos.
    SOS& operator=(SOS &&sos) noexcept;
    /// @brief Equality operator.
    /// @param sos  Class to compare to this class.
    /// @result True indicates that sos equals this class to within a given
    ///         tolerance.
    [[nodiscard]] bool operator==(const SOS &sos) const;
    /// @brief Inequality operator.
    /// @param[in] sos  Class to compare to this class.
    /// @result True indicates that sos does not equal this class to within a
    ///         given tolerance.
    [[nodiscard]] bool operator!=(const SOS &sos) const;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Default destructor.
    ~SOS();
    /// @brief Clears the class.
    void clear() noexcept;
    /// @}

    /// @brief Prints the SOS structure.
    /// @param[in] fout   File handle to print to.  If fout is NULL then this
    ///                   will print to stdout.
    void print(FILE *fout = stdout) const noexcept;

    /// @brief Sets the second order sections on the class.
    /// @param[in] ns   The number of second order sections.
    /// @param[in] bs   The numerator coefficients.  This has dimension
    ///                 [3 x ns] with leading order 3.  Futhermore, the
    ///                 \f$ b[3 i_s] \f$ cannot be 0 for
    ///                 \f$ i=1,2, \cdots n_s \f$.
    /// @param[in] as   The denominator coefficients.  This has dimension
    ///                 [3 x ns] with leading order 3.  Futhermore,
    ///                 the \f$ b[3 i_s] \f$ cannot be 0 for
    ///                 \f$ i=1,2, \cdots n_s \f$.
    /// @throws std::invalid_argument if the ns is less than 1 or the leading 
    ///         coefficient of as or bs is 0.
    void setSecondOrderSections(int ns,
                                const std::vector<double> &bs,
                                const std::vector<double> &as);
    /// @result A vector holding the numerator coefficients.  This has
    ///        dimension [3 x ns] with leading dimension 3 where ns is
    ///        given by getNumberOfSections().
    [[nodiscard]] std::vector<double> getNumeratorCoefficients() const noexcept;
    /// @result A vector holding the denominator coefficients.  This has
    ///         dimension [3 x ns] with leading dimension 3 where ns is
    ///         given by getNumberOfSections().
    [[nodiscard]] std::vector<double> getDenominatorCoefficients() const noexcept;
    /// @result The number of second order sections.
    [[nodiscard]] int getNumberOfSections() const noexcept;
    /// @brief Sets the tolerance when testing for equality.
    /// @param[in] tol   Tolerance.  Recommend 1.e-12.
    void setEqualityTolerance(double tol);
private:
    class SOSImpl;
    std::unique_ptr<SOSImpl> pImpl;
}; // End SOS
/// @result The text version of the sos filter.
std::ostream& operator<<(std::ostream &os, const SOS &sos);
}
#endif
