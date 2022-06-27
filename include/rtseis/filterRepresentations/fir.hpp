#ifndef RTSEIS_FILTERREPRESENTATIONS_FIR_HPP
#define RTSEIS_FILTERREPRESENTATIONS_FIR_HPP 1
#include <vector>
#include <ostream>
#include <memory>
namespace RTSeis::FilterRepresentations
{
/// @class FIR fir.hpp "include/rtseis/filterRepresentations/fir.hpp"
/// @brief Represents a finite impulse response filter in terms of filter taps.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_utils_fr
class FIR
{
public:
    /// @name Constructors
    /// @{
    /// @brief Default constructor.
    FIR();
    /// @brief Constructs an FIR class from the given filter taps.
    /// @param[in] pTaps   Filter taps to set.
    /// @throws std::invalid_argument if pTaps is empty.
    explicit FIR(const std::vector<double> &pTaps);
    /// @brief Copy constructor.
    /// @param[in] fir  Class from which to initialize this class.
    FIR(const FIR &fir);
    /// @brief Move constructor.
    /// @param[in,out] fir  Class to whose memory will be moved to this class.
    ///                     On exit fir's behavior is undefined.
    FIR(FIR &&fir) noexcept;
    /// @}
 
    /// @name Operators
    /// @{
    /// @brief Assignment operator.
    /// @param[in] fir   Class to copy.
    /// @result A deep copy of the input class.
    FIR& operator=(const FIR &fir);
    /// @brief Move operator.
    /// @param[in,out] fir  FIR class to move.  On exit fir's behavior
    ///                     is undefined.
    /// @result The memory on FIR moved to this.
    FIR& operator=(FIR &&fir) noexcept;
    /// @brief Equality operator.
    /// @param[in] fir  Class to compare to this class.
    /// @result True indicates that fir equals this class to within a given
    ///         tolerance.
    [[nodiscard]] bool operator==(const FIR &fir) const;
    /// @brief Inequality operator.
    /// @param[in] fir   Class to compare to this class.
    /// @result True indicates that fir does not equal this class to within a
    ///         given tolerance. 
    [[nodiscard]] bool operator!=(const FIR &fir) const;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Default destructor.
    ~FIR();
    /// @brief Clears the structure.
    void clear() noexcept;
    /// @}

    /// @brief Prints the FIR structure.
    /// @param[in] fout   File handle to print to.  If fout is NULL then this
    ///                   will print to stdout.
    void print(FILE *fout = stdout) const noexcept;

    /// @result Number of filter taps.
    [[nodiscard]] int getNumberOfFilterTaps() const noexcept;

    /// @name Set Filter Taps
    /// @{
    /// @brief Sets the filter taps.
    /// @param[in] n     The number of filter taps.
    /// @param[in] taps  The filter coefficients to set.  This is an array
    ///                  of dimension [n].
    /// @throws std::invalid_argument if n is less than 1 or taps is NULL.
    void setFilterTaps(size_t n, const double taps[]);
    /// @brief Sets the filter taps.
    /// @param[in] taps  The filter coefficients to set.
    /// @throws std::invalid_argument if taps is empty.
    void setFilterTaps(const std::vector<double> &taps);
    /// @}

    /// @name Get Filter Taps
    /// @{
    /// @result The filter coefficients.
    [[nodiscard]] std::vector<double> getFilterTaps() const noexcept;
    /// @}

    /// @brief Sets the tolerance when testing for equality.
    /// @param[in] tol  The maximum absolute tolerance between coefficients
    ///                 when checking for inequality.  Recommend 1.e-12.
    void setEqualityTolerance(double tol);
  
private:
    class FIRImpl;
    std::unique_ptr<FIRImpl> pImpl;
};
std::ostream& operator<<(std::ostream &os, const FIR &fir);
}
#endif
