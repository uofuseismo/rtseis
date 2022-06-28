#ifndef RTSEIS_FILTERFILTERIMPLEMENTATION_IIRFILTER_HPP
#define RTSEIS_FILTERFILTERIMPLEMENTATION_IIRFILTER_HPP 1
#include <memory>
#include "rtseis/enums.hpp"
#include "rtseis/filterImplementations/enums.hpp"
namespace RTSeis::FilterRepresentations
{
class BA;
}
namespace RTSeis::FilterImplementations
{
/// @class IIRFilter iirFilter.hpp "rtseis/filterImplementations/iirFilter.hpp"
/// @brief This is the core implementation for IIR filtering.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_filterImplemenations
template<RTSeis::ProcessingMode E = RTSeis::ProcessingMode::POST,
         class T = double>
class IIRFilter
{
public:
    /// @name Constructors
    /// @{
    /// @brief Default destructor.
    IIRFilter();
    /// @brief Copy constructor.
    /// @param[in] iir  IIR filter class from which to initialize this class.
    IIRFilter(const IIRFilter &iir);
    /// @brief Move constructor.
    /// @param[in,out] iir   The IIR filter class from which to initialize this
    ///                      class.  On exit, iir's behavior is undefined.
    IIRFilter(IIRFilter &&iir) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignent operator.
    /// @param[in] iir  IIR filter class to copy.
    /// @result A deep copy of the IIR filter class.
    IIRFilter &operator=(const IIRFilter &iir);
    /// @brief Move assignment operator.
    /// @param[in,out] iir  The IIR filter class to move to this.  On exit,
    ///                     IIR's behavior is undefined.
    /// @result The memory from iir moved to this.
    IIRFilter &operator=(IIRFilter &&iir) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Default destructor.
    ~IIRFilter();
    /// @}

    /// @brief Initializes the IIR filter.
    /// @param[in] ba  The IIR filter.
    /// @param[in] implementation  The filter implementation.
    /// @throws std::invalid_argument if the IIR filter is invalid.
    void initialize(const RTSeis::FilterRepresentations::BA &ba,
                    IIRDFImplementation implementation = IIRDFImplementation::DF2_FAST);
    /// @brief Initializes the IIR filter.
    /// @param[in] nb     Number of numerator coefficients.
    /// @param[in] b      The numerator coefficients.  This has
    ///                   dimension [nb].
    /// @param[in] na     Number of denominator coefficients.
    /// @param[in] a      The denominator coefficients.  This has
    ///                   dimension [na].
    /// @param[in] implementation  Defines the algorithmic filter
    ///                            implementation.
    /// @throws std::invalid_argument if any arguments are invalid.
    void initialize(int nb, const double b[],
                    int na, const double a[],
                    IIRDFImplementation implementation = IIRDFImplementation::DF2_FAST);
    /// @result True indicates that the module is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @result The length of the initial conditions array.
    [[nodiscard]] int getInitialConditionLength() const;
    /// @brief Sets the initial conditions for the filter.  This should
    ///        be called prior to filter application as it will reset
    ///        the filter.
    /// @param[in] nz   The IIR filter initial condition length.
    ///                 This should be equal to
    ///                 getInitialConditionLength().
    /// @param[in] zi   The initial conditions.  This has dimension [nz].
    /// @throws std::invalid_argument if nz is invalid or nz is positive
    ///         zi is NULL.
    /// @throws std::runtime_error if the class is not initialized.
    void setInitialConditions(int nz, const double zi[]);
    /// @brief Applies the IIR filter.  Note, the class must be
    ///        initialized prior to using this function.
    /// @param[in] n   The number of points in the signal.
    /// @param[in] x   The input signal to filter.  This has dimension [n].
    /// @param[out] y  The filtered signal.  This has dimension [n].
    /// @throws std;:invalid_argument if n is positive and x or y is NULL.
    /// @throws std::runtime_error if the class is not initialized.
    void apply(int n, const T x[], T *y[]);
    /// @brief Resets the initial conditions to those set in
    ///        \c setInitialConditions().
    /// @throws std::runtime_error if the class is not initialized.
    void resetInitialConditions();
    /// @brief Clears memory and resets the filter.  After applying
    ///        this function the filter must be re-initialized prior
    ///        to being applied to the data.
    void clear() noexcept;
private:
    class IIRFilterImpl;
    std::unique_ptr<IIRFilterImpl> pImpl;
};
}
#endif
