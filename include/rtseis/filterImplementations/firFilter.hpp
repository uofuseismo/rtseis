#ifndef RTSEIS_FILTERIMPLEMENTATIONS_FIRFILTER_HPP
#define RTSEIS_FILTERIMPLEMENTATIONS_FIRFILTER_HPP 1
#include <memory>
#include "rtseis/enums.hpp"
#include "rtseis/filterImplementations/enums.hpp"
namespace RTSeis::FilterRepresentations
{
class FIR;
}
namespace RTSeis::FilterImplementations
{
/// @class FIRFilter firFilter.hpp "rtseis/filterImplementations/firFilter.hpp"
/// @brief This is the core implementation for FIR filtering.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_filterImplemenations
template<RTSeis::ProcessingMode E = RTSeis::ProcessingMode::POST,
         class T = double>
class FIRFilter
{
public:
    /// @name Constructors
    /// @{
    /// @brief Default constructor.
    FIRFilter();
    /// @brief Copy constructor.
    /// @param[in] fir   FIR class from which to initialize.
    FIRFilter(const FIRFilter &fir);
    /// @brief Move constructor.
    /// @param[in,out] fir   FIR class from which to initialize this class.
    ///                      On exit, fir's behavior is undefined.
    FIRFilter(FIRFilter &&fir) noexcept; 
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy operator.
    /// @param[in] fir   FIR class to copy.
    /// @result A deep copy of the FIR class.
    FIRFilter& operator=(const FIRFilter &fir);
    /// @brief Move assignment operator.
    /// @param[in,out] fir   FIR class whose memory will be moved to this.
    ///                      On exit, fir's behavior is undefined.
    /// @result FIR's memory moved to this.
    FIRFilter& operator=(FIRFilter &&fir) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Default destructor.
    ~FIRFilter();
    /// @}

    /// @brief Initializes the FIR filter.
    /// @param[in] fir   The FIR filter.
    /// @throws std::invalid_argument if the FIR filter is invalid.
    void initialize(const RTSeis::FilterRepresentations::FIR &fir,
                    FIRImplementation implementation = FIRImplementation::DIRECT);
    /// @brief Initializes the FIR filter.
    /// @param[in] nb    Number of numerator coefficients.
    /// @param[in] b     Numerator coefficients.  This is an array of
    ///                  dimension [nb].
    /// @param[in] mode  The processing mode.  By default this
    ///                  is for post-processing.
    /// @param[in] implementation  Defines the implementation.
    ///                            The default is to use the direct form.
    /// @throws std::invalid_argument if any of the arguments are invalid.
    void initialize(int nb, const double b[],
                    FIRImplementation implementation = FIRImplementation::DIRECT);
    /// @result True indicates that the module is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @result The length of the initial condition array.
    /// @throws std::runtime_error if the class is not initialized.
    [[nodiscard]] int getInitialConditionLength() const;
    /// @brief Sets the initial conditions for the filter.  This should
    ///        be called prior to filter application as it will reset
    ///        the filter.
    /// @param[in] nz   The FIR filter initial condition length.
    ///                 This should be equal to
    ///                 getInitialConditionLength().
    /// @param[in] zi   The initial conditions.  This has dimension [nz].
    /// @throws std::invalid_argument if nz is invalid or nz is positive
    ///         and zi is NULL.
    /// @throws std::runtime_error if the class is not initialized.
    void setInitialConditions(int nz, const double zi[]);
    /// @brief Gets a copy of the initial conditions.
    /// @param[in] nz   The FIR filter initial condition length.
    ///                 This should be equal to
    ///                 getInitialConditionLength().
    /// @param[out] zi  The initial conditions.  This has dimension [nz].
    /// @throws std::invalid_argument if nz is invalid or nz is positive
    ///        and zi is NULL.
    /// @throws std::runtime_error if the class is not initialized.
    void getInitialConditions(int nz, double *zi[]) const;
    /// @brief Applies the FIR filter to the data.
    /// @param[in] n   Number of points in signals.
    /// @param[in] x   Signal to filter.  This has dimension [n].
    /// @param[out] y  The filtered signal.  This has dimension [n].
    /// @throws std::invalid_argument if n is positive and x or y is NULL.
    /// @throws std::runtime_error if the class is not initialized.
    void apply(int n, const T x[], T *y[]);
    /// @brief Resets the initial conditions on the source delay line to
    ///        the default initial conditions or the initial conditions
    ///        set when FIRFilter::setInitialConditions() was called.
    /// @throws std::runtime_error if the class is not initialized.
    void resetInitialConditions();
    /// @brief Clears the module and resets all parameters.
    void clear() noexcept;
private:
    class FIRImpl;
    std::unique_ptr<FIRImpl> pImpl;
};
}
#endif
