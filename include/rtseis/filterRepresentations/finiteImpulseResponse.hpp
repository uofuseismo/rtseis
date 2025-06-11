#ifndef RTSEIS_FILTER_REPRESENTATIONS_FINITE_IMPULSE_RESPONSE_HPP
#define RTSEIS_FILTER_REPRESENTATIONS_FINITE_IMPULSE_RESPONSE_HPP
#include <rtseis/vector.hpp>
namespace RTSeis::FilterRepresentations
{
/// @name FiniteImpulseResponse "finiteImpulseReseponse.hpp"
/// @brief A filter whose impulse response is of finite duration and is
///        represented by
/// \f[
///    y[n] = \sum_{i=0}^{N} b_i \cdot x[n - i].
/// \f].
/// where \f$ x[n] \f$ is the input signal, \f$ y[n] \f$ the output signal,
/// \f$ N \f$ the filter order, and \f$ b_i \f$ the i'th filter coefficient
/// of the FIR filter.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<class T>
class FiniteImpulseResponse
{
public:
    /// @brief Constructs an FIR filter from the input filter coefficients.
    /// @param[in] b   The FIR filter coefficients.
    /// @throws std::invalid_argument if the b is empty.
    explicit FiniteImpulseResponse(const RTSeis::Vector<T> &b);
    /// @brief Copy constructor.
    /// @param[in] firFilter  The FIR filter to copy to this.
    FiniteImpulseResponse(const FiniteImpulseResponse &firFilter); 
    /// @brief Move constructor.
    /// @param[in,out] firFilter  The FIR filter whose memory will be moved to 
    ///                           this.  On exit, firFilter's behavior is
    ///                           undefined.
    FiniteImpulseResponse(FiniteImpulseResponse &&firFilter) noexcept;
   
    /// @result The FIR filter coefficients.
    [[nodiscard]] RTSeis::Vector<T> getFilterCoefficients() const noexcept;
    /// @result A reference to the FIR filter coefficients.
    [[nodiscard]] const RTSeis::Vector<T> &getFilterCoefficientsReference() const noexcept;
    /// @result The order of the filter.
    [[nodiscard]] int getOrder() const noexcept; 

    /// @brief Copy assignment operator.
    /// @param[in] firFilter  The FIR filter to copy to this.
    FiniteImpulseResponse& operator=(const FiniteImpulseResponse &firFilter);
    /// @brief Move assignment operator.
    /// @param[in,out] firFilter  The FIR filter to copy to this.
    ///                           On exit, FIR filter's behavior is undefined. 
    FiniteImpulseResponse& operator=(FiniteImpulseResponse &&firFilter) noexcept;

    /// @brief Destructor.
    ~FiniteImpulseResponse();
    FiniteImpulseResponse() = delete;
private:
    class FiniteImpulseResponseImpl;
    std::unique_ptr<FiniteImpulseResponseImpl> pImpl;
};
}
#endif
