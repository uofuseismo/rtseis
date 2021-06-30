#ifndef RTSEIS_DECONVOLUTION_INSTRUMENTRESPONSE_HPP
#define RTSEIS_DECONVOLUTION_INSTRUMENTRESPONSE_HPP 1
#include <memory>
// Forward declarations
namespace RTSeis::FilterRepresentations
{
class BA;
class ZPK;
}
namespace RTSeis::Deconvolution
{
/// @class InstrumentResponse instrumentResponse.hpp "rtseis/deconvolution/instrumentResponse.hpp"
/// @brief A container for storing the instrument transfer function, usually
///        expressed as poles, zeros, and gain, and for computing the
///        corresponding instrument response.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class InstrumentResponse
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    InstrumentResponse(); 
    /// @brief Copy constructor.
    /// @param[in] response  The instrument response class from which to
    ///                      initialize this class.
    InstrumentResponse(const InstrumentResponse &response);
    /// @brief Move constructor.
    /// @param[in,out] response  The instrument response class from which to
    ///                          intialize this class.
    ///                          On exit, response's behavior is undefined.
    InstrumentResponse(InstrumentResponse &&response);
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] response  The instrument response to copy.
    /// @result A deep copy of the response.
    InstrumentResponse& operator=(const InstrumentResponse &response);
    /// @brief Move assignment operator.
    /// @param[in,out] response  The response whose memory is moved to this.
    ///                          On exit response's behavior is undefined.
    /// @result Contains the moved memory from response.
    InstrumentResponse& operator=(InstrumentResponse &&response) noexcept;
    /// @}

    //// @name Destructors
    /// @{
    /// @brief Destructor.
    ~InstrumentResponse();
    /// @brief Clears all memory and resets the class.
    void clear() noexcept;
    /// @}

    /// @brief Sets the sampling rate.
    /// @param[in] df   The sampling rate in Hz.
    /// @throws std::invalid_argument if df is not positive.
    void setSamplingRate(double df);
    /// @result The sampling rate in Hz.
    /// @throw std::runtime_error if \c haveSamplingRate() is false.
    [[nodiscard]] double getSamplingRate() const;
    /// @result True indicates that the sampling rate was set.
    [[nodiscard]] bool haveSamplingRate() const noexcept;

    /// @name Analog Transfer Function
    /// @{
    /// @brief Sets the analog digital response.
    /// @param[in] zpk  The analog poles, zeros, and gain defining the response.
    void setAnalogResponse(const RTSeis::FilterRepresentations::ZPK &zpk) noexcept;
    /// @brief Sets the analog digital response.
    /// @param[in] ba   The numerator and denominator coefficients defining the
    ///                 response.
    /// @throws std::invalid_argument if there are no numerator or denominator
    ///         coefficients or the first denominator coefficient is 0.
    void setAnalogResponse(const RTSeis::FilterRepresentations::BA &ba);
    /// @}

    /// @name Digital Transfer Function
    /// @{
    /// @brief Sets a digital instrument response.
    /// @param[in] zpk  The digital poles, zeros, and gain defining the response.
    void setDigitalResponse(const RTSeis::FilterRepresentations::ZPK &zpk) noexcept;
    /// @brief Sets the digital instrument response.
    /// @param[in] ba   The numerator and denominator coefficients defining the
    ///                 response.
    /// @throws std::invalid_argument if there are no numerator or denominator
    ///         coefficients or the first denominator coefficient is 0.
    void setDigitalResponse(const RTSeis::FilterRepresentations::BA &ba);
    /// @}

    /// @brief Computes the instrument response at the given frequencies.
    /// @param[in] nfreqs       The number of frequencies at which to compute
    ///                         the response.
    /// @param[in] frequencies  The frequencies in Hz at which to compute the
    ///                         response.  This is an array of dimension
    ///                         [nfreqs].
    /// @param[in] response     The instrument response tabulated at each
    ///                         frequency.  This is an array of dimension
    ///                         [nfreqs].
    /// @throws std::invalid_argument if nfreqs is positive and frequencies is
    ///         NULL or response is NULL.
    /// @throws std::runtime_error if the response was not set.
    /// @sa \c haveResponse()
    void computeResponse(int nfreqs, const double frequencies[],
                         std::complex<double> response[]);
    /// @result True indicates that the response was set.
    /// @sa \c setAnalogResponse(), \c setDigitalResponse()
    [[nodiscard]] bool haveResponse() const noexcept;
    /// @result True indicates that this is an analog response.
    /// @throws std::runtime_error if the response is not yet set.
    /// @sa \c haveResponse()
    [[nodiscard]] bool isAnalogResponse() const;
     
private:
    class InstrumentResponseImpl;
    std::unique_ptr<InstrumentResponseImpl> pImpl;
};
}
#endif
