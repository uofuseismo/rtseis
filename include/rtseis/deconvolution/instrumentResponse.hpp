#ifndef RTSEIS_DECONVOLUTION_INSTRUMENTRESPONSE_HPP
#define RTSEIS_DECONVOLUTION_INSTRUMENTRESPONSE_HPP 1
#include <memory>
#include <vector>
// Forward declarations
namespace RTSeis::FilterRepresentations
{
class BA;
class ZPK;
}
/// @brief
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
    virtual ~InstrumentResponse();
    /// @brief Clears all memory and resets the class.
    void clear() noexcept;
    /// @}

    /// @name Sampling Rate
    /// @{
    /// @brief Sets the sampling rate.
    /// @param[in] df   The sampling rate in Hz.
    /// @throws std::invalid_argument if df is not positive.
    void setSamplingRate(double df);
    /// @result The sampling rate in Hz.
    /// @throw std::runtime_error if \c haveSamplingRate() is false.
    [[nodiscard]] double getSamplingRate() const;
    /// @result True indicates that the sampling rate was set.
    [[nodiscard]] bool haveSamplingRate() const noexcept;
    /// @}

    /// @name Analog Transfer Function
    /// @{
    /// @brief Sets the transfer function.
    /// @param[in] zpk  The analog poles, zeros, and gain defining the 
    ///                 transfer function.
    void setAnalogTransferFunction(
        const RTSeis::FilterRepresentations::ZPK &zpk) noexcept;
    /// @brief Sets the instrument's analog transfer function.
    /// @param[in] ba   The numerator and denominator coefficients defining the
    ///                 transfer function.
    /// @throws std::invalid_argument if there are no numerator or denominator
    ///         coefficients or the first denominator coefficient is 0.
    void setAnalogTransferFunction(
        const RTSeis::FilterRepresentations::BA &ba);
    /// @}

    /// @name Digital Transfer Function
    /// @{
    /// @brief Sets the instrument's digital transfer function.
    /// @param[in] zpk  The digital poles, zeros, and gain defining the response.
    void setDigitalTransferFunction(
        const RTSeis::FilterRepresentations::ZPK &zpk) noexcept;
    /// @param[in] ba   The digital transfer function.
    /// @throws std::invalid_argument if there are no numerator or denominator
    ///         coefficients or all the denominator coefficients are zero.
    void setDigitalTransferFunction(
        const RTSeis::FilterRepresentations::BA &ba);
    /// @}

    /// @result True indicates that the response was set.
    /// @sa \c setAnalogTransferFunction(), \c setDigitalTransferFunction()
    [[nodiscard]] bool haveTransferFunction() const noexcept;
    /// @result True indicates that this is an analog response.
    /// @throws std::runtime_error if the response is not yet set.
    /// @sa \c haveTransferFunction()
    [[nodiscard]] bool isAnalogTransferFunction() const;

    /// @brief Computes the instrument response at the given frequencies.
    /// @param[in] frequencies  The frequencies in Hz at which to compute the
    ///                         response.
    /// @result The instrument response tabulated at each frequency.
    virtual std::vector<std::complex<double>> 
        compute(const std::vector<double> &frequencies) const;
    /// @brief Computes the instrument response at the given frequencies.
    /// @param[in] nFrequencies  The number of frequencies at which to compute
    ///                          the response.
    /// @param[in] frequencies   The frequencies in Hz at which to compute the
    ///                          response.  This is an array of dimension
    ///                          [nFrequencies].
    /// @param[in] response      The instrument response tabulated at each
    ///                          frequency.  This is an array of dimension
    ///                          [nFrquencies].
    /// @throws std::invalid_argument if nFrequencies is positive and
    ///         frequencies is NULL or response is NULL.
    /// @throws std::runtime_error if the transfer function or sampling rate
    ///         was not set.
    /// @sa \c haveTransferFunction() or \c haveSamplingRate() is false.
    virtual void compute(int nFrequencies, const double frequencies[],
                         std::complex<double> *response[]) const;
private:
    class InstrumentResponseImpl;
    std::unique_ptr<InstrumentResponseImpl> pImpl;
};
}
#endif
