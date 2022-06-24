#ifndef RTSEIS_AMPLITUDE_TAUPPARAMETERS_HPP
#define RTSEIS_AMPLITUDE_TAUPPARAMETERS_HPP
#include <memory>
#include "rtseis/amplitude/enums.hpp"
namespace RTSeis::FilterRepresentations
{
class SOS;
}
namespace RTSeis::Amplitude
{
/// @class TauPParameters "tauPParameters.hpp" "rtseis/amplitude/tauPParameters.hpp"
/// @brief Sets the parameters for the Berkeley \f$ \tau_p \f$ method.
///        Note, a lot of pertinent details can be found in the supplement of
///        https://doi.org/10.1029/2008GL035611.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class TauPParameters
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    TauPParameters();
    /// @brief Copy constructor.
    /// @param[in] parameters  The parameters class from which to initialize
    ///                        this class.
    TauPParameters(const TauPParameters &parameters);
    /// @brief Move constructor.
    /// @param[in,out] parameters  The parameters class from which to initialize
    ///                            this class.  On exit, the behavior of
    ///                            parameters is undefined.
    TauPParameters(TauPParameters &&parameters) noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] parameters  The parameters class to copy to this.
    /// @result A deep copy of the parameters class.
    TauPParameters& operator=(const TauPParameters &parameters);
    /// @brief Move assignment.
    /// @param[in,out] parameters  The parameters class to copy to this.
    ///                            On exit, parameters's behavior is undefined.
    TauPParameters& operator=(TauPParameters &&parameters) noexcept;
    /// @}

    /// @name Required Parameters 
    /// @{

    /// @brief Allows the Wood-Anderson filter distinguish between acceleration
    ///        or velocity filtering.
    /// @param[in] units  The input units.
    void setInputUnits(InputUnits units) noexcept;
    /// @result The units to which the input signal is proportional.
    /// @throws std::runtime_error if \c haveInputUnits() is false.
    [[nodiscard]] InputUnits getInputUnits() const;
    /// @result True indicates that the input units were set.
    [[nodiscard]] bool haveInputUnits() const noexcept;

    /// @result After dividing by this number the input signal will be
    ///         proportional to meters/second or meters/second/second.
    /// @throws std::invalid_argument if the gain is 0.
    void setSimpleResponse(double gain);
    /// @result The instrument simple response to make it proportional to
    ///         meters/second or meters/second/second.
    [[nodiscard]] double getSimpleResponse() const;
    /// @result True indicates the simple response was set. 
    [[nodiscard]] bool haveSimpleResponse() const noexcept;

    /// @param[in] samplingRate  The signal sampling rate in Hz.
    /// @note This this will change the lowpass filter and smoothing parameter.
    void setSamplingRate(double samplingRate);
    /// @result The sampling rate in Hz.
    /// @throws std::runtime_error if the sampling rate is not set.
    [[nodiscard]] double getSamplingRate() const;
    /// @result True indiciates the sampling rate was set.
    [[nodiscard]] bool haveSamplingRate() const noexcept;
    /// @}

    /// @name Optional Parameters
    /// @{

    /// @brief Stable integration requires the removal of long-period signals. 
    /// @param[in] sos  The (likely highpass) filter to apply to the
    ///                 signal prior. 
    /// @note When the sampling rate is set this will default to a 
    ///       Butterworth 4-pole lowpass filter with corner at 0.075 Hz.
    ///       Initially, a 2-pole filter was suggested, however, 
    ///       Shieh, Wu, and Allen suggest a 4, 5, or 6 pole filter lowers
    ///       the variance in the max(tau_p) estimates.
    void setAccelerationFilter(const FilterRepresentations::SOS &sos);
    /// @result The (highpass) acceleration filter.
    [[nodiscard]] FilterRepresentations::SOS getAccelerationFilter() const noexcept;

    /// @brief Prior to computing tauP a causal 3 Hz low-pass Butterworth
    ///        filter with 2 poles is applied. 
    void setVelocityFilter(const FilterRepresentations::SOS &sos);
    /// @result The (lowpass) velocity filter representation.
    [[nodiscard]] FilterRepresentations::SOS getVelocityFilter() const noexcept;

    /// @brief Sets the smoothing constant for updating the squared
    ///        velocity and acceleration signals:
    ///        \[
    ///           V_i = \alpha V_{i - 1} + v_i^2
    ///        \]
    ///        \[
    ///           V_i = \alpha V_{i - 1} + a_i^2
    ///        \]
    /// @throws std::invalid_argument if $\f \alpha \f$ is not in the range
    ///         [0,1).
    /// @note After setting the sampling rate \f$ \alpha \f$ will be set to
    ///       \f$ \alpha = 1 - \Delta T \f$ where \f$ \Delta T \f$ is the
    ///       sampling period.
    void setSmoothingParameter(double alpha); 
    /// @result The smoothing parameter.
    [[nodiscard]] double getSmoothingParameter() const noexcept;

    /// @brief Accelerometers are filtered with a simple high-pass RC filter
    ///        \[
    ///           a_i - Q a_{i-1}= (x_i - x_{i-1}) (1 + Q)/2G
    ///        \]
    ///        where \f$ x_i \f$ is the input signal (proportional to
    ///        acceleration, \f$ G \f$ the gain for the acceleratometer,
    ///        and \f$ Q \f$ the filter constant.  The accelerometers are
    ///        then integrated to velocity using the trapezoid rule
    ///        \[
    ///           v_i - Q v_{i-1} = (x_i + x_{i-1}) (\delta T)/2 (1 + Q)/2
    ///        \]
    ///        and simultaneously high-pass filtered.
    void setFilterConstantQ(double constant);
    /// @result The filter constant.
    [[nodiscard]] double getFilterConstantQ() const noexcept;

    /// @brief Sets the detrend strategy to apply prior to tapering.
    /// @param[in] detrendType   The detrend strategy.
    void setDetrendType(DetrendType detrendType);
    /// @result The detrend 
    [[nodiscard]] DetrendType getDetrendType() const noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class.
    void clear() noexcept;
    /// @brief Destructor. 
    ~TauPParameters();
    /// @}
private:
    class TauPParametersImpl;
    std::unique_ptr<TauPParametersImpl> pImpl;
};
}
#endif
