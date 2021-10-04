#ifndef RTSEIS_AMPLITUDE_TIMEDOMAINWOODANDERSONPARAMETERS_HPP
#define RTSEIS_AMPLITUDE_TIMEDOMAINWOODANDERSONPARAMETERS_HPP 1
#include <memory>
namespace RTSeis::Amplitude
{
/// @brief Defines the units to which the recording instrument is 
///        sensitive.
enum class InputUnits
{
    Velocity,     /*!< The recording instrument is sensitive to ground
                       velocity.  This is most likely an "HH?" or "BH?"
                       channel code. */
    Acceleration  /*!< The recording instrument is sensitive to ground
                       acceleration.  This is most likely an "EN?"
                       channel code.  */
};
/// @brief Defines the gain constant for the Wood-Anderson filter.
enum class WoodAndersonGain
{
    WA_2800, /*!< This uses the classic gain of 2800. */
    WA_2080, /*!< This uses the `corrected' gain of 2080.  In this
                  case, the resulting signal will be divided by 2800
                  and multiplied by 2080. */
};
/// @brief Defines the detrend strategy for post-processing mode.
/// @note This is the first operation applied to the signal.
enum class DetrendType
{
    None,       /*!< Do not remove the mean or trend from the signal
                     prior to tapering. */
    RemoveMean, /*!< Remove the mean from the signal */
    RemoveTrend /*!< Remove a best-fitting line from the signal. */
};
/// @brief Defines the tapering window for post-processing.
/// @note This is the operation applied after detrending and before
///       filtering.
enum class WindowType
{
    None,   /*!< Do not apply a taper to the signal prior to filtering. */
    Cosine  /*!< Apply a cosine taper prior to filtering */
};

/// @class TimeDomainWoodAndersonParameters timeDomainWoodAndersonParameters.hpp "rtseis/amplitude/timeDomainWoodAndersonParameters.hpp"
/// @brief The parameters that control the time-domain Wood-Anderson
///        filtering described by Kanamori et al.'s 1999 paper:
///        Continuous monitoring of ground-motion parameters.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_amplitude
class TimeDomainWoodAndersonParameters
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    TimeDomainWoodAndersonParameters();
    /// @brief Copy constructor.
    /// @param[in] parameters  The Wood Anderson parameters class from which
    ///                        to initialize this class.
    TimeDomainWoodAndersonParameters(
        const TimeDomainWoodAndersonParameters &parameters);
    /// @brief Move constructor.
    /// @param[in,out] parameters  The Wood Anderson parameters class from
    ///                            which to initialize this class.  On exit,
    ///                            parameters's behavior is undefined.
    TimeDomainWoodAndersonParameters(
        TimeDomainWoodAndersonParameters &&parameters) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] parameters  The parameters class to copy to this.
    /// @result A deep copy of the input parameters.
    TimeDomainWoodAndersonParameters& operator=
        (const TimeDomainWoodAndersonParameters &parameters);
    /// @brief Move assignment operator.
    /// @param[in,out] parameters  The parameters class whose memory will be
    ///                            moved to this.  On exit, parameters's
    ///                            behavior is undefined.
    /// @result The memory from parameters moved to this.
    TimeDomainWoodAndersonParameters& operator=
        (TimeDomainWoodAndersonParameters &&parameters) noexcept;
    /// @}

    /// @name Required Parameters:
    /// @{
    /// @brief Allows the Wood-Anderson filter distinguish between acceleration
    ///        or velocity filtering.
    /// @param[in] units  The input units.
    void setInputUnits(InputUnits units) noexcept;
    /// @result The units to which the input signal is proportional.
    /// @throws std::runtime_error if \c haveInputUnits() is false.
    InputUntis getInputUnits() const noexcept;
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

    /// @brief Sets the sampling rate for this instrument.
    /// @param[in] samplingRate  The sampling rate in Hz.  This filter uses
    ///                          pre-optimzed constants so only a narrow
    ///                          set of sampling rates are supported.
    /// @throws std::invalid_argument if \c isSamplingRateSupported() is false.
    void setSamplingRate(double samplingRate);
    /// @result The sampling rate in Hz.
    /// @throws std::runtime_error if \c haveSamplingRate() is false.
    [[nodiscard]] double getSamplingRate() const;
    /// @result True indicates the sampling rate is set.
    [[nodiscard]] bool haveSamplingRate() const noexcept;
    /// @result True indicates the sampling rate is supported.
    [[nodiscard]] bool isSamplingRateSupported() const noexcept;

    /// @brief Sets the Wood-Anderson gain convention used by the filter.
    /// @param[in] gain  The Wood-Anderson gain convention.
    /// @note The default is the historic WA_2800.
    void setWoodAndersonGain(WoodAndersonGain gain) noexcept;
    /// @result The Wood-Anderson gain convention.
    [[nodiscard]] WoodAndersonGain getWoodAndersonGain() const noexcept;
    /// @}

    /// @name Additional for Post-Processing Options
    /// @{
    /// @brief Sets the taper to apply to the signal after detrending and
    ///        prior to filtering.
    /// @param[in] percentage  The percentage of the signal to which the taper
    ///                        will be applied.  For example, if 5 percent then
    ///                        the initial 2.5 percent and final 2.5 percent of 
    ///                        the signal will be tapered with the given window.
    ///                        This must be in the range [0,100].
    /// @param[in] window      The window to apply.  If this is WindowType::None
    ///                        then no filtering will be performed.
    void setTaper(double percentage = 5, WindowType window = WindowType::Cosine);
    /// @result The percentage to taper the signal prior to filtering.
    [[nodiscard]] double getTaperPercentage() const noexcept;
    /// @result The window type.
    [[nodiscard]] Taper getWindowType() const noexcept;
    
    /// @brief Sets the detrend strategy to apply prior to tapering.
    /// @param[in] detrendType   The detrend strategy.
    void setDetrendType(DetrendType detrendType);
    /// @result The detrend 
    [[nodiscard]] DetrendType getDetrendType() const noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Resets the class and releases memory
    void clear() noexcept;
    /// @brief Destructor.
    ~TimeDomainWoodAndersonParameters();
    /// @}
private:
    class TimeDomainWoodAndersonParametersImpl;
    std::unique_ptr<TimeDomainWoodAndersonParametersImpl> pImpl;
};
}
#endif
