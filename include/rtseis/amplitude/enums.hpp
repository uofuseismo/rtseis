#ifndef RTSEIS_AMPLITUDE_ENUMS_HPP
#define RTSEIS_AMPLITUDE_ENUMS_HPP 1
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
    Sine    /*!< Apply a sine taper prior to filtering */
};
enum class HighPassFilter
{
    No, 
    Yes 
};
}
#endif
