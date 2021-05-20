#ifndef RTSEIS_UTILITIES_TRANSFORMS_SPECTROGRAM_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_SPECTROGRAM_HPP 1
#include <memory>
#include <vector>
#include "rtseis/utilities/transforms/enums.hpp"

namespace RTSeis::Utilities::Transforms
{
class SlidingWindowRealDFTParameters;
/// @class Spectrogram "spectrogram.hpp" "rtseis/utilities/transforms/spectrogram.hpp"
/// @brief Computes the spectrogram of a time series.  This is accomplished
///        by splitting the time series into (potentially overlapping) segments
///        and Fourier transforming each segment.
/// @author Ben Baker
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
/// @date July 2019
template<class T = double>
class Spectrogram 
{
public:
    /// @name Constructors
    /// @{
    /// @brief Default constructor.
    Spectrogram();
    /// @brief Copy constructor.
    /// @param[in] spectrogram  The spectrogram class from which to initialize
    ///                         this class.
    Spectrogram(const Spectrogram &spectrogram);
    /// @brief Move constructor.
    /// @param[in,out] spectrogram  The spectrogram class from which to
    ///                             initialize this class.  On exit,
    ///                             spectrogram's behavior is undefined.
    Spectrogram(Spectrogram &&wspectrogram) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] spectrogram  The class to copy to this.
    /// @result A deep copy of the spectrogram class.
    /// @throws std::runtime_error if there is an initialization error.
    ///         This is indicative of an internal error.
    Spectrogram& operator=(const Spectrogram &spectrogram);
    /// @brief Move assignment operator.
    /// @param[in,out] spectrogram  The class to move to this.
    ///                             On exit spectrogram's behavior is undefined.
    /// @result The memory that was moved from spectrogram to this.
    Spectrogram& operator=(Spectrogram &&spectrogram) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Default destructor.
    ~Spectrogram();
    /// @brief Releases memory on the class.
    void clear() noexcept;
    /// @}

    /// @name Step 1: Initialization
    /// @{
    /// @brief Initializes the spectrogram transform.
    /// @param[in] parameters    The sliding window DFT parameters.
    /// @param[in] samplingRate  The sampling rate in Hz.
    /// @throws std::invalid_argument if parameters.isValid() is false
    ///         or the sampling rate is not positive.
    void initialize(const SlidingWindowRealDFTParameters &parameters,
                    double samplingRate = 1.0);
    /// @result True indicates that the class is inititalized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @brief Returns the expected number of time domain samples in the signal
    ///        to transform.
    /// @result The number of samples the signal to transform to contain should
    ///         contain.
    /// @throws std::runtime_error if the class is not inititalized.
    [[nodiscard]] int getNumberOfSamples() const;
    /// @result The number of frequencies.
    /// @throws std::runtime_error if the class is not intitialized.
    [[nodiscard]] int getNumberOfFrequencies() const;
    /// @result The number of windowsin which a transform was computed.
    /// @throws std::runtime_error if the class is not intitialized.
    [[nodiscard]] int getNumberOfTransformWindows() const;
    /// @result The frequencies (Hz) at which the power spectral density
    ///         was estimated.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] std::vector<T> getFrequencies() const;
    /// @brief Gets the frequencies at which the power spectral density was
    ///        estimated. 
    /// @param[in] nFrequencies  The number of frequencies.  This must match the
    ///                          result of \c getNumberOfFrequencies().
    /// @param[out] frequencies  The frequencies (Hz) at which the spectral
    ///                          density was estimated.  This is an array of
    ///                          dimension [nFrequencies].
    /// @throws std::invalid_argument if nFrequencies is invalid or frequencies
    ///         is NULL.
    /// @throws std::runtime_error if \c isInitialized() is false.
    /// @sa \c getNumberOfFrequencies()
    void getFrequencies(int nFrequencies, T *frequencies[]) const;
    /// @result The time windows (seconds) of each transform.  For example,
    ///         the it'th transform was computed for a time window given by
    ///         [times[it], times[it+1]).
    /// @throws std::runtime_error if \c isInitialized() is false.
    /// @throws std::invalid_argument if samplingRate is not positive.
    [[nodiscard]] std::vector<T> getTimeWindows() const;
    /// @brief Gets the time windows in which each transform was computed.
    /// @param[in] nTimes        The number of times.  This should be equal
    ///                          to \c getNumberTransformWindows() + 1.
    /// @param[out] times        The times (seconds) defining each bin.
    ///                          This is an array of dimension [nTimes].
    /// @throws std::invalid_argument if the sampling rate is not positive,
    ///         nTimes is incorrect, or times is NULL.
    /// @throws std::runtime_error if \c isInitialized() is false.
    void getTimeWindows(int nTimes, T *times[]) const;
    /// @}

    /// @name Step 2: Transform
    /// @{
    /// @brief Computes the spectrogram of a signal.
    /// @param[in] nSamples   The number of samples in the signal.  This must
    ///                       equal the result of \c getNumberOfSamples().
    /// @param[in] x          The signal to transform.
    /// @throws std::invalid_argument if nSamples or x is invalid.
    /// @throws std::runtime_error if the class is not initialized.
    /// @sa \c isInitialized(), \c getNumberOfSamples()
    void transform(int nSamples, const T x[]);
    /// @retval True indicates that the transform has been computed.
    [[nodiscard]] bool haveTransform() const noexcept;
    /// @}

    /// @name Step 3: Get Results
    /// @{
    /// @result The amplitude spectrogram.  This is an 
    ///         [\c getNumberOfTransformWindows x \c getNumberOfFrequencie) ]
    ///         matrix with leading dimension \c getNumberOfFrequencies().
    /// @throws std::runtime_error if \c haveTransform() is false.
    std::vector<T> getAmplitude() const;
    /// @result A pointer to the amplitude spectrum.  This is an
    ///         [\c getNumberOfTransformWindows x \c getNumberOfFrequencies ]
    ///         matrix with leading dimension \c getNumberOfFrequencies().
    /// @throws std::runtime_error if \c haveTransform() is false.
    const T *getAmplitudePointer() const;
    /// @result The phase spectrogram in radians.  This is an 
    ///         [\c getNumberOfTransformWindows x \c getNumberOfFrequencie) ]
    ///         matrix with leading dimension \c getNumberOfFrequencies().
    /// @throws std::runtime_error if \c haveTransform() is false.
    std::vector<T> getPhase() const;
    /// @result A pointer to the phase spectrum in radians.  This is an
    ///         [\c getNumberOfTransformWindows x \c getNumberOfFrequencies ]
    ///         matrix with leading dimension \c getNumberOfFrequencies().
    /// @throws std::runtime_error if \c haveTransform() is false.
    const T *getPhasePointer() const;
    /// @}
private:
    class SpectrogramImpl;
    std::unique_ptr<SpectrogramImpl> pImpl;
};
}
#endif
