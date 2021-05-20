#ifndef RTSEIS_UTILITIES_TRANSFORMS_SLIDINGWINDOWREALDFT_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_SLIDINWINDOWGREALDFT_HPP 1
#include <memory>
#include <vector>
#include <complex>
#include "rtseis/enums.hpp"
namespace RTSeis::Utilities::Transforms
{
// Forward declaration
class SlidingWindowRealDFTParameters;
/// @brief This is a utility class that slides a window over a real time series
///        and computes the DFT in each window.  This is the basis of many 
///        higher order spectral analysis methods such as the STFT, spectrogram,
///        and Welch's method.  
/// @note This class is not intended for use by a large audience.  Because it is
///       foundational to higher-order methods whose results are more desirable
///       this interface emphasizes efficiency at the sake of clarity.
/// @author Ben Baker, University of Utah
/// @copyright Ben Baker distributed under the MIT license.
/// @date May 2019.
template<class T = double>
class SlidingWindowRealDFT
{
public:
    /// @name Constructors
    /// @{
    /// @brief Default constructor.
    SlidingWindowRealDFT();
    /// @brief Copy constructor.
    /// @param[in] swdft  The sliding window real DFT class to initialize from.
    /// @throws std::runtime_error if the class fails to initialize. 
    ///         This is indicative of an internal error. 
    SlidingWindowRealDFT(const SlidingWindowRealDFT &swdft);
    /// @brief Move constructor.
    /// @param[in,out] swdft  The sliding window real DFT class to initialize
    ///                       from. On exit, swdft's behavior will be undefined.
    SlidingWindowRealDFT(SlidingWindowRealDFT &&swdft) noexcept;
    /// @}
    
    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] swdft  The class to copy.
    /// @result A deep copy of the sliding window real DFT class.
    /// @throws std::runtime_error if the class fails to initialize.
    ///         This is indicative of an internal error.
    SlidingWindowRealDFT& operator=(const SlidingWindowRealDFT &swdft);
    /// @brief Move assignment operator.
    /// @param[in,out] swdft  The class to move.  On exit swdft's behavior is
    ///                       undefined.
    /// @result The memory that was moved from swdft to this.
    SlidingWindowRealDFT& operator=(SlidingWindowRealDFT &&swdft) noexcept;
    /// @result The frequencies (Hz) at which the transform is computed.
    /// @param[in] samplingRate  The sampling rate (Hz) of the signal.
    /// @throws std::runtime_error if \c isInitialized() is false.
    /// @throws std::invalid_argument if samplingRate is not positive.
    [[nodiscard]] std::vector<T> getFrequencies(double samplingRate) const;
    /// @brief Gets the frequencies at which the power spectral density was
    ///        estimated. 
    /// @param[in] samplingRate  The sampling rate (Hz) of the signal.
    /// @param[in] nFrequencies  The number of frequencies.  This must match the
    ///                          result of \c getNumberOfFrequencies().
    /// @param[out] frequencies  The frequencies (Hz) at which the spectral
    ///                          density was estimated.  This is an array of
    ///                          dimension [nFrequencies].
    /// @throws std::invalid_argument if the sampling rate is not positive,
    ///         nFrequencies is invalid or frequencies is NULL.
    /// @throws std::runtime_error if \c isInitialized() is false.
    /// @sa \c getNumberOfFrequencies()
    void getFrequencies(double samplingRate, int nFrequencies, T *frequencies[]) const;
    /// @param[in] samplingRate  The sampling rate in Hz.
    /// @result The time windows (seconds) of each transform.  For example,
    ///         the it'th transform was computed for a time window given by
    ///         [times[it], times[it+1]).
    /// @throws std::runtime_error if \c isInitialized() is false.
    /// @throws std::invalid_argument if samplingRate is not positive.
    [[nodiscard]] std::vector<T> getTimeWindows(double samplingRate) const;
    /// @brief Gets the time windows in which each transform was computed.
    /// @param[in] samplingRate  The sampling rate (Hz) of the signal.
    /// @param[in] nTimes        The number of times.  This should be equal
    ///                          to \c getNumberOfTransformWindows() + 1.
    /// @param[out] times        The times (seconds) defining each bin.
    ///                          This is an array of dimension [nTimes].
    /// @throws std::invalid_argument if the sampling rate is not positive,
    ///         nTimes is incorrect, or times is NULL.
    /// @throws std::runtime_error if \c isInitialized() is false.
    void getTimeWindows(double samplingRate, int nTimes, T *times[]) const;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Default destructor.
    ~SlidingWindowRealDFT();
    /// @brief Releases memory on the class.
    void clear() noexcept;
    /// @}

    /// @brief Initializes the sliding window real discrete Fourier transform.
    /// @param[in] parameters   The sliding window real DFT parameters.
    ///                         This must be valid.
    /// @throws std::invalid_argument if parameters.isValid() is false.
    void initialize(const SlidingWindowRealDFTParameters &parameters);
    /// @result True indicates that the class is inititalized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @result The number of frequencies.
    /// @throws std::runtime_error if the class is not intitialized.
    [[nodiscard]] int getNumberOfFrequencies() const;
    /// @result The number of sliding time windows for which a
    ///         DFT was computed.  This is the number of columns in the
    ///        output matrix.
    /// @throws std::runtime_error if the class is not intitialized.
    [[nodiscard]] int getNumberOfTransformWindows() const;
    /// @brief Returns the expected number of samples in the time series.
    /// @throws std::runtime_error if the class is not initialized.
    [[nodiscard]] int getNumberOfSamples() const;

    /// @brief Computes the sliding window DFT of the real signal.
    /// @param[in] nSamples  The number of samples in the signal.
    ///                      This must match the result of 
    ///                      \c getNumberOfSamples().
    /// @param[in] x         The signal to transform.  This is an array whose
    ///                      dimension is [nSamples].
    /// @throws std::invalid_argument if any arguments are invalid.
    /// @throws std::runtime_error if the class is not initalized.
    /// @sa \c getNumberOfTransformWindow(), \c getNumberOfFrequencies()
    void transform(int nSamples, const T x[]);

    /// @brief Gets a pointer to the transform in the iWindow'th window.
    /// @param[in] iWindow  The window of the given transform.  This must
    ///                     be in the range 
    ///                     [0, \c getNumberOfTransformWindows()-1]. 
    /// @throws std::invalid_argument if iWindow is out of bounds.
    /// @throws std::runtime_error if \c transform() has not yet been called. 
    /// @sa \c transform(), \c getNumberOfTransformWindows()
    [[nodiscard]]
    const std::complex<T> *getTransform(int iWindow) const;
    /*! 
     * @brief Gets a pointer to the transform in the iWindow'th window.
     * @param[in] iWindow  The window of the given transform.  This must
     *                     be in the range
     *                     [0, \c getNumberOfTransformWindows()-1]. 
     * @throws std::invalid_argument if iWindow is out of bounds.
     * @throws std::runtime_error if the class precision is DOUBLE
     *         or the \c transform() has not yet been called. 
     * @sa \c transform(), \c getNumberOfTransformWindows(), \c getPrecision()
     */
    //[[nodiscard]]
    //const std::complex<float> *getTransform32f(const int iWindow) const;
private:
    class SlidingWindowRealDFTImpl;
    std::unique_ptr<SlidingWindowRealDFTImpl> pImpl;
};
}
#endif
