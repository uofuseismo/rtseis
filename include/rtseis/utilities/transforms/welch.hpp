#ifndef RTSEIS_UTILITIES_TRANSFORMS_WELCH_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_WELCH_HPP 1
#include <memory>
#include <vector>
#include "rtseis/utilities/transforms/enums.hpp"

namespace RTSeis::Utilities::Transforms
{
class SlidingWindowRealDFTParameters;
/// @class Welch "welch.hpp" "rtseis/utilities/transforms/welch.hpp"
/// @brief Estimates the power spectral density using Welch's method.
///        This amounts to first dividing the data into overlapping segments.
///        Next, a modified periodogram is computed in each segment. 
///        Finally, the modified periodogram for all segments are averaged. 
/// @author Ben Baker
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
/// @date July 2019
template<class T = double>
class Welch
{
public:
    /// @name Constructors
    /// @{
    /// @brief Default constructor.
    Welch();
    /// @brief Copy constructor.
    /// @param[in] welch  The Welch class from which to initialize this class.
    Welch(const Welch &welch);
    /// @brief Move constructor.
    /// @param[in,out] welch  The Welch class from which to initialize this class.
    ///                       On exit, welch's behavior is undefined.
    Welch(Welch &&welch) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] welch  The class to copy.
    /// @result A deep copy of the welch class.
    /// @throws std::runtime_error if there is an initialization error.
    ///         This is indicative of an internal error.
    Welch& operator=(const Welch &welch);
    /// @brief Move assignment operator.
    /// @param[in,out] welch  The class to move to this.
    ///                       On exit welch's behavior is undefined.
    /// @result The memory that was moved from welch to this.
    Welch& operator=(Welch &&welch) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Default destructor.
    ~Welch();
    /// @brief Releases memory on the class.
    void clear() noexcept;
    /// @}

    /// @name Step 1: Initialization
    /// @{
    /// @brief Initializes the Welch transform.
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
    /// @}

    /// @name Step 2: Transform
    /// @{
    /// @brief Computes the Welch transform of a signal.
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
    /// @result The power spectral density estimate at each frequency.  If
    ///         the input signal has units of \f$ Volts \f$ then this
    ///         has units of \f$ \frac{Volts^2}{Hz} \f$.
    /// @throws std::runtime_error if \c haveTransform() is false.
    std::vector<T> getPowerSpectralDensity() const; 
    /// @brief Gets the power spectral density estimate.
    /// @param[in] nFrequencies  The number of frequencies. This must match the
    ///                          result of \c getNumberOfFrequencies().
    /// @param[out] psd          The power spectral density.  This is an array
    ///                          of dimension [nFrequencies].  If the input 
    ///                          signal has units of \f$ Volts \f$ then this
    ///                          has units of \f$ \frac{Volts^2}{Hz} \f$.
    /// @throws std::invalid_argument if nFrequencies is invalid or psd is NULL.
    /// @throws std::runtime_error if the class is not initialized or the 
    ///         transform has not been computed.
    /// @sa \c getNumberOfFrequencies()
    /// @sa \c isInitialized()
    /// @sa \c haveTransform()
    void getPowerSpectralDensity(int nFrequencies, T *psd[]) const;
    /// @result The power spectrum at each frequency.  If the input signal
    ///         signal has units of \f$ Volts \f$ then this has units of
    ///         \f$ Volts^2 \f$.
    /// @throws std::runtime_error if \c haveTransform() is false.
    std::vector<T> getPowerSpectrum() const;
    /// @brief Gets the power spectrum.
    /// @param[in] nFrequencies    The number of frequencies. This must match the
    ///                            result of \c getNumberOfFrequencies().
    /// @param[out] powerSpectrum  The power spectrum.  This is an array of
    ///                            of dimension [nFrequencies].  If the input
    ///                            signal has units of \f$ Volts \f$ then this
    ///                            has units of \f$ Volts^2 \f$.
    /// @throws std::invalid_argument if nFrequencies is invalid or 
    ///         powerSpectrum is NULL.
    /// @throws std::runtime_error if the class is not initialized or the
    ///         transform has not been computed.
    /// @sa \c getNumberOfFrequencies()
    /// @sa \c isInitialized()
    /// @sa \c haveTransform()
    void getPowerSpectrum(int nFrequencies, T *powerSpectrum[]) const;
    /// @}
private:
    class WelchImpl;
    std::unique_ptr<WelchImpl> pImpl;
};
}
#endif
