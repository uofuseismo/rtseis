#ifndef RTSEIS_UTILITIES_TRANSFORMS_WELCH_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_WELCH_HPP 1
#include <memory>

namespace RTSeis::Utilities::Transforms
{
/*!
 * @brief Estimates the power spectral density using Welch's method.
 *        This amounts to first dividing the data into overlapping segments.
 *        Next, a modified periodogram is computed in each segment. 
 *        Finally, the modified periodogram for all segments are averaged. 
 * @author Ben Baker, University of Utah
 * @copyright Ben Baker distributed under the MIT license.
 * @date July 2019
 */
class Welch
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    Welch();
    /*!
     * @brief Copy constructor.
     * @param[in] welch  The Welch class from which to initialize this class.
     */
    Welch(const Welch &welch);
    /*!
     * @brief Move constructor.
     * @param[in,out] welch  The Welch class from which to initialize this class.
     *                       On exit, welch's behavior is undefined.
     */
    Welch(Welch &&welch) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] welch  The class to copy.
     * @result A deep copy of the welch class.
     * @throws std::runtime_error if there is an initialization error.
     *         This is indicative of an internal error.
     */
    Welch& operator=(const Welch &welch);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] welch  The class to move to this.
     *                       On exit welch's behavior is undefined.
     * @result The memory that was moved from welch to this.
     */
    Welch& operator=(Welch &&welch) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~SlidingWindowRealDFT();
    /*!
     * @brief Releases memory on the class.
     */
    void clear() noexcept;
    /*! @} */ 

    /*!
     * @brief Initializes the Welch transform.
     * @param[in] nSamples            The number of samples in the signal.
     *                                This must be positive.
     * @param[in] nSamplesPerSegment  The number of samples in each DFT segment.
     * @param[in] dftLength           The length of the DFT.  If the number of
     *                                samples in each segment results in an
     *                                inefficient transform length then this
     *                                will zero pad each window to the desired
     *                                DFT length.  This must be greater than or
     *                                equal to nSamplesPerSegment.
     * @param[in] nSamplesInOverlap   As the window slides along this defines
     *                                the number of samples overlapping between
     *                                the current and prior window.  This must
     *                                be in the range [0, nSamplesPerSegment-1].
     *                                For a Hanning window nSamplesInOverlap 
     *                                being equal to half nSamplesPerSegment
     *                                represents a reasonable tradeoff between
     *                                accureately estimating the signal power
     *                                while not over counting the data. 
     *                                Narrower windows may require this to 
     *                                correspond to a greater percentage of
     *                                nSamplesInSegment.
     *                                In the limit, if this is 0 then this
     *                                method is equivalent to Bartlett's method.
     * @param[in] windowLength        The length of the window.  If 0 then no
     *                                window function will be applied.
     *                                Otherwise, this must equal
     *                                nSamplesPerSegment.
     * @param[in] window              The window to apply to each segment.
     *                                If windowLength is 0 then this is ignored.
     *                                Otherwise, it is an array of dimension
     *                                [windowLength].
     * @param[in] samplingRate        The sampling rate in Hz.
     * @param[in] detrendType         Defines the detrend strategy to be
     *                                applied to each signal prior to windowing.
     * @param[in] precision           The precision of the underlying DFT.
     * @throws std::invalid_argument if any parameters are incorrect.
     */
    void initialize(const int nSamples,
                    const int nSamplesPerSegment,
                    const int dftLength,
                    const int nSamplesInOverlap,
                    const int windowLength,
                    const double window[],
                    const double samplingRate = 1.0,
                    const SlidingWindowDetrendType detrendType,
                    const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Flag indicating whether or not the class is initialized.
     * @result True indicates that the class is inititalized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Returns the number of frequencies.
     * @result The number of frequencies.
     * @throws std::runtime_error if the class is not intitialized.
     */
    int getNumberOfFrequencies() const;
    /*!
     * @brief Gets the frequencies at which the power spectral density was
     *        estimated. 
     * @param[in] nFrequencies  The number of frequencies.  This must match the
     *                          result of \c getNumberOfFrequencies().
     * @param[out] frequencies  The frequencies (Hz) at which the spectral
     *                          density was estimated.  This is an array of
     *                          dimensino [nFrequencies].
     * @throws std::invalid_argument if nFrequencies is invalid or frequencies
     *         is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c getNumberOfFrequencies()
     * @sa \c isInitialized()
     */
    void getFrequencies(const int nFrequencies, double *frequencies[]);
    /*!
     * @brief Gets the power spectral density estimate.
     * @param[in] nFrequencies  The number of frequencies. This must match the
     *                          result of \c getNumberOfFrequencies().
     * @param[out] psd          The power spectral density.  This is an array
     *                          of dimension [nFrequencies].  If the input 
     *                          signal has units of \f$ Volts \f$ then this
     *                          has units of \f$ \frac{Volts^2}{Hz} \f$.
     * @throws std::invalid_argument if nFrequencies is invalid or psd is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c getNumberOfFrequencies()
     * @sa \c isInitialized()
     */
    void getPowerSpectralDensity(const int nFrequencies, double *psd[]);
    /*!
     * @brief Gets the power spectrum.
     * @param[in] nFrequencies    The number of frequencies. This must match the
     *                            result of \c getNumberOfFrequencies().
     * @param[out] powerSpectrum  The power spectrum.  This is an array of
     *                            of dimension [nFrequencies].  If the input
     *                            signal has units of \f$ Volts \f$ then this
     *                            has units of \f$ Volts^2 \f$.
     * @throws std::invalid_argument if nFrequencies is invalid or 
     *         powerSpectrum is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c getNumberOfFrequencies()
     * @sa \c isInitialized()
     */
    void getPowerSpectrum(const int nFrequencies, double *powerSpectrum[]);
private:
    class WelchImpl;
    std::unique_ptr<WelchImpl> pImpl;
};
}
#endif
