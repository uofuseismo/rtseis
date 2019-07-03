#ifndef RTSEIS_UTILITIES_TRANSFORMS_SLIDINGWINDOWREALDFT_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_SLIDINWINDOWGREALDFT_HPP 1
#include <memory>
#include <complex>
#include "rtseis/enums.h"
#include "rtseis/utilities/transforms/enums.hpp"

namespace RTSeis::Utilities::Transforms
{
/*!
 * @brief This is a utility class that slides a window over a real time series
 *        and computes the DFT in each window.  This is the basis of many 
 *        higher order spectral analysis methods such as the STFT, spectrogram,
 *        and Welch's method.  
 * @note This class is not intended for use by a large audience.  Because it is
 *       foundational to higher-order methods whose results are more desirable
 *       this interface emphasizes efficiency at the sake of clarity.
 * @author Ben Baker, University of Utah
 * @copyright Ben Baker distributed under the MIT license.
 * @date May 2019.
 */
class SlidingWindowRealDFT
{
public:
    /*! @name Constructors
     * @{
     */
    /*! 
     * @brief Default constructor.
     */
    SlidingWindowRealDFT();
    /*!
     * @brief Copy constructor.
     * @param[in] swdft  The sliding window real DFT class to initialize from.
     * @throws std::runtime_error if the class fails to initialize. 
     *         This is indicative of an internal error. 
     */
    SlidingWindowRealDFT(const SlidingWindowRealDFT &swdft);
    /*!
     * @brief Move constructor.
     * @param[in,out] swdft  The sliding window real DFT class to initialize
     *                       from.  On exit, swdft's behavior will be undefined.
     */
    SlidingWindowRealDFT(SlidingWindowRealDFT &&swdft) noexcept;
    /*! @} */
 
    
    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] swdft  The class to copy.
     * @result A deep copy of the sliding window real DFT class.
     * @throws std::runtime_error if the class fails to initialize.
     *         This is indicative of an internal error.
     */
    SlidingWindowRealDFT& operator=(const SlidingWindowRealDFT &swdft);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] swdft  The class to move.  On exit swdft's behavior is
     *                       undefined.
     * @result The memory that was moved from swdft to this.
     */
    SlidingWindowRealDFT& operator=(SlidingWindowRealDFT &&swdft) noexcept;
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
     * @brief Initializes the sliding-window real discrete Fourier Transform.
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
     * @param[in] windowLength        The length of the window.  If 0 then no
     *                                window function will be applied.
     *                                Otherwise, this must equal 
     *                                nSamplesPerSegment.
     * @param[in] window              The window to apply to each segment.
     *                                If windowLength is 0 then this is ignored.
     *                                Otherwise, it is an array of dimension
     *                                [windowLength].
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
     * @brief Returns the number of sliding time windows for which a
     *        DFT was computed.  This is the number of columns in the
     *        output matrix.
     * @result The number of windows.
     * @throws std::runtime_error if the class is not intitialized.
     */
    int getNumberOfTransformWindows() const;
    /*!
     * @brief Returns the expected number of samples in the time series.
     * @result The expected number of samples.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getNumberOfSamples() const;
    /*!
     * @brief Gets the precision of the underlying Fourier transform.
     * @result The precision of the underlying transform.
     * @throws std::runtime_error if the class is not initialized.
     */
    RTSeis::Precision getPrecision() const;

    /*!
     * @brief Computes the sliding window DFT of the real signal.
     * @param[in] nSamples  The number of samples in the signal.
     *                      This must match the result of 
     *                      \c getNumberOfSamples().
     * @param[in] x         The signal to transform.  This is an array whose
     *                      dimension is [nSamples].
     * @throws std::invalid_argument if any arguments are invalid.
     * @throws std::runtime_error if the class is not initalized.
     * @sa \c getNumberOfTransformWindow(), \c getNumberOfFrequencies()
     */
    void transform(const int nSamples, const double x[]);

    /*!
     * @brief Gets a pointer to the transform in the iWindow'th window.
     * @param[in] iWindow  The window of the given transform.  This must
     *                     be in the range 
     *                     [0, \c getNumberOfTransformWidnwos()-1]. 
     * @throws std::invalid_argument if iWindow is out of bounds.
     * @throws std::runtime_error if the class precision is FLOAT 
     *         or the \c transform() has not yet been called. 
     * @sa \c transform(), \c getNumberOfTransformWindows(), \c getPrecision()
     */
    const std::complex<double> *getTransform64f(const int iWindow) const;
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
    const std::complex<float> *getTransform32f(const int iWindow) const;
private:
    class SlidingWindowRealDFTImpl;
    std::unique_ptr<SlidingWindowRealDFTImpl> pImpl;
};
}
#endif
