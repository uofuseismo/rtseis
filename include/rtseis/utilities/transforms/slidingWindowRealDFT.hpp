#ifndef RTSEIS_UTILITIES_TRANSFORMS_SLIDINGWINDOWREALDFT_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_SLIDINWINDOWGREALDFT_HPP 1
#include <memory>
#include <complex>
#include "rtseis/enums.h"
#include "rtseis/utilities/transforms/enums.hpp"

namespace RTSeis
{
namespace Utilities
{
namespace Transforms
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
    /*! 
     * @brief Default constructor.
     */
    SlidingWindowRealDFT();

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
     *                                will zero pad zero pad each window to 
     *                                the desired DFT length.  This must be
     *                                greater than or equal to
     *                                nSamplesPerSegment.
     * @param[in]
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
     * @brief Returns the number of time samples in the sliding window DFT.
     *        This is the number of columns in the output matrix.
     * @result The number of time samples or, equivalently, columsn in the
     *         output matrix.
     * @throws std::runtime_error if the class is not intitialized.
     */
    int getNumberOfTimeSamples() const;
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
     * @param[in] nWork     The workspace allocated to y.  This must be at
     *                      least \c getNumberOfFrequencies() x 
     *                            \c getNumberOfTimeSamples().
     * @param[out] y        This is a row-major matrix containing the
     *                      DFT's in each window.  Here, the rows correspond
     *                      to frequency and the columns correspond to time.
     *                      While this has dimension [nwork] only the first
     *                      [\c getNumberOfFrequencies() x
     *                       \c getNumberOfTimeSamples()] samples are valid.
     * @throws std::invalid_argument if any arguments are invalid.
     * @throws std::runtime_error if the class is not initalized.
     * @sa \c getNumberOfTimeSamples(), \c getNumberOfFrequencies()
     */
    void transform(const int nSamples, const double x[]);

    /*!
     * @brief Gets a pointer to the transform in the iWindow'th window.
     * @param[in] iWindow  The window of the given transform.  This must
     *                     be in the range [0, \c getNumberOfColumns()-1]. 
     * @throws std::invalid_argument if the class precision is DOUBLE or
     *         iWindow is out of bounds.
     * @sa \c getNumberOfColumns(), \c getPrecision()
     */
    const std::complex<double> *getTransform64f(const int iWindow) const;
    /*! 
     * @brief Gets a pointer to the transform in the iWindow'th window.
     * @param[in] iWindow  The window of the given transform.  This must
     *                     be in the range [0, \c getNumberOfColumns()-1]. 
     * @throws std::invalid_argument if the class precision is DOUBLE or
     *         iWindow is out of bounds.
     * @sa \c getNumberOfColumns(), \c getPrecision()
     */
    const std::complex<float> *getTransform32f(const int iWindow) const;
private:
    class SlidingWindowRealDFTImpl;
    std::unique_ptr<SlidingWindowRealDFTImpl> pImpl;
};

}
}
}

#endif
