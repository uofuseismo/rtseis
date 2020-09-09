#ifndef RTSEIS_UTILITIES_TRANSFORMS_CONTINUOUSWAVELET_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_CONTINUOUSWAVELET_HPP 1
#include <memory>
#include "rtseis/enums.h"
namespace RTSeis::Utilities::Transforms
{
namespace Wavelets
{
class IContinuousWavelet;
}
/*!
 * @class ContinuousWavelet continuousWavelet.hpp "include/rtseis/utilities/transforms/continuousWavelet.hpp"
 * @brief Computes the continuous wavelet transform of a signal.  
 * @ingroup rtseis_utils_transforms
 * @sa Hilbert
 */
template<class T = double>
class ContinuousWavelet
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    ContinuousWavelet();
    /*!
     * @brief Copy constructor.
     * @param[in] cwt   The continuous wavelet transform class from which to 
     *                  initialize this class.
     */ 
    ContinuousWavelet(const ContinuousWavelet &cwt);
    /*!
     * @brief Move constructor.
     * @param[in,out] cwt  The continous wavelet transform class from which to
     *                     initialize this class.  On exit, cwt's behavior is
     *                     undefined.
     */
    ContinuousWavelet(ContinuousWavelet &&cwt) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] cwt   The continuous wavelet transform class to copy to this.
     * @result A deep copy of the continuous wavelet transform.
     */
    ContinuousWavelet& operator=(const ContinuousWavelet &cwt);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] cwt  The continuous wavelet transform class whose memory
     *                     will be moved to this.  On exit, cwt's behavior is
     *                     undefined.
     * @result The memory from cwt moved to this.
     */
    ContinuousWavelet& operator=(ContinuousWavelet &&cwt) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~ContinuousWavelet();
    /*!
     * @brief Releases all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the 
     * @param[in] nSamples   The number of samples in the signal to transform.
     * @param[in] nScales    The number of scales.  This is the scaleogram 
     *                       equivalent to the number of frequencies in an
     *                       spectrogram.
     * @param[in] scales     These are the dimensionless scales.  This is an 
     *                       array of dimension [n].  Note, if you are using
     *                       the Morlet wavelet then scale, s, is related to
     *                       frequency, f in Hz, by
     *                       \f$ s = \frac{\omega_0 f_s}{2 \pi f} \f$
     *                       where \f$ f_s \f$ is the sampling rate in Hz.
     * @param[in] wavelet    The wavelet to evaluate.  This must have an 
     *                       evaluate method which, for a given scale,
     *                       returns the wavelet's nSamples time domain samples.
     * @param[in] samplingRate  The sampling rate in Hz.
     * @throws std::invalid_argument nSamples is less than 1,
     *         nScales is less than 1, samplingPeriod is not positive,
     *         or any scale is not positive.
     */
    void initialize(int nSamples,
                    int nScales, const double scales[],
                    const Wavelets::IContinuousWavelet &wavelet,
                    double samplingRate = 1);

    /*!
     * @result The number of samples in the CWT.
     * @throws std::runtime_error if \c isInitialized() is false.
     */
    [[nodiscard]] int getNumberOfSamples() const;
    /*!
     * @result The number of scales in the CWT.
     * @throws std::runtime_error if \c isInitialized() is false.
     */
    [[nodiscard]] int getNumberOfScales() const;
    /*!
     * @result The sampling rate in Hz.
     * @throws std::runtime_error if \c isInitialized() is false.
     */
    [[nodiscard]] double getSamplingRate() const;
    /*!
     * @result The e-folding factors (in Hz) at each sample.
     *         This has dimension [\c getNumberOfSamples].
     * @throws std::runtime_error if the class is not initialized.
     */
    //[[nodiscard]] std::vector<T> getConeOfInfluence() const;
    /*!
     * @result True indicates that the class is initialized.
     */
    [[nodiscard]] bool isInitialized() const noexcept;

    void transform(int n, const T x[]);
private:
    class ContinuousWaveletImpl;
    std::unique_ptr<ContinuousWaveletImpl> pImpl;
};
}
#endif
