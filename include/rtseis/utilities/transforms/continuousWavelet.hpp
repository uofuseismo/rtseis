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
 * @brief Computes the continuous wavelet transform of a signal
 *        \f[
 *           \frac{1}{\sqrt{|a|}} 
 *           \int s(t) \psi^* \left ( \frac{t}{a} \right ) \, dt
 *        \f]
 *        where \f$ a \f$ is the scale.
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

    /*! @name Initialization
     * @{
     */
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
    /*! @} */

    /*! @name Transform
     * @{
     */
    /*!
     * @brief Computes the continuous wavelet transform of the given signal. 
     * @param[in] n   The number of samples in x.  This must match
     *                \c getNumberOfSamples().
     * @param[in] x   The signal to transform.  This is an array whose dimension
     *                is [n].
     * @throws std::invalid_argument if n is the wrong size or x is NULL.
     * @throws std::runtime_error if \c isInitialized() is false. 
     */
    void transform(int n, const T x[]);
    /*! @} */

    /*! @name Results
     * @{
     */
    /*!
     * @result True indicates the continuous wavelet transform has been
     *         computed.
     */
    [[nodiscard]] bool haveTransform() const noexcept;
    /*!
     * @brief Gets the CWT.  
     * @param[in] nSamples  The number of samples in the CWT. 
     *                      This must match \c getNumberOfSamples().
     * @param[in] nScales   The number of scales in the CWT.
     *                      This must match \c getNumberOfScales().
     * @param[out] cwt      The CWT.  This is an [nScales x nSamples] matrix
     *                      stored in row major order.  cwt[0,:] corresponds
     *                      to the first scale.
     * @throws std::runtime_error if \c haveTransform() is false.
     * @throws std::invalid_argument if nSamples or nScales is the wrong size
     *         or cwt is NULL.
     */
    void getTransform(int nSamples, int nScales, std::complex<T> *cwt[]) const;
    /*!
     * @brief Gets the amplitude of the complex CWT.
     * @param[in] nSamples  The number of samples in the CWT.
     *                      This must match \c getNumberOfSamples().
     * @param[in] nScales   The number of scales in the CWT.
     *                      This must match \c getNumberOfScales().
     * @param[out] amplitude   The amplitude of the CWT.  This is an
     *                         [nScales x nSamples] matrix stored in row major
     *                         order.
     * @throws std::runtime_error if \c haveTransform() is false.
     * @throws std::invalid_argument if nSamples or nScales is the wrong size
     *         or amplitude is NULL.
     * @sa \c getTransform()
     */
    void getAmplitudeTransform(int nSamples, int nScales, T *amplitude[]) const;
    /*!
     * @brief Gets the phase of the complex CWT.
     * @param[in] nSamples  The number of samples in the CWT.
     *                      This must match \c getNumberOfSamples().
     * @param[in] nScales   The number of scales in the CWT.
     *                      This must match \c getNumberOfScales().
     * @param[out] phase    The phase angle of the CWT in radians.  This is an
     *                      [nScales x nSamples] matrix stored in row major
     *                      order.
     * @throws std::runtime_error if \c haveTransform() is false.
     * @throws std::invalid_argument if nSamples or nScales is the wrong size
     *         or phase is NULL.
     * @sa \c getTransform()
     */
    void getPhaseTransform(int nSamples, int nScales, T *phase[]) const;
    /*! @} */
private:
    class ContinuousWaveletImpl;
    std::unique_ptr<ContinuousWaveletImpl> pImpl;
};
}
#endif
