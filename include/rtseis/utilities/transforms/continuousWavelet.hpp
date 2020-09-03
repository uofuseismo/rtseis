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
     * @param[in] nSamples        The number of samples in the signal to
     *                            transform.
     * @param[in] wavelet         The wavelet to evaluate.
     *                            \c wavelet.isInitialized() must be true. 
     * @param[in] scales          The scales at which to evaluate the CWT.
     * @param[in] samplingPeriod  The signal sampling rate in seconds.
     *                            Note, this will override the sampling period
     *                            set in wavelet.
     * @throws std::invalid_argument nSamples is less than 1, samplingPeriod is not positive,
     */
    void initialize(int nSamples,
                    const Wavelets::IContinuousWavelet &wavelet,
                    std::vector<double> &scales,
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
     * @result The e-folding factors (in Hz) at each sample.
     *         This has dimension [\c getNumberOfSamples].
     * @throws std::runtime_error if the class is not initialized.
     */
    [[nodiscard]] std::vector<T> getConeOfInfluence() const;
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
