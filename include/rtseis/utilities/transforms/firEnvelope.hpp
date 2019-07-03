#ifndef RTSEIS_UTILITIES_TRANSFORMS_FIRENVELOPE_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_FIRENVELOPE_HPP 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis::Utilities::Transforms
{
/*!
 * @class FIREnvelope firEnvelope.hpp "include/rtseis/utilities/transforms/firEnvelope.hpp"
 * @brief Computes the envelope of a signal using.  This works by designing an
 *        FIR Hilbert transfomer.
 * @note The post-processing modes and real-time modes yield different results
 *       because the post-processing mode will accomodate the mean of the data
 *       and compensate for the phase delay of the FIR filter.  Moreover,
 *       if running in real-time it is important to remove the mean by first
 *       applying a high-pass filter.
 * @sa Envelope
 */
class FIREnvelope
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    FIREnvelope();
    /*!
     * @brief Copy constructor.
     * @param[in] firEnvelope  Class from which to initialize this class.
     */
    FIREnvelope(const FIREnvelope &firEnvelope);
    /*!
     * @brief Move constructor. 
     * @param[in] firEnvelope  Class from which to initialize this class.
     *                         On exit firEnvelope's behavior will be undefined.
     */
    FIREnvelope(FIREnvelope &&firEnvelope) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] firEnvelope  FIREnvelope class to copy.
     * @result A deep copy of the firEnvelope.
     */
    FIREnvelope& operator=(const FIREnvelope &firEnvelope);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] firEnvelope  FIR envelope class whose memory is to be
     *                             moved to this class.
     *                             On exit firEnvelope's behavior is undefined.
     * @result Contains the moved memory from firEnvelope.
     */
    FIREnvelope& operator=(FIREnvelope &&firEnvelope) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~FIREnvelope();
    /*!
     * @brief Releases all memory on the class and resets all parameters.
     *        The class must be re-initialized prior to usage.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Checks if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;

    /*!
     * @brief Initializes the FIR-based envelope.
     * @param[in] ntaps      The number of filter taps.  This must be positive.
     *                       If ntaps is odd then a Type III FIR filter will be 
     *                       used.  This can be computationally advantageous
     *                       as the real-part of the FIR filter is simply a 
     *                       unit delay of ntaps/2 samples.  However, the filter
     *                       response at the Nyquist frequency will be 0.
     *                       If ntaps is even then a Type IV FIR filter will be
     *                       used.  This is less computationally efficient as
     *                       both the real and imaginary parts of the Hilbert
     *                       transform filter are to be applied.  However,
     *                       this will better characterize the transformer
     *                       response at the Nyquist frequency. 
     * @param[in] mode       The processing mode.
     * @param[in] precision  Controls the underlying precision of the envelope
     *                       calculation.
     * @throws std::invalid_argument if ntaps is not positive.
     */
    void initialize(const int ntaps,
                    const RTSeis::ProcessingMode mode=RTSeis::POST_PROCESSING,
                    const RTSeis::Precision precision=RTSeis::DOUBLE);

    /*!
     * @brief Gets the length of the initital condition array.
     * @result The length of the initial condition array.
     * @throws std::runtime_error if the the class is not inititalized.
     */
    int getInitialConditionLength() const;
    /*!
     * @brief Sets the initial conditions.
     * @param[in] nz   The length of the initial conditions.
     *                 This must equal \c getInitialConditionLength().
     * @param[in] zi   The initial conditions.  This is an array of
     *                 dimension [nz].
     * @throws std::invalid_argument if the nz is not correct of zi is NULL. 
     * @throws std::runtime_error if the class is not inititalized.
     */
    void setInitialConditions(const int nz, const double zi[]);

    /*!
     * @brief Computes the envelope of the signal.
     * @param[in] n   The number of samples in the signal.
     * @param[in] x   The signal of which to take the envelope.  This is
     *                an array whose dimension is [n].
     * @param[out] y  The upper envelope of the x.  This is an array 
     *                whose dimension is [n].
     * @throws std::invalid_argument if x or y are NULL.
     * @throws std::runtime_error if the class is not inititalized.
     */
    void transform(const int n, const double x[], double y[]);
    /*! @copydoc transform */
    void transform(const int n, const float x[], float y[]);
    /*!
     * @brief Resets the initial conditions of the underlying FIR filters.
     *        This may be useful after a gap is encountered.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void resetInitialConditions();
private:
    class FIREnvelopeImpl;
    std::unique_ptr<FIREnvelopeImpl> pImpl;
};
} // RTSeis
#endif
