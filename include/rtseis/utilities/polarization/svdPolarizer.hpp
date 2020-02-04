#ifndef RTSEIS_UTILITIES_POLARIZATION_SVDPOLARIZER_HPP
#define RTSEIS_UTILITIES_POLARIZATION_SVDPOLARIZER_HPP
#include "rtseis/enums.h"
#include <memory>
namespace RTSeis::Utilities::Polarization
{
template<class T>
/*!
 * @brief Performs the recursive SVD polarization proposed by Rosenberger
 *        "Realtime ground-motion analysis: Distinguishing P and S arrivals
 *         in a noisy environment", 2010.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
class SVDPolarizer
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    SVDPolarizer(); 
    /*!
     * @brief Copy constructor.
     * @param[in] polarizer  Initializes this class from the polarizer class.
     */
    SVDPolarizer(const SVDPolarizer &polarizer);
    /*!
     * @brief Move constructor.
     * @param[in,out] polarizer  The polarizer class from which to initialize
     *                           this class.  On exit, polarizer's behavior
     *                           is undefined.
     */
    SVDPolarizer(SVDPolarizer &&polarizer) noexcept; 
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] polarizer   The polarizer class to copy to this.
     * @result A deep copy of the polarizer class.
     */
    SVDPolarizer& operator=(const SVDPolarizer &polarizer);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] polarizer  The polarizer class whose memory will be copied
     *                           to this.  On exit, polarizer's behavior is
     *                           undefined.
     * @result The memory from polarizer moved to this.
     */
    SVDPolarizer& operator=(SVDPolarizer &&polarizer) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~SVDPolarizer();
    /*!
     * @brief Releases all memory and resets the class.
     */
    void clear() noexcept;
    /*! @}*/

    /*!
     * @brief Initializes the SVD polarizer.
     * @param[in] decayFactor  This is a forgetting factor in the range of 
     *                         (0,1) for the polarizer.  For example, for a
     *                         P-wave one could compute the number of samples in
     *                         the P-wave window, \f$ N_p \f$, and set
     *                         decayFactor to \f$ \frac{N_p - 1}{N_p} \f$.
     *                         Similarly, for an S-wave, one could compute the
     *                         number of samples in the S-wave window,
     *                         \f$ N_s \f$, and set the decayFactor to
     *                         \f$ \frac{N_s 1}{N_s} \f$.
     * @param[in] noise        This is tolerance above which a the magnitude
     *                         of project vector must exceed to consistute
     *                         a rank update.  This can be a numerical
     *                         precision or instrument noise floor and has the
     *                         same units as the signals which are being polarized.
     * @param[in] mode         Distinguishes beteween real-time and post-processing.
     */
    void initialize(T decayFactor,
                    T noise,
                    RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING);
    /*!
     * @brief Initializes the SVD polarizer but here the noise floor is set
     *        to near epsilon.
     * @copydoc initialize
     */
    void initialize(T decayFactor,
                    RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING);
    /*! 
     * @brief Determines whether or not the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;

    /*!
     * @brief Computes the Karhunen-Loeve transform of the input data, the
     *        cosine of the incidence angle w.r.t. z, and the rectilinearity
     *        at each sample.
     * @param[in] npts  The number of samples.
     * @param[in] z     The trace on the vertical channel.  This is an array
     *                  whose dimension is [npts].
     * @param[in] n     The trace on the north channel.  This is an array whose
     *                  dimension is [npts].
     * @param[in] e     the trace on the east channel.  This is an array whose
     *                  dimension is [npts].
     * @param[out] klz  The Karhunen-Loeve transform of the vertical channel.
     *                  It is projected into the reduced rank subspace then
     *                  projected back to the ZNE frame.
     *                  This is an array whose dimension is [npts].
     * @param[out] kln  The Karhunen-Loeve transform of the north channel.
     *                  It is projected into the reduced rank subspace then
     *                  projected back to the ZNE frame.
     *                  This is an array whose dimension is [npts].
     * @param[out] kle  The Karhunen-Loeve transform of the east channel.
     *                  It is projected into the reduced rank subspace then
     *                  projected back to the ZNE frame.
     *                  This is an array whose dimension is [npts].
     * @param[out] cosIncidenceAngle  The cosine of the incidence angle w.r.t.
     *                                to the vertical channel.  This will be
     *                                close to 1 for for particle motion close
     *                                to the vertical which is indicative of a
     *                                P wave.  Conversely, this will be close
     *                                to 0 for an S-wave.  This is an array
     *                                whose dimension is [npts].
     * @param[out] rectilinearity  The rectilinearity.  This is defined as 
     *                             \f$ 1 - \frac{sigma_1}{sigma_0} \f$ where 
     *                             \f$ sigma_1 \f$ is the second largest
     *                             singular value and \f$ sigma_0 \f$ is the 
     *                             largest singular value.  This is an array
     *                             whose dimension is [npts].
     * @throws std::runtime_error if the class is not initialized.
     * @throws std::invalid_argument if any array is NULL.
     * @sa \c isInitialized()
     */
    void polarize(int npts, const T z[], const T n[], const T e[],
                  T *klz[], T *kln[], T *kle[],
                  T *cosIncidenceAngle[], T *rectilinearity[]);
    /*! 
     * @brief Computes the Karhunen-Loeve transform of the input data.
     * @param[in] npts  The number of samples.
     * @param[in] z     The trace on the vertical channel.  This is an array
     *                  whose dimension is [npts].
     * @param[in] n     The trace on the north channel.  This is an array whose
     *                  dimension is [npts].
     * @param[in] e     the trace on the east channel.  This is an array whose
     *                  dimension is [npts].
     * @param[out] klz  The Karhunen-Loeve transform of the vertical channel.
     *                  It is projected into the reduced rank subspace then
     *                  projected back to the ZNE frame.
     *                  This is an array whose dimension is [npts].
     * @param[out] kln  The Karhunen-Loeve transform of the north channel.
     *                  It is projected into the reduced rank subspace then
     *                  projected back to the ZNE frame.
     *                  This is an array whose dimension is [npts].
     * @param[out] kle  The Karhunen-Loeve transform of the east channel.
     *                  It is projected into the reduced rank subspace then
     *                  projected back to the ZNE frame.
     *                  This is an array whose dimension is [npts].
     * @throws std::runtime_error if the class is not initialized.
     * @throws std::invalid_argument if any array is NULL.
     * @sa \c isInitialized()
     */
    void polarize(int npts, const T z[], const T n[], const T e[],
                  T *klz[], T *kln[], T *kle[]);
    /*!
     * @brief Computes the cosine of the incidence angle w.r.t. z, and the
     *        rectilinearity at each sample.
     * @param[in] npts  The number of samples.
     * @param[in] z     The trace on the vertical channel.  This is an array
     *                  whose dimension is [npts].
     * @param[in] n     The trace on the north channel.  This is an array whose
     *                  dimension is [npts].
     * @param[in] e     the trace on the east channel.  This is an array whose
     *                  dimension is [npts].
     * @param[out] cosIncidenceAngle  The cosine of the incidence angle w.r.t.
     *                                to the vertical channel.  This will be
     *                                close to 1 for for particle motion close
     *                                to the vertical which is indicative of a
     *                                P wave.  Conversely, this will be close
     *                                to 0 for an S-wave.  This is an array
     *                                whose dimension is [npts].  This will be
     *                                bounded to [0,1].  If you want the 
     *                                angle of incidence then you must apply 
     *                                an arccos to this.
     * @param[out] rectilinearity  The rectilinearity.  This is defined as 
     *                             \f$ 1 - \frac{sigma_1}{sigma_0} \f$ where 
     *                             \f$ sigma_1 \f$ is the second largest
     *                             singular value and \f$ sigma_0 \f$ is the 
     *                             largest singular value.  This is an array
     *                             whose dimension is [npts].
     * @throws std::runtime_error if the class is not initialized.
     * @throws std::invalid_argument if any array is NULL.
     * @sa \c isInitialized()
     */
    void polarize(int npts, const T z[], const T n[], const T e[],
                  T *cosIncidenceAngle[], T *rectilinearity[]);
    /*!
     * @brief Resets the initial conditions on the source delay line
     *        to the default conditions or the initial conditions
     *        set when SVDPolarizer::setInitialConditions() was called.
     * @throws std::runtime_error if the class is not initialized.
     */
    void resetInitialConditions();
private:
    class SVDPolarizerImpl;
    std::unique_ptr<SVDPolarizerImpl> pImpl;
};

/*!
 * @brief Enhances a P-wave by modulating the signal using:
 *        \f$ I_n R_n x_n \f$ to the n'th sample.  Here, \f$ I_n \f$ is the
 *        cosine of the incidence angle, \f$ R_n \f$ is the rectilinearity,
 *        and \f$ x_n \f$ is the trace.  This can also be the KL-transformed
 *        trace.
 * @param[in] npts    The number of samples.
 * @param[in] z       The vertical trace to modulate.  This is an array whose
 *                    dimension is [npts].
 * @param[in] n       The north trace to modulate.  This is an array whose 
 *                    dimension is [npts].
 * @param[in] e       The east trace to modulate.  This is an array whose
 *                    dimension is [npts].
 * @param[in] cosIncidenceAngle  The cosine of the incidence angle at each
 *                               sample.  This is an array whose dimension
 *                               is [npts].
 * @param[in] rectlinearity      The rectilinearity at each sample.  This is
 *                               an array whose dimension is [nts]. 
 * @param[out] pz     The P-modulated vertical trace.  This is an array whose
 *                    dimension is [npts].
 * @param[out] pn     The P-modulated north trace.  This is an array whose
 *                    dimension is [npts].
 * @param[out] pe     The P-modulated east trace.  This is an array whose
 *                    dimension is [npts].
 * @throws std::invalid_argument of if any array is NULL.
 */
template<typename T>
void modulateP(int npts,
               const T z[], const T n[], const T e[],
               const T cosIncidenceAngle[], const T rectilinearity[],
               T *pz[], T *pn[], T *pe[]);
/*!
 * @brief Enhances the S-wave by modulating the signal using:
 *        \f$ (1 - I_n) R_n x_n \f$ to the n'th sample.  Here, \f$ I_n \f$ is
 *        the cosine of the incidence angle, \f$ R_n \f$ is the rectilinearity,
 *        and \f$ x_n \f$ is the trace.  This can also be the KL-transformed
 *        trace.
 * @param[in] npts    The number of samples.
 * @param[in] z       The vertical trace to modulate.  This is an array whose
 *                    dimension is [npts].
 * @param[in] n       The north trace to modulate.  This is an array whose 
 *                    dimension is [npts].
 * @param[in] e       The east trace to modulate.  This is an array whose
 *                    dimension is [npts].
 * @param[in] cosIncidenceAngle  The cosine of the incidence angle at each
 *                               sample.  This is an array whose dimension
 *                               is [npts].
 * @param[in] rectlinearity      The rectilinearity at each sample.  This is
 *                               an array whose dimension is [nts]. 
 * @param[out] sz     The S-modulated vertical trace.  This is an array whose
 *                    dimension is [npts].
 * @param[out] sn     The S-modulated north trace.  This is an array whose
 *                    dimension is [npts].
 * @param[out] se     The S-modulated east trace.  This is an array whose
 *                    dimension is [npts].
 * @throws std::invalid_argument of if any array is NULL.
 */
template<typename T>
void modulateS(int npts,
               const T z[], const T n[], const T e[],
               const T cosIncidenceAngle[], const T rectilinearity[],
               T *sz[], T *sn[], T *se[]);
}
#endif
