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
     /*!
      * @brief Defines the polarization type.
      */
     enum PolarizationType
     {
        P_POLARIZE, /*!< Denotes a P-polarization. */
        S_POLARIZE  /*!< Denotes an S-polarization. */
     };
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor
     */
    SVDPolarizer(); 
    /*! @} */

    /*! @name Destructors
     * @{
     */
    ~SVDPolarizer();
    /*! @}*/

    /*!
     * @brief Initializes the SVD polarizer.
     * @param[in] pDecayFactor  This is a forgetting factor in the range of 
     *                          (0,1) for the P-wave polarizer.  For example,
     *                          one could compute the number of samples in the
     *                          P-wave window, Np, and set this to (Np - 1)/Np.
     * @param[in] sDecayFactor  This is a forgetting factor in the range of
     *                          (0,1) for the S-wave polarizer.  For example,
     *                          one could compute the number of samples in the
     *                          S-wave window, Ns, and set this to (Ns - 1)/Ns.
     */
    void initialize(const double decayFactor,
                    const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING);
    /*!
     * @brief Computes the incidence angle w.r.t. z, the rectilinearity,
     *        and rotated data at each sample.
     * @param[in] n     The number of samples.
     * @param[in] zne   The multiplexed 3-component seismogram.  This is a 
     *                  [n x 3] array in row major format.  Each row corresponds
     *                  is a length 3 tuple giving the response on the
     *                  (vertical, north, east) channels.
     * @param[out] cosIncidenceAngle  The cosine of the incidence angle w.r.t.
     *                                to the vertical channel.  This will be
     *                                close to 1 for for particle motion close
     *                                to the vertical which is indicative of a
     *                                P wave.  Conversely, this will be close
     *                                to 0 for an S-wave.  This is an array
     *                                whose dimension is [n].
     * @param[out] rectilinearity  The rectilinearity.  This is defined as 
     *                             \f$ 1 - \frac{sigma_1}{sigma_0} \f$ where 
     *                             \f$ sigma_1 \f$ is the second largest
     *                             singular value and \f$ sigma_0 \f$ is the 
     *                             largest singular value.  This is an array
     *                             whose dimension is [n].
     * @throws std::invalid_argument if n is positive and zne is NULL. 
     */
    void polarize(const int n, const T zne[],
                  T incidenceAngle[], T rectilinearity[], T r[]);
    /*!
     * @brief Polarizes the 3-component seismgoram.
     * @param[in] n    The number of input samples.
     * @param[in] z1   The vertical channel'response.  This is an array of dimension [n].
     * @param[in] n2   The north channel's response.  This is an array of dimension [n].
     * @param[in] e3   The east channel's response.  This is an array of dimension [n].
     * @param[out] rz  The P-polarized vertical channel. This is an array of dimension [n].
     * @param[out] rn  The P-polarized north channel.  This is an array of dimension [n].
     * @param[out] re  The P-polarized east channel.  This is an array of dimension [n].
     * @throws std::runtime_error if the class is not inititalized.
     */
    void polarize(const int n,
                  const T z1[], const T n2[], const T e3[],
                  T rz[], T rn[], T re[]);
//                  SVDPolarizer::PolarizationType type);
private:
    class SVDPolarizerImpl;
    std::unique_ptr<SVDPolarizerImpl> pImpl;
};

template<typename T>
void pPolarize(const int n, const T zne[],
               const T inc[], const T rect[], const T r[],
               T pPolarized[]);
template<typename T>
void sPolarize(const int n, const T zne[],
               const T inc[], const T rect[], const T r[],
               T sPolarized[]);

}
#endif
