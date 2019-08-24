#ifndef RTSEIS_UTILITIES_POLARIZATION_EIGENPOLARIZER_HPP
#define RTSEIS_UTILITIES_POLARIZATION_EIGENPOLARIZER_HPP
#include "rtseis/enums.h"
#include <memory>
namespace RTSeis::Utilities::Polarization
{
/*!
 * @brief Polarizes a three-component seismogram using the method described
 *        by Jurkevics, 1988.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T = double>
class EigenPolarizer
{
public:
    /*! @name Constructors
     *  @{
     */
    /*!
     * @brief Constructor
     */
    EigenPolarizer();
    /*!
     * @brief Copy constructor.
     * @param[in] polarizer  The polarization class from which to
     *                       initialize this class.
     */
    EigenPolarizer(const EigenPolarizer &polarizer);
    /*!
     * @brief Move constructor.
     * @param[in,out] polarizer  The polarization class fromw which to 
     *                           initialize this class.  On exit, polarizer's
     *                           behavior is undefined.
     */
    EigenPolarizer(EigenPolarizer &&polarizer) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] polarizer  The polarizer to copy.
     * @result A deep copy of the input eigen-polarizer.
     */
    EigenPolarizer& operator=(const EigenPolarizer &polarizer);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] polarizer  The polarizer whose memory will be moved to
     *                           this.  On exit, polarizer's behavior is
     *                           undefined.
     * @result The memory from polarizer moved to this.
     */
    EigenPolarizer& operator=(EigenPolarizer &&polarizer) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~EigenPolarizer();
    /*!
     * @brief Clears all memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the polarization class.
     * @param[in] nSamples  The number of samples expected in each signal.
     * @throws std::invalid_argument if nSamples is not positive.
     */
    void initialize(const int nSamples);

    /*!
     * @brief Sets the input signals.
     * @param[in] nSamples  The number of samples in the input signal.
     *                      This must match \c getNumberOfSamples().
     * @param[in] vertical  The vertical response where $+Z$ is positive up.
     *                      This is an array whose dimension is [nSamples].
     * @param[in] north     The north response where $+N$ is positive north.
     *                      This is an array whose dimension is [nSamples].
     * @param[in] east      The east component where $+E$ is postiive east.
     *                      This is an array whose dimension is [nSamples]. 
     * @throws std::invalid_argument if nSamples is inconsistent or the
     *         vertical, north, or east is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized(), \c getNumberOfSamples()
     */
    void setSignals(const int nSamples,
                    const T vertical[],
                    const T north[],
                    const T east[]);

    /*!
     * @brief Gets the radial and transverse signals.  These are rotated
     *        around the back-azimuth estimated from the eigendecomposition
     *        of the cross-variance matrix.
     * @param[in] nSamples     The number of samples.  This must equal the
     *                         expected number of samples in
     *                         \c getNumberOfSamples().
     * @param[out] radial      The radial channel.  This is an array whose
     *                         dimension is [nSamples].
     * @param[out] transverse  The transverse channel.  This is an array
     *                          whose dimension [nSamples].
     * @throws std::invalid_argument if nSamples or radial or transverse
     *         is NULL.
     * @throws std::runtime_error if the input signals were not set.
     * @sa \c getBackAzimuth(), \c haveInputSignals(), \c \getNumberOfSamples()
     */
    void getRadialTransverseSignals(const int nSamples,
                                    T *radial[],
                                    T *transverse[]) const;
    /*!
     * @brief Gets the longitudinal, radial, and transverse signals.  These are
     *        rotated around the incidence angle and back-azimuth estimated
     *        from the eigendecomposition of the cross-variance matrix.
     * @param[in] nSamples       The number of samples.  This must equal the
     *                           expected number of samples in
     *                           \c getNumberOfSamples().
     * @param[out] longitudinal  The longitudinal channel.  This is an array
     *                           whose dimension is [nSamples].
     * @param[out] radial        The radial channel.  This is an array whose
     *                           dimension is [nSamples].
     * @param[out] transverse    The transverse channel.  This is an array
     *                           whose dimension [nSamples].
     * @throws std::invalid_argument if nSamples is invalid or longitudinal,
     *         radial, or transverse is NULL.
     * @throws std::runtime_error if the input signals were not set.
     * @sa \c getBackAzimuth(), \c haveInputSignals(), \c \getNumberOfSamples()
     */
    void getLongitudinalRadialTransverseSignals(const int nSamples,
                                                T *longitudinal[],
                                                T *radial[],
                                                T *transverse[]) const; 
    /*! @name Derived Metrics
     * @{
     */
    /*!
     * @brief Gets the rectilinearity:
     *        \f[
     *           1 - \frac{\lambda_2 + \lambda_3}/{2 \lambda_1},
     *        \f] 
     *        For pure body waves this will be close to $1$. Here,
     *        $\lambda_1$ is the largest eigenvalue and $lambda_3$
     *        is the smallest eigenvalue. 
     * @result The rectilinearity which is in the range [0,1].
     * @throws std::runtime_error if the input signals are not set.
     * @sa \c haveInputSignals() 
     */
    T getRectilinearity() const;
    /*!
     * @brief Gets the apparent source-to-receiver azimuth.  
     * @param[in] wantRadians  If true then the result is in radians.
     *                         Otherwise, the result is in degrees.
     * @result The apparent azimuth measured positive east from north.  This
     *         will be in the range \f$ [0,2*\pi) \f$ if wantRadians is true
     *         (default) or \f$ [0,360) \f$ if wantRadians is false.
     * @throws std::runtime_error if the is not initialized.
     * @details The azimuth is calculated from the largest eigenvector by
     *          computing \f$ \frac{\pi}{2} - \atan2(u_N, u_E) \f$.
     *          where \f$ u_N \f$ is the north component of the eigenvector
     *          and \f$ u_E \f$ the east component.
     * @sa \c haveInputSignals()
     */
    T getAzimuth(const bool wantRadians = true) const;
    /*!
     * @brief Gets the apparent receiver-to-source azimuth, i.e., the
     *        backazimuth.
     * @param[in] wantRadians  If true then the result is in radians.
     *                         Otherwise, the result is in degrees.
     * @result The apparent backazimuth measured positive east from north.  This
     *         will be in the range \f$ [0,2*\pi) \f$ if wantRadians is true
     *         (default) or \f$ [0,360) \f$ if wantRadians is false.
     * @throws std::runtime_error if the is not initialized.
     * @sa \c haveInputSignals(), \c getAzimuth()
     */
    T getBackAzimuth(const bool wantRadians = true) const;
    /*!
     * @brief Gets the apparent incidence angle.
     * @result The apparent incidence angle.  This will be in the range
     *         \f$ [0, \pi/2] \f$ if wantRadians is true (default)
     *         or \f$ [0,90] \f$ if wantRadians is false. 
     * @throws std::runtime_error if the is not initialized.
     * @details The incidence angle is calculated from the largest eigenvector
     *          by computing \f$ \acos(|u_Z|) \f$ where \f$ |u_Z| \f$ is the
     *          the magnitude of the z component of the eigenvector.
     * @note There is an ambiguity for downhole instruments in that the 
     *       incidence angle is confined to [0,90] degrees.
     */
    T getIncidenceAngle(const bool wantRadians = true) const;
    /*! @} */
    /*! @name Properties
     * @{
     */
    /*!
     * @brief Gets the number of samples in the input signals.
     * @throws std::runtime_error if the class is not yet initialized.
     * @sa \c isInitialized()
     */
    int getNumberOfSamples() const;
    /*!
     * @brief Flag indicating that the input signals have been set.
     * @result True indicates that the input signals have been set.
     */
    bool haveInputSignals() const noexcept;
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;
    /*! @} */
private:
    class EigenPolarizerImpl;
    std::unique_ptr<EigenPolarizerImpl> pImpl;
};
}
#endif
