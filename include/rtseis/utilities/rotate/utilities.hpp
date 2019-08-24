#ifndef RTSEIS_UTILITIES_ROTATE_UTILITIES_HPP
#define RTSEIS_UTILITIES_ROTATE_UTILITIES_HPP
/*!
 * @brief Utilities for rotating seismograms.
 * @author Ben Baker (University of Utah) distributed under the MIT license.
 */
namespace RTSeis::Utilities::Rotate
{
/*!
 * @brief Rotates a (north,east) channel pair to (radial,transverse).
 * @param[in] nSamples     The number of samples in the seismograms.
 * @param[in] backAzimuth  The receiver to source azimuth in radians.
 * @param[in] north        The north channel.  This is an array whose dimension
 *                         is [nSamples].
 * @param[in] east         The east channel.  This is an array whose dimension
 *                         is [nSamples].
 * @param[out] radial      The radial channel.  This is an array whose dimension
 *                         is [nSamples].
 * @param[out] transverse  The transverse channel.  This is an array whose
 *                         dimension is [nSamples].
 * @throws std::invalid_argument if nSamples is positive and north, east,
 *         radial, or transverse is NULL.
 */
template<typename T>
void northEastToRadialTransverse(const int nSamples,
                                 const T backAzimuth,
                                 const T north[],
                                 const T east[],
                                 T *radial[],
                                 T *transverse[]);
/*!
 * @brief Rotates a (radial,transverse) channel pair to (north,east)
 * @param[in] nSamples     The number of samples in the seismograms.
 * @param[in] backAzimuth  The receiver to source azimuth in radians.
 * @param[in] radial       The radial channel.  This is an array whose dimension
 *                         is [nSamples].
 * @param[in] transverse   The transverse channel.  This is an array whose
 *                         dimension is [nSamples].
 * @param[out] north       The north channel.  This is an array whose dimension
 *                         is [nSamples].
 * @param[out] east        The east channel.  This is an array whose dimension
 *                         is [nSamples].
 * @throws std::invalid_argument if nSamples is positive and vertical, north,
 *         east, longitudinal, radial, or transverse is NULL.
 */
template<typename T>
void radialTransverseToNorthEast(const int nSamples,
                                 const T backAzimuth,
                                 const T radial[],
                                 const T transverse[],
                                 T *north[],
                                 T *east[]);

/*!
 * @brief Rotates a (vertical,north,east) channel triplet to 
 *        (longitudinal,radial,transverse).
 * @param[in] nSamples        The number of samples in the seismograms.
 * @param[in] backAzimuth     The receiver to source azimuth in radians.
 * @param[in] incidenceAngle  The angle of incidence to the station in radians.
 *                            This is measured positive from vertical. 
 * @param[in] vertical        The vertical channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[in] north           The north channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[in] east            The east channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[out] longitudinal   The longitudinal channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[out] radial         The radial channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[out] transverse     The transverse channel.  This is an array whose
 *                            dimension is [nSamples].
 * @throws std::invalid_argument if nSamples is positive and north, east,
 *         radial, or transverse is NULL.
 */
template<typename T>
void verticalNorthEastToLongitudinalRadialTransverse(
    const int nSamples,
    const T backAzimuth,
    const T incidenceAngle,
    const T vertical[],
    const T north[],
    const T east[],
    T *longitudinalIn[],
    T *radialIn[],
    T *transverseIn[]
    );

/*!
 * @brief Rotates a (longitudinal,radial,transverse) channel triplet to
 *        (vertical,north,east).
 * @param[in] nSamples        The number of samples in the seismograms.
 * @param[in] backAzimuth     The receiver to source azimuth in radians.
 * @param[in] incidenceAngle  The angle of incidence to the station in radians.
 *                            This is measured positive from vertical.
 * @param[in] longitudinal    The longitudinal channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[in] radial          The radial channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[in] transverse      The transverse channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[out] vertical       The vertical channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[out] north          The north channel.  This is an array whose
 *                            dimension is [nSamples].
 * @param[out] east           The east channel.  This is an array whose
 *                            dimension is [nSamples].
 * @throws std::invalid_argument if nSamples is positive and vertical, north,
 *         east, longitudinal, radial, or transverse is NULL.
 */
template<typename T>
void longitudinalRadialTransverseToVerticalNorthEast(
    const int nSamples,
    const T backAzimuth,
    const T incidenceAngle,
    const T longitudinal[],
    const T radial[],
    const T transverse[],
    T *verticalIn[],
    T *northIn[],
    T *eastIn[]
    );

}
#endif
