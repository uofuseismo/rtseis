#ifndef RTSEIS_UTILITIES_TRANSFORMS_DFTUTILITIES_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_DFTUTILITIES_HPP 1
#include <cmath>
#include <complex>
#include <vector>

/*!
 * @defgroup rtseis_utils_transforms_utils Utility Functions
 * @brief Utility routines accompanying the DFT.
 * @ingroup rtseis_utils_transforms
 */
namespace RTSeis::Utilities::Transforms::DFTUtilities
{
/*! @name Unwrap
 * @{
 */
/*! 
 * @brief Unwraps the phase, p, by changing the absolute jumps
 *        greater than \f$ \pi \f$ to their \f$ 2 \pi \f$ 
 *        complement.  If tolerance is set, then the jump is
 *        \f$ 2 \cdot tol \f$.
 * @param[in] p    The array of phase angles in radians to unwrap.
 * @param[in] tol  The jump tolerance specified in radians which
 *                 must be positive.  The default is \f$ \pi \f$.
 * @result The unwrapped phase angles in radians.  This will have
 *         dimension [p.size()].
 * @throws std::invalid_argument if the arguments are wrong.
 * @ingroup rtseis_utils_transforms_utils
 */
std::vector<double> unwrap(const std::vector<double> &z,
                           const double tol = M_PI);
/*! 
 * @brief Unwraps the phase, p, by changing the absolute jumps
 *        greater than \f$ \pi \f$ to their \f$ 2 \pi \f$ 
 *        complement.  If tolerance is set, then the jump is
 *        \f$ 2 \cdot tol \f$.
 * @param[in] n    The number of points in the phase signals.
 * @param[in] p    The array of phase angles in radians to unwrap.
 *                 This has dimension [n].
 * @param[in] q    The unwrapped phase angles in radians.  This
 *                 has dimension [n].
 * @param[in] tol  The jump tolerance specified in radians which
 *                 must be positive.  The default is \f$ \pi \f$.
 * @throws std::invalid_argument if the arguments are incorrect.
 * @ingroup rtseis_utils_transforms_utils
 */
void unwrap(const int n, const double p[], double *q[], 
            const double tol = M_PI);
/*! @} */

/*! @name Phase
 * @{
 */
/*!
 * @brief Computes the phase angle, i.e., the angle between the
 *        imaginary and real parts of z 
 *        \f[ \phi
 *          = \arctan \left ( \frac{ \Im \{z \} }
 *                                 { \Re \{z \} } \right ) \f].
 * @param[in] z     The array of complex numbers of which to compute the angle.
 * @param[in] lwantDeg  If true then phi is given in degrees.  
 *                      If false then phi is given in radians.
 *                      This is the default.
 * @result The phase angle between the real and imaginary parts of z.
 *         This will have dimension [z.size()].  By default the phase angle
 *         is in radians though, if lwantDeg is true, the phase angle
 *         would then be given in degrees.
 * @throws std::invalid_argument if the arguments are incorrect.
 * @ingroup rtseis_utils_transforms_utils
 */
std::vector<double> phase(const std::vector<std::complex<double>> &z,
                          const bool lwantDeg = false);
/*!
 * @brief Computes the phase angle, i.e., the angle between the
 *        imaginary and real parts of z 
 *        \f[ \phi
 *          = \arctan \left ( \frac{ \Im \{z \} }
 *                                 { \Re \{z \} } \right ) \f].
 * @param[in] n     Number of points in array z.
 * @param[in] z     The array of complex numbers of which to compute
 *                  the angle.  This has dimension [n].
 * @param[out] phi  The phase angle between the real and imaginary 
 *                  parts of z.  This has dimension [n].  By
 *                  default this will be in radians but it can
 *                  optionally be given in degrees.
 * @param[in] lwantDeg  If true then phi is given in degrees.  
 *                      If false then phi is given in radians.
 *                      This is the default.
 * @throws std::invalid_argument if hte arguments are incorrect.
 * @ingroup rtseis_utils_transforms_utils
 */
void phase(const int n, const std::complex<double> z[], double *phi[],
           const bool lwantDeg = false);
/*! @} */

/*! @name Real to Complex DFT Frequencies
 * @{
 */
/*!
 * @brief Computes the DFT sampling frequencies corresponding to the DFT of a 
 *        real time domain signal.
 * @param[in] nSamples        The number of time domain samples.  This must be
 *                            positive.
 * @param[in] samplingPeriod  The sampling period in seconds.  This must
 *                            be positive. 
 * @result The frequencies (Hz) at which the DFT was computed.
 * @throws std::invalid_argument if nSamples or samplingPeriod is not positive.
 * @sa DFTRealToComplex
 */
std::vector<double>
realToComplexDFTFrequencies(const int nSamples, const double samplingPeriod);
/*!
 * @brief Computes the DFT sampling frequencies corresponding to the DFT of a
 *        a real time domain signal.
 * @param[in] nSamples        The number of time domain samples.  This must be
 *                            positive.
 * @param[in] samplingPeriod  The sampling period in seconds.  This must be
 *                            positive.
 * @param[in] lengthFreqs     The length of the frequencies.  This must be
 *                            at least nSamples/2 + 1. 
 * @param[out] freqs          The frequencies (Hz) at which the real to complex
 *                            DFT was computed.  This is an array of dimension
 *                            [lengthFrequencies] however only the first
 *                            nSamples/2 + 1 samples are defined..
 * @throws std::invalid_argument if nSamples or samplingPeriod is not positive,
 *         lengthFreqs is too small, or freqs is NULL.
 * @sa DFTRealToComplex
 */
void realToComplexDFTFrequencies(const int nSamples,
                                 const double samplingPeriod,
                                 const int lengthFreqs,
                                 double *freqs[]);
/*! @} */


/*!
 * @brief Finds the next number, n2, such that n2 is a power of 2 and
 *        n2 is greater than or equal to n.
 * @param[in] n  Non-negative number of which to find the next power
 *               of 2.
 * @result On successful exit this is a number that is a power of 2
 *         and is greater than or equal to n.
 * @throws std::invalid_argument if n is negative.
 * @throws std::runtime_error if n is too large and there is an overflow.
 * @ingroup rtseis_utils_transforms_utils
 */
int nextPowerOfTwo(const int n);

/*! @name Shuffle
 * @{
 */
/*! 
 * @brief Shifts the zero-frequency component of a transform to the
 *        center of the spectrum.
 * @param[in] x   The transform to shuffle.
 * @result The shuffled transform.  This has the same size as x but now
 *         the zero frequency is at the center of the spectrum.
 */
template<class T> std::vector<T> fftShift(const std::vector<T> &x);
/*!
 * @brief Shifts the zero-frequency component of a transform to the center
 *        of the spectrum.
 * @param[in] nSamples  The number of samples in x.
 * @param[in] x         The transform to shuffle.  This is an array whose
 *                      dimension is [nSamples].
 * @param[out] y        This is the (potentially symmetric) transform where
 *                      the zero-frequency of the transform is now the center
 *                      element.  This is an array whose dimension is [n]. 
 * @throws std::invalid_argument if nSamples is not positive or x or y is NULL.
 */
template<typename T>
void fftShift(const int nSamples, const T x[], T *y[]);
/*! } */

}
#endif
