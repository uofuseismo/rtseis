#ifndef RTSEIS_UTILITIES_TRANSFORMS_WAVELETS_IWAVELETS_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_WAVELETS_IWAVELETS_HPP 1
#include <memory>
#include <complex>
namespace RTSeis::Utilities::Transforms::Wavelets
{
/*!
 * @class Wavelet wavelet.hpp "include/rtseis/utilities/transforms/wavelets/wavelet.hpp"
 * @brief An abstract base class defining a wavelet.
 * @note This work is based on "A Practical Guide to Wavelet Analysis"
 *       by Torrence and Compo, 1998.
 * @ingroup rtseis_utils_transforms_wavelets
 * @author Ben Baker (University of Utah) distributed under the MIT license.
 */
class IContinuousWavelet
{
public:
    /*!
     * @brief Destructor.
     */
    virtual ~IContinuousWavelet() = default;
    /*!
     * @result A copy of the the continuous wavelet class.
     */
    virtual std::unique_ptr<IContinuousWavelet> clone() const = 0; 
    /*!
     * @brief Evaluates the wavelet at the scaling factors.
     * @param[in] n   The number of scaling factors.
     * @param[in] s   The scaling factors in seconds/radian.
     */
    //virtual void evaluate(int n, double s[]) const = 0;
    /*! @copydoc evaluate */
    //virtual void evaluate(int n, float s[]) const = 0;
    /*!
     * @brief Sets the sampling period in evaluating the wavelet and cone of
     *        influence.
     * @param[in] dt   The sampling period in seconds.
     * @throws std::invalid_argument if dt is not positive.
     */
    //virtual void setSamplingPeriod(const double dt) = 0;
    /*!
     * @brief Gets the sampling period (seconds).
     */
    //virtual double getSamplingPeriod() const noexcept = 0;
    /*!
     * @result The e-folding time.
     * @param[in] s   The wavelet's scale in seconds/radian.
     * @note Basically what's going on here is that to perform the Fourier
     *       transform we assume that the signal is periodic.  This means,
     *       that the wavelet will `see' data outside of the window.
     *       This manifests as an edge effect.  The e-folding factor is 
     *       a heuristic where the wavelet's amplitude spectrum decays below
     *       some `acceptable' amount.
     */
    //[[nodiscard]] virtual double getEFoldingTime(double s) const = 0;
    /*!
     * @result The wavelength in seconds/radian.
     */
    //[[nodiscard]] virtual double getWavelength(double s) const = 0;
    /*!
     * @brief Evaluates the wavelet at the given scales.
     * @param[in] n   The number of samples at which to evaluate the wavelet.
     * @param[in] s   The dimensionless scale at which to evalute the
     *                wavelet.
     * @param[out] daughter  The daughter wavelet evaluated at the n samples.
     *                       This is an array whose dimension is [n]. 
     * @result The wavelet daughter evaluted at the given scales. 
     * @throws std::runtime_error if the class is not initialized.
     */  
    virtual void evaluate(int n, double scale,
                          std::complex<double> *daughter[]) const = 0;
    virtual void evaluate(int n, float scale,
                          std::complex<float> *daughter[]) const = 0;
    /*!
     * @result The cone of influence scalar for the given nominal wavenumber.
     *         The cone of influence at a time is given by scalar*time
     *         where time is measured from window's start.
     */
    //virtual double computeConeOfInfluenceScalar() const noexcept = 0;
};
}
#endif
