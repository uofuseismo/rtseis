#ifndef RTSEIS_UTILITIES_TRANSFORMS_WAVELETS_MORLET_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_WAVELETS_MORLET_HPP 1
#include <memory>
#include "rtseis/utilities/transforms/wavelets/iwavelets.hpp"
namespace RTSeis::Utilities::Transforms::Wavelets
{
/*!
 * @class Morlet morlet.hpp "include/rtseis/utilities/transforms/wavelets/morlet.hpp"
 * @brief Computes the Morlet wavelet.
 * @ingroup rtseis_utils_transforms_wavelets
 * @author Ben Baker (University of Utah) distributed under the MIT license.
 */
class Morlet  : public IContinuousWavelet
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    Morlet();
    /*!
     * @brief Copy constructor.
     * @param[in] morlet  The Morlet class from which to initialize this class.
     */ 
    Morlet(const Morlet &morlet);
    /*!
     * @brief Move constructor.
     * @param[in,out] morlet  The Morlet class from from which to initialize
     *                        this class.  On exit, morlet's behavior is
     *                        undefined.
     */
    Morlet(Morlet &&morlet) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] morlet   The class to copy to this.
     * @result A deep copy of morlet.
     */
    Morlet& operator=(const Morlet &morlet);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] morlet  The class whose memory will be moved to this.
     *                        On exit, morlet's behavior is undefined. 
     * @result The memory from morlet moved to this.
     */
    Morlet& operator=(Morlet &&morlet) noexcept;
    /*!
     * @result A deep copy of this class. 
     */
    std::unique_ptr<IContinuousWavelet> clone() const override;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    virtual ~Morlet();
    /*!
     * @brief Resets the class to the defaults.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Parameters
     * @{
     */
    /*!
     * @brief Sets the wavelet parameter. 
     * @param[in] omega0   The wavelet parameter.  This represents a tradeoff
     *                     between time and frequency resolution.  
     *                     The default is 6.
     * @throws std::invalid_argument if this is not positive.
     */
    void setParameter(double omega0);
    /*!
     * @result The nominal wavenumber.
     */
    double getParameter() const noexcept;
    /*! @} */

    /*!
     * @result The e-folding time for the derivative of a Gaussian wavelet.
     *         This has units of seconds/radian.
     */
    //[[nodiscard]] double getEFoldingTime(double s) const override;
    //[[nodiscard]] double getWavelength(double s) const override;
    /*!
     * @brief Evalutes the Morlet wavelet
     *        \f [
     *           W = \pi^{-1/4}
     *               e^{i \frac{\omega_0 x}{2 s} }
     *               e^{-\frac{1}{2} \left ( \frac{x}{s} \right )^2}
     *        \f ]  
     *        where \f$ x \f$ is the sample and \f$ \omega_0 \f$ is the 
     *        center frequency (nominally 6).
     * @param[in] n       The number of samples at which to evalute the wavelet.
     * @param[in] scale   Nominally, this is the unitless scale.
     *                    However, you may find it more convenient to note that
     *                    \f$ s = \frac{F_c}{2 \pi f} \f$ where s is the scale,
     *                    \f$ F_c \f$ is the center frequency of the wavelet
     *                    (nominally 6), and \f$ f \f$ is the frequency of
     *                    of interest in Hz.  In this case, you would define the
     *                    scale as \f$ s = \frac{F_c f_s}{2 \pi f} \f$
     *                    where \f$ f_s \f$ is the sampling rate in Hz.
     * @param[out] daughter  The daughter wavelet evaluated at the n samples.
     *                       This is an array whose dimension is [n].
     *                       This will be centered at \f$ x = 0 \f$ which 
     *                       corresponds to index n/2 for n odd.
     */ 
    void evaluate(int n, double scale,
                  std::complex<double> *w[]) const override;
    /*! @copydoc evaluate() */
    void evaluate(int n, float scale,
                  std::complex<float> *daughter[]) const override;
    /*!
     * @result The cone of influence scalar for the given nominal wavenumber.
     *         The cone of influence at a time is given by scalar*time 
     *         where time is measured from window's start.
     */
    //[[nodiscard]] double computeConeOfInfluenceScalar() const noexcept override;
    /*! @} */
private:
    class MorletImpl;
    std::unique_ptr<MorletImpl> pImpl;
};
}
#endif
