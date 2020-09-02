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
     * @brief Sets the nominal wavenumber.
     * @param[in] k0   The wavenumber.  The default is 6.
     * @throws std::invalid_argument if this is not positive.
     */
    void setWaveNumber(double k0);
    /*!
     * @result The nominal wavenumber.
     */
    double getWaveNumber() const noexcept;

    /*!
     * @brief Sets the sampling period used in evaluating the wavelet and 
     *        cone of influence.
     * @param[in] dt   The sampling period in seconds.
     * @throws std::invalid_argument if dt is not positive.
     */
    void setSamplingPeriod(double dt) override;
    /*!
     * @result The sampling period in seconds.
     */
    double getSamplingPeriod() const noexcept override;
    /*! @} */

    /*!
     * @result The e-folding time for the derivative of a Gaussian wavelet.
     *         This has units of seconds/radian.
     */
    //[[nodiscard]] double getEFoldingTime(double s) const override;
    //[[nodiscard]] double getWavelength(double s) const override;
    void evaluate(int n, double scale,
                  const double k[], std::complex<double> *w[]) const override;
    void evaluate(int n, float scale,
                  const float k[],
                  std::complex<float> *daughter[]) const override;
    /*!
     * @result The cone of influence scalar for the given nominal wavenumber.
     *         The cone of influence at a time is given by scalar*time 
     *         where time is measured from window's start.
     */
    [[nodiscard]] double computeConeOfInfluenceScalar() const noexcept override;
    /*! @} */
private:
    class MorletImpl;
    std::unique_ptr<MorletImpl> pImpl;
};
}
#endif
