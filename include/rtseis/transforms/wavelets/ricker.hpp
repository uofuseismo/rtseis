#ifndef RTSEIS_UTILITIES_TRANSFORMS_WAVELETS_RICKER_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_WAVELETS_RICKER_HPP 1
#include <memory>
#include "rtseis/utilities/transforms/wavelets/iwavelets.hpp"
namespace RTSeis::Utilities::Transforms::Wavelets
{
/*!
 * @class Ricker ricker.hpp "include/rtseis/utilities/transforms/wavelets/ricker.hpp"
 * @brief Computes the Ricker wavelet (this is the second derivative of a
 *        Gaussian wavelet).
 * @ingroup rtseis_utils_transforms_wavelets
 * @author Ben Baker (University of Utah) distributed under the MIT license.
 */
class Ricker : public IContinuousWavelet
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    Ricker();
    /*!
     * @brief Copy constructor.
     * @param[in] ricker  The derivative of a Gaussian class from which to
     *                    initialize this class.
     */ 
    Ricker(const Ricker &ricker);
    /*!
     * @brief Move constructor.
     * @param[in,out] ricker  The derivative of a Gaussian class from which to
     *                        initialize this class.  On exit, ricker's behavior
     *                        is undefined.
     */
    Ricker(Ricker &&ricker) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] ricker   The class to copy to this.
     * @result A deep copy of ricker.
     */
    Ricker& operator=(const Ricker &ricker);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] ricker   The class whose memory will be moved to this.
     *                         On exit, ricker's behavior is undefined.
     * @result The memory from ricker moved to this.
     */
    Ricker& operator=(Ricker &&ricker) noexcept;
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
    virtual ~Ricker();
    /*!
     * @brief Resets the class to the defaults.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Parameters
     * @{
     */
    /*!
     * @brief Sets the sampling period in evaluating the wavelet and cone of
     *        influence.
     * @param[in] dt   The sampling period in seconds.
     * @throws std::invalid_argument if dt is not positive.
     */
    void setSamplingPeriod(double dt) override;
    /*!
     * @result The sampling period in seconds.
     */
    double getSamplingPeriod() const noexcept override;
    /*! @} */

    /*! @name Evaluation
     * @{
     */
    /*!
     * @result The e-folding time for the derivative of a Gaussian wavelet.
     *         This has units of seconds/radian.
     */
    //[[nodiscard]] double getEFoldingTime(double s) const override;
    //[[nodiscard]] double getWavelength(double s) const override;
    void evaluate(int n, double scale,
                  const double k[],
                  std::complex<double> *daughter[]) const override;
    void evaluate(int n, float scale,
                  const float k[],
                  std::complex<float> *daughter[]) const override;
    /*!
     * @result The cone of influence scalar for the second derivative of
     *         a Gaussian.  The cone of influence at a time is given by
     *         scalar*time where time is measured from window's start.
     */
    [[nodiscard]] double computeConeOfInfluenceScalar() const noexcept override;
    /*! @} */
private:
    class RickerImpl;
    std::unique_ptr<RickerImpl> pImpl;
};
}
#endif
