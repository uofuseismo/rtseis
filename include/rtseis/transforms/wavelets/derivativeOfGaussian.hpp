#ifndef RTSEIS_UTILITIES_TRANSFORMS_WAVELETS_DOG_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_WAVELETS_DOG_HPP 1
#include <memory>
#include "rtseis/utilities/transforms/wavelets/iwavelets.hpp"
namespace RTSeis::Utilities::Transforms::Wavelets
{
/*!
 * @class DerivativeOfGaussian derivativeOfGaussian.hpp "include/rtseis/utilities/transforms/wavelets/derivativeOfGaussian.hpp"
 * @brief Computes the derivative of a Gaussian wavelet.
 * @ingroup rtseis_utils_transforms_wavelets
 * @author Ben Baker (University of Utah) distributed under the MIT license.
 */
class DerivativeOfGaussian : public IContinuousWavelet
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    DerivativeOfGaussian();
    /*!
     * @brief Copy constructor.
     * @param[in] dog  The derivative of a Gaussian class from which to
     *                 initialize this class.
     */ 
    DerivativeOfGaussian(const DerivativeOfGaussian &dog);
    /*!
     * @brief Move constructor.
     * @param[in,out] dog  The derivative of a Gaussian class from which to
     *                     initialize this class.  On exit, dog's behavior
     *                     is undefined.
     */
    DerivativeOfGaussian(DerivativeOfGaussian &&dog) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] dog   The class to copy to this.
     * @result A deep copy of dog.
     */
    DerivativeOfGaussian& operator=(const DerivativeOfGaussian &dog);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] dog   The class whose memory will be moved to this.
     *                      On exit, dog's behavior is undefined.  Bad dog.
     * @result The memory from dog moved to this.
     */
    DerivativeOfGaussian& operator=(DerivativeOfGaussian &&dog) noexcept;
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
    virtual ~DerivativeOfGaussian();
    /*!
     * @brief Resets the class to the defaults.
     */
    void clear() noexcept;
    /*! @} */

    /*! @name Parameters
     * @{
     */
    /*!
     * @brief Initializes the derivative of the Gaussian wavelet.
     * @param[in] order   The order of the derivative.
     *                    For example, 2 is a Ricker wavelet. 
     * @param[in] dt      The sampling period in seconds.
     * @throws std::invalid_argument if order is negative or dt is not positive.
     */
    void setOrder(int order);
    /*!
     * @result The derivative order.
     */
    int getOrder() const noexcept;

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
                  const double k[],
                  std::complex<double> *daughter[]) const override;
    void evaluate(int n, float scale,
                  const float k[],
                  std::complex<float> *daughter[]) const override;
    /*!
     * @result The cone of influence scalar for the given derivative order.
     *         The cone of influence at a time is given by scalar*time
     *         where time is measured from window's start.
     */
    [[nodiscard]] double computeConeOfInfluenceScalar() const noexcept override;
    /*! @} */
private:
    class DerivativeOfGaussianImpl;
    std::unique_ptr<DerivativeOfGaussianImpl> pImpl;
};
}
#endif
