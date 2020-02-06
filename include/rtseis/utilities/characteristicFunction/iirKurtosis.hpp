#ifndef RTSEIS_UTILITIES_CHARACTERISTICFUNCTION_IIRKURTOSIS_HPP
#define RTSEIS_UTILITIES_CHARACTERISTICFUNCTION_IIRKURTOSIS_HPP
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis::Utilities::CharacteristicFunction
{
template<class T> class IIRKurtosisImpl;
namespace PostProcessing
{
/*!
 * @brief This is the IIR post-processing implementation for computing the
 *        kurtosis of a signal.  This can be a useful characteristic function
 *        when creating a picker.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 * @note For more details see: Testing the normality of gravitational wave
 *       data with a low cost recursive estimate of the kurtosis - 
 *       E Chassande-Mottin. 
 */ 
template<class T = double>
class IIRKurtosis
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    IIRKurtosis();
    /*!
     * @brief Copy constructor.
     * @param[in] iirKurtosis   The IIR kurtosis class from which to 
     *                          initialize this class.
     */
    IIRKurtosis(const IIRKurtosis &iirKurtosis);
    /*!
     * @brief Move constructor.
     * @param[in,out] iirKurtosis  The IIR kurtosis class from which to 
     *                             initialize this class.  On exit,
     *                             iirKurtosis's behavior is undefined.
     */
    IIRKurtosis(IIRKurtosis &&iirKurtosis) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] iirKurtosis  The IIR kurtosis filter to copy to this.
     * @result A deep copy of the IIR kurtosis filter.
     */
    IIRKurtosis& operator=(const IIRKurtosis &iirKurtosis);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] iirKurtosis  The IIR kurtosis filter whose memory will
     *                             be moved to this.  On exit, iirKurtosis's
     *                             behavior is undefined.
     * @result The memory from iirKurtosis moved to this.
     */
    IIRKurtosis& operator=(IIRKurtosis &&iirKurtosis) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~IIRKurtosis();
    /*!
     * @brief Resets the class and releases memory.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the IIR kurtosis filter. 
     * @param[in] c1  This is the pole in the IIR moving average.  For 
     *                stability it is required that \f$ |c_1| < 1 \f$.
     *                One strategy is to set 
     *                \f$ c_1 = \frac{\Delta T}{W} \f$ where \f$ \Delta T \f$
     *                is the sampling period and \f$ W \f$ is the window length.
     * @throws std::invalid_argument if \f$ |c_1| \ge 1 \f$.
     */
    void initialize(T c1);
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Sets the initial conditions.
     * @param[in] mu1    The mean of the first order moment.
     * @param[in] mu2    The mean of the second order moment.
     * @param[in] k4bar  The initial unbiases kurtosis value.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void setInitialConditions(T mu1, T mu2, T k4bar);
    /*!
     * @brief Applies the kurtosis filter.
     * @param[in] npts   The number of samples in the signal.
     * @param[in] x      The signal whose kurtosis will be computed at each
     *                   sample.  This is an array whose dimension is [npts].
     * @param[out] y     The kurtosis of x.  This is an array whose dimension
     *                   is [npts].
     * @throws std::invalid_argument if x or y is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void apply(int nx, const T x[], T *y[]);
private:
    std::unique_ptr<IIRKurtosisImpl<T>> pImpl;
};
} // End namespace on postprocessing

namespace RealTime
{
/*!
 * @brief This is the IIR real-time implementation for computing the kurtosis
 *        of a signal.  This can be a useful characteristic function when
 *        when creating a picker.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 * @note For more details see: Testing the normality of gravitational wave
 *       data with a low cost recursive estimate of the kurtosis - 
 *       E Chassande-Mottin. 
 */ 
template<class T = double>
class IIRKurtosis
{
public:
    /*! @name Constructors
     * @{
     */
    /*! 
     * @brief Default constructor.
     */
    IIRKurtosis();
    /*! 
     * @brief Copy constructor.
     * @param[in] iirKurtosis   The IIR kurtosis class from which to 
     *                          initialize this class.
     */
    IIRKurtosis(const IIRKurtosis &iirKurtosis);
    /*! 
     * @brief Move constructor.
     * @param[in,out] iirKurtosis  The IIR kurtosis class from which to 
     *                             initialize this class.  On exit,
     *                             iirKurtosis's behavior is undefined.
     */
    IIRKurtosis(IIRKurtosis &&iirKurtosis) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] iirKurtosis  The IIR kurtosis filter to copy to this.
     * @result A deep copy of the IIR kurtosis filter.
     */
    IIRKurtosis& operator=(const IIRKurtosis &iirKurtosis);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] iirKurtosis  The IIR kurtosis filter whose memory will
     *                             be moved to this.  On exit, iirKurtosis's
     *                             behavior is undefined.
     * @result The memory from iirKurtosis moved to this.
     */
    IIRKurtosis& operator=(IIRKurtosis &&iirKurtosis) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*! 
     * @brief Default destructor.
     */
    ~IIRKurtosis();
    /*! 
     * @brief Resets the class and releases memory.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the IIR kurtosis filter. 
     * @param[in] c1  This is the pole in the IIR moving average.  For 
     *                stability it is required that \f$ |c_1| < 1 \f$.
     *                One strategy is to set 
     *                \f$ c_1 = \frac{\Delta T}{W} \f$ where \f$ \Delta T \f$
     *                is the sampling period and \f$ W \f$ is the window length.
     * @throws std::invalid_argument if \f$ |c_1| \ge 1 \f$.
     */
    void initialize(T c1);
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Sets the initial conditions.
     * @param[in] mu1    The mean of the first order moment.
     * @param[in] mu2    The mean of the second order moment.
     * @param[in] k4bar  The initial unbiases kurtosis value.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void setInitialConditions(T mu1, T mu2, T k4bar);
    /*!
     * @brief Applies the kurtosis filter.
     * @param[in] npts   The number of samples in the signal.
     * @param[in] x      The signal whose kurtosis will be computed at each
     *                   sample.  This is an array whose dimension is [npts].
     * @param[out] y     The kurtosis of x.  This is an array whose dimension
     *                   is [npts].
     * @throws std::invalid_argument if x or y is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void apply(int nx, const T x[], T *y[]);
    /*!
     * @brief Resets the filter's default initial conditions or the initial
     *        conditions set in \c setInitialConditions().  This is useful when
     *        dealing with a gap.
     * @throws std::runtime_error if the class is not set.
     * @sa \c isInitialized(), \c setInitialConditions()
     */
    void resetInitialConditions();
private:
    std::unique_ptr<IIRKurtosisImpl<T>> pImpl;
};
} // End namespace on RealTime
}
#endif
