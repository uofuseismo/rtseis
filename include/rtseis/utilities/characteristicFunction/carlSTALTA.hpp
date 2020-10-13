#ifndef RTSEIS_UTILITIES_CHARATERISTICFUNCTION_CARLSTALTA_HPP
#define RTSEIS_UTILITIES_CHARATERISTICFUNCTION_CARLSTALTA_HPP
#include <memory>
#include "rtseis/enums.hpp"
namespace RTSeis::Utilities::CharacteristicFunction
{
/*!
 * @brief Implements Carl Johnson's short-term-average /long-term-average
 *        (STA-LTA) characteristic function:
 *        \f[
 *          y[n] = \frac{1}{N_{sta}} \sum_{i=0}^{N_{sta}-1} |x[n-i]-LTA[n-i]|
 *               - R \frac{1}{N_{lta}} \sum_{i=0}^{N_{lta}-1} |x[n-i]-LTA[n-i]|
 *               - |STA[n] - LTA[n]|
 *               - Q
 *        \f] 
 *        where \f$ N_{sta} \f$ is the number of samples in the short term
 *        average window, \f$ N_{lta} \f$ is the number of samples in the
 *        long term average window, STA is the short-term average a point,
 *        LTA is the long-term average, \f$ R \f$ is the `ratio', and
 *        \f$ Q \f$ is the `quiet' parameter.   Effectively, for a fixed
 *        threshold, the sensitivity increases for any values that make
 *        \f$ y[n] \f$ less.  So as \f$ \R \f$ decreases the detector
 *        becomes more sensitive.  Likewise, as \f$ Q \f$ increases the
 *        detector becomes more sensitive.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 * @ingroup  rtseis_utils_characteristicFunction
 * @note This is not precisely what is done in Earthworm.  Carl appears to hold
 *       constant the LTA in the rectified terms (the first two terms) for a 
 *       window of duraction of STA.  For this reason I expect this 
 *       implementation produces a characteristic function with a little more
 *       variance.
 */
template<class T = double,
         RTSeis::ProcessingMode E = RTSeis::ProcessingMode::POST_PROCESSING>
class CarlSTALTA
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    CarlSTALTA();
    /*!
     * @brief Copy constructor.
     * @param[in] stalta  The STA/LTA class from which to initialize this class.
     */
    CarlSTALTA(const CarlSTALTA &stalta);
    /*!
     * @brief Move constructor.
     * @param[in,out] stalta  The STA/LTA class from which to initialize this
     *                        class.  On exit, stalta's behavior is undefined.
     */
    CarlSTALTA(CarlSTALTA &&stalta) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] stalta   The STA/LTA class to copy to this.
     * @result A deep copy of stalta.
     */
    CarlSTALTA& operator=(const CarlSTALTA &stalta);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] stalta  The STA/LTA class whose memory will be moved to
     *                        this.  On exit, stalta's behavior is undefined.
     * @result The memory from stalta moved to this.
     */
    CarlSTALTA& operator=(CarlSTALTA &&stalta) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~CarlSTALTA();
    /*! 
     * @brief Resets the class and releases all memory.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the classic STA/LTA.
     * @param[in] nSTA   The number of samples in the short-term average window.
     *                   This must be at least 2.
     * @param[in] nLTA   The number of samples in the long-term average window.
     *                   This must be greater than nsta.
     * @param[in] ratio  This is the ratio that scales the rectified long-term
     *                   average.  To start try 2.3.
     * @param[in] quiet  This is the shift to subtract from the characteristic
     *                   function.  To start, try 4.
     * @throws std::invalid_argument if nSTA or nLTA is too small.
     */
    void initialize(int nSTA, int nLTA, double ratio, double quiet);
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Gets the initial condition length.
     * @result The initial condition lengths.  Here, result.first corresponds
     *         to the length of the numerator initial conditions and
     *         result.second corresponds to the length of the denominator
     *         initial conditions.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    std::pair<int, int> getInitialConditionLength() const;
    /*!
     * @brief Sets the initial conditions.
     * @param[in] nzSTA   The length of the STA initial conditions.
     *                    This should be \c getInitialConditionLength().first.
     * @param[in] zSTA    The initial conditions for the STA.
     *                    This has dimension [nzSTA].
     * @param[in] zSTAR   The initial conditions for the rectified STA.
     *                    This has dimension [nzSTA].
     * @param[in] nzLTA   The length of the LTA initial conditions. 
     *                    This should be \c getInitialConditionLength().second.
     * @param[in] zLTA    The initial conditions for the LTA.
     *                    This has dimension [nzLTA].
     * @param[in] zLTAR   The initial conditions for the rectified LTA.
     *                    This has dimension [nzLTA].
     * @throws std::runtime_error if the class is not initialized.
     * @throws std::invalid_argument if the array dimension are incorrect
     *         or any arrays are NULL.
     * @sa \c isInitialized(), \c getInitialConditionLength().
     */
    void setInitialConditions(int nzSTA,
                              const double zSTA[],
                              const double zSTAR[],
                              int nzLTA,
                              const double zLTA[],
                              const double zLTAR[]);
    /*!
     * @brief Applies the STA/LTA filter.
     * @param[in] nx   The number of samples in the input signal.
     * @param[in] x    The signal to filter.  This is an array whose 
     *                 dimension is [n].
     * @param[out] y   The classic STA/LTA charactersistic function.
     *                 This is an array whose dimension [n].
     * @throws std::invalid_argument if any of the arrays are NULL.
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
    class CarlSTALTAImpl;
    std::unique_ptr<CarlSTALTAImpl> pImpl;
};
}
#endif
