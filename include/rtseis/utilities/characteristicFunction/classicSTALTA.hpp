#ifndef RTSEIS_UTILITIES_CHARATERISTICFUNCTION_CLASSICSTALTA_HPP
#define RTSEIS_UTILITIES_CHARATERISTICFUNCTION_CLASSICSTALTA_HPP
#include <memory>
namespace RTSeis::Utilities::CharacteristicFunction
{
template<class T> class ClassicSTALTAImpl;
namespace PostProcessing
{
/*!
 * @brief Implements the post-processing variant of the classic
 *        short-term-average to long-term-average
 *        (STA/LTA):
 *        \f[
 *          y[n] = \frac{\sum_{i=0}{N_{sta}-1} x[n-i]^2  }
 *                      {\sum_{i=0}{N_{lta}-1} x[n-i]^2 }
 *        \f] 
 *       where \f$ N_{sta} \f$ is the number of samples in the short term
 *       average window and \f$ N_{lta} \f$ is the number of samples in the
 *       long term average window.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 * @ingroup  rtseis_utils_characteristicFunction
 */
template<class T = double>
class ClassicSTALTA
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    ClassicSTALTA();
    /*!
     * @brief Copy constructor.
     * @param[in] stalta  The classic STA/LTA class from which to initialize
     *                    this class.
     */
    ClassicSTALTA(const ClassicSTALTA &stalta);
    /*!
     * @brief Move constructor.
     * @param[in,out] stalta  The classic STA/LTA class from which to initialize
     *                        this class.  On exit, stalta's behavior is
     *                        undefined.
     */
    ClassicSTALTA(ClassicSTALTA &&stalta) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] stalta   The classic STA/LTA class to copy to this.
     * @result A deep copy of stalta.
     */
    ClassicSTALTA& operator=(const ClassicSTALTA &stalta);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] stalta  The classic STA/LTA class whose memory will be
     *                        moved to this.  On exit, stalta's behavior is
     *                        undefined.
     * @result The memory from stalta moved to this.
     */
    ClassicSTALTA& operator=(ClassicSTALTA &&stalta) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~ClassicSTALTA();
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
     * @throws std::invalid_argument 
     */
    void initialize(int nSTA, int nLTA);
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
     * @param[in] nzNum   The length of the numerator initial conditions
     *                    array.
     * @param[in] zNum    The initial conditions in the numerator.  This is
     *                    an array whose dimension is [nzNum].
     * @param[in] nzDen   The length of the denominator initial conditions
     *                    array.
     * @param[in] zDen    The initial conditions in the denominator.  This is
     *                    an array whose dimension is [nzDen].
     * @throws std::invalid_argument if the array lengths are inconsistent with
     *         \c getInitialConditionLength() or either of the arrays is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized(), \c getInitialConditionLength()
     */
    void setInitialConditions(int nzNum, const double zNum[],
                              int nzDen, const double zDen[]);
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
private:
    std::unique_ptr<ClassicSTALTAImpl<T>> pImpl;
};
}

namespace RealTime
{
/*!
 * @brief Implements the post-processing variant of the classic
 *        short-term-average to long-term-average
 *        (STA/LTA):
 *        \f[
 *          y[n] = \frac{\sum_{i=0}{N_{sta}-1} x[n-i]^2  }
 *                      {\sum_{i=0}{N_{lta}-1} x[n-i]^2 }
 *        \f] 
 *       where \f$ N_{sta} \f$ is the number of samples in the short term
 *       average window and \f$ N_{lta} \f$ is the number of samples in the
 *       long term average window.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 * @ingroup  rtseis_utils_characteristicFunction
 */
template<class T = double>
class ClassicSTALTA
{
public:
    /*! @name Constructors
     * @{
     */
    /*! 
     * @brief Default constructor.
     */
    ClassicSTALTA();
    /*!
     * @brief Copy constructor.
     * @param[in] stalta  The classic STA/LTA class from which to initialize
     *                    this class.
     */
    ClassicSTALTA(const ClassicSTALTA &stalta);
    /*!
     * @brief Move constructor.
     * @param[in,out] stalta  The classic STA/LTA class from which to initialize
     *                        this class.  On exit, stalta's behavior is
     *                        undefined.
     */
    ClassicSTALTA(ClassicSTALTA &&stalta) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*! 
     * @brief Copy assignment operator.
     * @param[in] stalta   The classic STA/LTA class to copy to this.
     * @result A deep copy of stalta.
     */
    ClassicSTALTA& operator=(const ClassicSTALTA &stalta);
    /*! 
     * @brief Move assignment operator.
     * @param[in,out] stalta  The classic STA/LTA class whose memory will be
     *                        moved to this.  On exit, stalta's behavior is
     *                        undefined.
     * @result The memory from stalta moved to this.
     */
    ClassicSTALTA& operator=(ClassicSTALTA &&stalta) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~ClassicSTALTA();
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
     * @throws std::invalid_argument 
     */
    void initialize(int nSTA, int nLTA);
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
     * @param[in] nzNum   The length of the numerator initial conditions
     *                    array.
     * @param[in] zNum    The initial conditions in the numerator.  This is
     *                    an array whose dimension is [nzNum].
     * @param[in] nzDen   The length of the denominator initial conditions
     *                    array.
     * @param[in] zDen    The initial conditions in the denominator.  This is
     *                    an array whose dimension is [nzDen].
     * @throws std::invalid_argument if the array lengths are inconsistent with
     *         \c getInitialConditionLength() or either of the arrays is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized(), \c getInitialConditionLength()
     */
    void setInitialConditions(int nzNum, const double zNum[],
                              int nzDen, const double zDen[]);
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
    std::unique_ptr<ClassicSTALTAImpl<T>> pImpl;
};
}
}
#endif
