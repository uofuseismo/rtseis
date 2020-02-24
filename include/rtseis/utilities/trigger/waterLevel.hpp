#ifndef RTSEIS_UTILITIES_TRIGGER_WATERLEVEL_HPP
#define RTSEIS_UTILITIES_TRIGGER_WATERLEVEL_HPP
#include <memory>
#include <vector>

namespace RTSeis::Utilities::Trigger::PostProcessing
{
/*!
 * @brief Defines the waterlevel-based trigger.  Effectively, this begins
 *        a trigger window when the characteristic function exceeds some
 *        waterlevel (tolerance) and finishes the trigger window when the
 *        characteristic function drops below some waterlevel (tolerance).
 */ 
template<class T = double>
class WaterLevel
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    WaterLevel();
    /*!
     * @brief Copy constructor.
     * @param[in] trigger  The trigger class from which to initialize
     *                     this class.
     */
    WaterLevel(const WaterLevel &trigger);
    /*!
     * @brief Move constructor.
     * @param[in,out] trigger  The waterlevel trigger class from which to
     *                         intialize this class.  On exit, trigger's
     *                         behavior is undefined.
     */ 
    WaterLevel(WaterLevel &&trigger) noexcept;
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] trigger   The waterlevel trigger class to copy to this.
     * @result A deep copy of trigger.
     */
    WaterLevel& operator=(const WaterLevel &trigger);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] trigger  The waterlevel trigger class whose memory will
     *                         be moved to this.  On exit, trigger's behavior
     *                         is undefined.
     * @result The memory from trigger moved to this.
     */
    WaterLevel& operator=(WaterLevel &&trigger) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor.
     */
    ~WaterLevel();
    /*!
     * @brief Releases memory and resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the trigger class.
     * @param[in] onTolerance   When the characteristic function first exceeds
     *                          this tolerance the trigger window commences.
     * @param[in] offTolerance  When the characteristic function first drops
     *                          below this tolerance the trigger window
     *                          finalizes.
     */
    void initialize(double onTolerance, double offTolerance);

    /*!
     * @brief Applies the triggering algorithm to the data.
     * @param[in] nSamples  The number of samples in the signal.
     * @param[in] x         The characteristic function from which to compute
     *                      the triggers.  This is an array whose dimension is
     *                      [nSamples].
     * @throws std::invalid_argument if nSamples is positive and x is NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void apply(const int npts, const T x[]);

    /*!
     * @brief Determines the number of trigger windows.
     * @result The number of trigger windows. 
     */
    int getNumberOfWindows() const noexcept;
    /*!
     * @brief Gets the trigger windows.
     * @param[in] nWindows  The length of the windows array.  This must equal
     *                      \c getNumberOfTriggerWindows().
     * @param[out] windows  The trigger windows.  Here, windows[i].first
     *                      defines the starting sample of the trigger and
     *                      windows[i].second defines the end sample of the
     *                      trigger window.  This is an array whose dimension
     *                      is [nWindows].
     * @throws std::invalid_argument if nWindows is invalid or nWindows is
     *         positive and windows is NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @sa \c isInitialized()
     */
    void getWindows(int nWindows, std::pair<int, int> *windows[]) const;
    std::vector<std::pair<int, int>> getWindows() const;
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept; 

private:
    class WaterLevelImpl;
    std::unique_ptr<WaterLevelImpl> pImpl;
};

}
#endif
