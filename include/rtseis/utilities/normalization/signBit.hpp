#ifndef RTSEIS_UTILITIES_NORMALIZATION_SIGNBIT_HPP
#define RTSEIS_UTILITIES_NORMALIZATION_SIGNBIT_HPP 1
#include <memory>

namespace RTSeis
{
namespace Utilities
{
/*!
 * @defgroup rtseis_utils_normalization Normalization
 * @brief These are the core real-time and post-processing 
 *        utilities for normalizing a signal.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils
 */
namespace Normalization
{
/*!
 * @class SignBit signBit.hpp "include/rtseis/utilities/normalization/signBit.hpp"
 * @brief This is the core implementation for sign-bit signal normalization.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_normalization
 */
class SignBit
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    SignBit();
    /*
     * @brief Copy constructor.
     * @param[in] signBit  Sign-bit class from which to initialize.
     */
    SignBit(const SignBit &signBit);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] signBit  Sign-bit normalization class to copy.
     * @result A deep copy of the SignBit class.
     */
    SignBit& operator=(const SignBit &signBit);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~SignBit();
    /*! @} */

    /*!
     * @brief Initializes the class.
     * @result 0 indicates success.
     */
    void initialize() noexcept;
    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Sets the initial conditions for the sign-bit normalization.
     *        Note, the sign operation is applied to each sample
     *        so there are no initial conditions to set.
     * @throws std::runtime_error if the class is not inititalized.
     */
    void setInitialConditions();
    /*!
     * @brief Applies the sign-bit normalization.
     * @param[in] nx   Number of points in signal.
     * @param[in] x    Signal to sign-bit normalized.  This is an array of
     *                 dimension [nx].
     * @param[out] y   The sign-bit normalized signal.  This has
     *                 dimension [nx].
     * @throws std::invalid_argument if x or y is NULL and nx is positive.
     * @throws std::runtime_error if the class is not initialized.
     */
    void apply(const int nx, const double x[], double y[]);
    /*! @copydoc apply */
    void apply(const int nx, const float  x[], float  y[]);
    /*!
     * @brief Resets the filter to the initial conditions specified
     *        in setInitialConditions().  Note, the sign operation
     *        is applied each sample so there are no initial
     *        parameters to restore.
     * @throws std::runtime_error if the class is not initialized.
     */
    void resetInitialConditions();
    /*!
     * @brief Clears variables in class and restores defaults.
     *        This class will have to be re-initialized to use again.
     */
    void clear() noexcept;
private:
    class SignBitImpl;
    std::unique_ptr<SignBitImpl> pSignBit_;
}; // End SignBit


}; // End Normalization
}; // End Utilities
}; // End RTSeis
#endif 

