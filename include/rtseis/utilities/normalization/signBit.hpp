#ifndef RTSEIS_UTILS_NORMALIZATION_SB_HPP
#define RTSEIS_UTILS_NORMALIZATION_SB_HPP 1
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
        /*!
         * @brief Default constructor.
         */
        SignBit(void);
        /*
         * @brief Copy constructor.
         * @param[in] signBit  Sign-bit class from which to initialize.
         */
        SignBit(const SignBit &signBit);
        /*!
         * @brief Copy assignment operator.
         * @param[in] signBit  Sign-bit normalization class to copy.
         * @result A deep copy of the SignBit class.
         */
        SignBit& operator=(const SignBit &signBit);
        /*!
         * @brief Default destructor.
         */
        ~SignBit(void);
         /*!
          * @brief Initializes the class.
          * @result 0 indicates success.
          */
         int initialize(void);
        /*!
         * @brief Determines if the module is initialized.
         * @retval True indicates that the module is initialized.
         * @retval False indicates that the module is not initialized.
         */
        bool isInitialized(void) const;
        /*!
         * @brief Sets the initial conditions for the sign-bit normalization.
         *        Note, the sign operation is applied to each sample
         *        so there are no initial conditions to set.
         * @result 0 indicates success.
         */
        int setInitialConditions(void);
        /*!
         * @brief Applies the sign-bit normalization.
         * @param[in] nx   Number of points in signal.
         * @param[in] x    Signal to sign-bit normalized.  This is an array of
         *                 dimension [nx].
         * @param[out] y   The sign-bit normalized signal.  This has
         *                 dimension [nx].
         * @result 0 indicates success.
         */
        int apply(const int nx, const double x[], double y[]);
        /*! @copydoc apply */
        int apply(const int nx, const float  x[], float  y[]);
        /*!
         * @brief Resets the filter to the initial conditions specified
         *        in setInitialConditions().  Note, the sign operation
         *        is applied each sample so there are no initial
         *        parameters to restore.
         * @result 0 indicates success.
         */
        int resetInitialConditions(void);
        /*!
         * @brief Clears variables in class and restores defaults.
         *        This class will have to be re-initialized to use again.
         */
        void clear(void);
    private:
        class SignBitImpl;
        std::unique_ptr<SignBitImpl> pSignBit_;

}; // End SignBit


}; // End Normalization
}; // End Utilities
}; // End RTSeis
#endif 

