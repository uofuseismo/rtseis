#ifndef RTSEIS_UTILITIES_TRANSFORMS_ENVELOPEANALYTIC_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_ENVELOPEANALYTIC_HPP 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{
namespace Transforms
{

/*!
 * @class EnvelopeAnalytic envelopeAnalytic.hpp "include/rtseis/utilities/transforms/envelopeAnalytic.hpp"
 * @brief Computes the envelope of a signal using.  This is for post-processing
 *        only and works by computing the Hilbert transform to obtain
 *        the analytic signal.  Then, by taking the absolute value of the 
 *        analytic signal the Hilbert transform is obtained.
 * @sa Hilbert
 */
class EnvelopeAnalytic
{
public:
    /*!
     * @name Constructors
     */
    /*!
     * @brief Default constructor.
     */
    EnvelopeAnalytic(); 
    /*!
     * @brief Copy constructor.
     * @param[in] envelope  Class from which to initialize.
     */
    EnvelopeAnalytic(const EnvelopeAnalytic &envelope);
    /*! @} */

    /*! @name Destructors
     */
    /*!
     * @brief Destructor.
     */
    ~EnvelopeAnalytic();
    /*!
     * @brief Resets the module and releases all memory.
     */
    void clear() noexcept;
    /*! @} */

    /*! 
     * @brief Checks if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() noexcept;
    /*! 
     * @brief Gets the length of the envelope.
     * @result The length of the envelope.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getTransformLength();

    /*!
     * @brief Initializes the envelope.
     * @param[in] n   The number of samples in the signal whose envelope is
     *                to be computed.
     * @param[in] precision  Defines the underlying precision of the envelope
     *                       computation. 
     * @throws std::invalid_argument if n is not positive.
     */
    void initialize(const int n,
                    const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Computes the envelope of the signal.
     * @param[in] n  The number of samples in the signal to transform.
     *               This cmust equal \c getTransformLength().
     * @throws std::invalid_argument if n is incorrect or x or yupper are NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void transform(const int n, const double x[], double yupper[]);
};

}
}
}

#endif
