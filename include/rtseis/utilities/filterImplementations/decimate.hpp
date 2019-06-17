#ifndef RTSEIS_UTILITIES_FILTER_DECIMATE_HPP
#define RTSEIS_UTILITIES_FILTER_DECIMATE_HPP
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{
namespace FilterImplementations
{

/*!
 * @class Decimate decimate.hpp "include/rtseis/utilities/filterImplementations/decimate.hpp"
 * @brief Lowpass filters then downsamples a signal.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
class Decimate
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Constructor.
     */
    Decimate();
    /*!
     * @brief Copy constructor.
     * @param[in] decimate  The decimation class from which to initialize
     *                      this class.
     */
    Decimate(const Decimate &decimate); 
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] decimate  The decimate class to copy.
     * @result A deep copy of the decimate class.
     */
    Decimate& operator=(const Decimate &decimate);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Destructor
     */
    ~Decimate();
    /*!
     * @brief Resets the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the decimator.
     * @param[in] downFactor         The downsampling factor.  This will retain
     *                               every (downFactor-1)'th sample.  This must
     *                               be positive.
     * @param[in] filterLength       The length of the FIR filter.  This must
     *                               be at least 5.
     * @param[in] lremovePhaseShift  If true then this will remove the phase
     *                               shift introduced by the FIR filter.
     *                               This is relevant when the operation mode
     *                               is for post-processing.
     * @param[in] mode               The processing mode.
     * @param[in] precision          Defines the precision of the underlying
     *                               calculations.
     * @throws std::invalid_argument if the downFactor is not positive, the
     *         filter length is too small.
     * @note This will design an Hamming window-based filter whose cutoff
     *       frequency is 1/downFactor.  Additionally, when post-processing
     *       and removing the phase shift, the algorithm will increase
     *       the filter length so that it's group delay + 1 is evenly
     *       divisible by the downsampling factor.
     */
    void initialize(const int downFactor,
                    const int filterLength = 30,
                    const bool lremovePhaseShift = true,
                    const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                    const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Determines if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() const noexcept;

    /*!
     * @brief Applies the decimator to the data.
     * @param[in] nx       The number data points in x.
     * @param[in] x        The signal to decimate.
     * @param[in] ny       The maximum number of samples in y.  One can
     *                     estimate ny by using estimateSpace(). 
     * @param[out] nyDown  The number of defined decimated points in y.
     * @param[out] y       The decimated signal.  This has dimension
     *                     [ny] however  only the first [nyDown] points
     *                     are defined.
     * @result 0 indicates success.
     * @throws std::invalid_argument if x or y is NULL.
     * @throws std::runtime_error if the module is not initialized.
     */
    void apply(const int nx, const double x[],
               const int ny, int *nyDown, double *y[]);
    /*! @copydoc apply */
    void apply(const int nx, const float x[],
               const int ny, int *nyDown, float *y[]);
};

}
}
}

#endif
