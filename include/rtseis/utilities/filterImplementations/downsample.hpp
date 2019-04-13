#ifndef RTSEIS_UTILS_FILTER_DOWNSAMPLE_HPP
#define RTSEIS_UTILS_FILTER_DOWNSAMPLE_HPP 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{

/*!
 * @defgroup rtseis_utils_filters Filter Implementations
 * @brief These are the core real-time and post-processing 
 *        filter implementations to be used by higher-level modules.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils
 */    
namespace FilterImplementations
{

/*!
 * @class Downsample downsample.hpp "include/rtseis/utilities/filterImplementations/downsample.hpp"
 * @brief This is the core implementation for downsampling a signal.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_filters
 */
class Downsample
{
public:
    /*! @name Constructors
     * @{
     */ 
    /*!
     * @brief Default constructor.
     */
    Downsample(void);
    /*!
     * @brief Copy constructor.
     * @param[in] downsample  Downsampling class from which to initialize.
     */
    Downsample(const Downsample &downsample);
    /*!
     * @brief Move constructor.
     * @param[in,out] downsample  Downsampling class to move to this class.
     *                            On exit downsample's behavior is undefined.
     */
    Downsample(Downsample &&downsample);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] downsample  Downsampling class to copy.
     * @result A deep copy of the downsampling class.
     */
    Downsample& operator=(const Downsample &downsample);
    /*!
     * @brief Move assignment operator.
     * @param[in,out] downsample  On exit this will no longer be usable. 
     * @result A downsample class whose memory was moved from the input class.
     */
    Downsample& operator=(Downsample &&downsample);
    /*! @} */

    /*!
     * @brief Default destructor.
     */
    ~Downsample(void);

    /*!
     * @brief Initializes the downsampler.
     * @param[in] downFactor  The downsampling factor.  This will retain
     *                        every (downFactor-1)'th sample.  This must be
     *                        positive.
     * @param[in] mode  The processing mode.  By default this
     *                  is for post-processing.
     * @param[in] precision   The precision of the filter.  By default
     *                        this is double precision.
     * @throws std::invalid_argument if downFactor is invalid.
     */
    void initialize(const int downFactor,
                   const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                   const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    bool isInitialized(void) const noexcept;
    /*!
     * @brief Estimates the space required to hold the downsampled signal.
     * @param[in] n   The length of the signal to downsample.  This must
     *                be non-negative.
     * @result The number of points required to store the output signal.
     * @throws std::runtime_error if the module is not initialized.
     * @throws std::invalid_argument if n is negative.
     */
    int estimateSpace(const int n) const;
    /*!
     * @brief Gets the downsampling factor.
     * @result The downsampling factor.
     */
    int getDownsampleFactor(void) const noexcept;
    /*!
     * @brief Sets the initial conditions of the downsampler which is the phase.
     * @param[in] phase  Phase of downsampler.  This must be in the 
     *                   range [0, getDownsampleFactor()].
     * @throws std::invalid_argument if phase is out of range.
     * @throws std::runtime_error if the class is not initialized.
     */
    void setInitialConditions(const int phase);
    /*!
     * @{
     * @brief Applies the downsampler to the data.
     * @param[in] nx       The number data points in x.
     * @param[in] x        The signal to downsample.
     * @param[in] ny       The maximum number of samples in y.  One can
     *                     estimate ny by using estimateSpace(). 
     * @param[out] nyDown  The number of defined downsampled points in y.
     * @param[out] y       The downsampled signal.  This has dimension
     *                     [ny] however  only the first [nyDown] points
     *                     are defined.
     * @result 0 indicates success.
     * @throws std::invalid_argument if x or y is NULL.
     * @throws std::runtime_error if the module is not initialized.
     */
    void apply(const int nx, const double x[],
               const int ny, int *nyDown, double y[]);
    void apply(const int nx, const float x[],
               const int ny, int *nyDown, float y[]);
    /*! @} */
    /*!
     * @brief Resets the initial conditions to the phase set in 
     *        setInitialConditions.  If setInitialConditions was not
     *        called then this will set the phase to 0.
     * @throws std::runtime_error if the class is not initialized.
     */
    void resetInitialConditions(void);
    /*! 
     * @brief Clears the module and resets all parameters.
     */
    void clear(void);
private:
    class DownsampleImpl;
    std::unique_ptr<DownsampleImpl> pDownsample_; 
}; // End downsample
}; // End filters
}; // End utilities
}; // End RTSeis
#endif
