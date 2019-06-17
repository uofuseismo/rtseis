#ifndef RTSEIS_UTILS_FILTER_MRFIR_HPP
#define RTSEIS_UTILS_FILTER_MRFIR_HPP 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{
namespace FilterImplementations
{
/*!
 * @class MultiRateFIRFilter multiRateFIRFilter.hpp "include/rtseis/utilities/filterImplementations/multiRateFIRFilter.hpp"
 * @brief Implements the multi-rate finite impulse response filters.  
 *        This allows for upsampling and/or downsampling while filtering.
 * @note  Realize that this module behaves slightly differently
 *        than Matlab.  When the upsample factor is greater than 1
 *        then the FIR filter will be implicitly multiplied 
 *        by the upsampling factor.  Compare this with Matlab where
 *        the user must gain the FIR filter prior to calling upfirdn.
 * @ingroup rtseis_utils_filters
 */
class MultiRateFIRFilter
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.  Note, this class is not yet
     *        inititalized and cannot be used.
     */
    MultiRateFIRFilter();
    /*!
     * @brief Copy constructor.
     * @param[in] firmr  Multi-rate filter class from which to 
     *                   initialize this class.
     */
    MultiRateFIRFilter(const MultiRateFIRFilter &firmr);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy operator.
     * @param[in] firmr  The multi-rate filter class to copy.
     * @result A deep copy of hte multi-rate filter.
     */
    MultiRateFIRFilter &operator=(const MultiRateFIRFilter &firmr);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~MultiRateFIRFilter();
    /*! @} */

    /*!
     * @{
     * @brief Initializes the multi-rate filtering.
     * @param[in] upFactor    The upsampling factor.  This will insert
     *                        upFactor - 1 zeros to the signal prior to
     *                        filtering.  This must be positive.  There
     *                        will be no upsampling effect if upFactor
     *                        is 1.
     * @param[in] downFactor  The downsampling factor.  This will remove
     *                        remove downFactor - 1 samples beteween
     *                        each upsampled output point.  This must be
     *                        positive.   There will be no downsampling
     *                        effect if downFactor is 1.
     * @param[in] nb          The number of FIR filter coefficients.
     * @param[in] b           The FIR filter coefficients.  This is
     *                        an array of dimension [nb].
     * @param[in] mode        The processing mode.  By default this
     *                        is for post-processing.
     * @param[in] precision   The precision of the filter.  By default
     *                        this is double precision.
     * @result 0 indicates success.
     */
    int initialize(const int upFactor, const int downFactor,
                   const int nb, const double b[],
                   const RTSeis::ProcessingMode mode
                       = RTSeis::ProcessingMode::POST_PROCESSING,
                   const RTSeis::Precision precision
                       = RTSeis::Precision::DOUBLE);
    /*!
     * @copydoc initialize()
     * @param[in] chunkSize   This is an optional tuning parameter.
     */
    int initialize(const int upFactor, const int downFactor,
                   const int nb, const double b[],
                   const int chunkSize,
                   const RTSeis::ProcessingMode mode 
                       = RTSeis::ProcessingMode::POST_PROCESSING,
                   const RTSeis::Precision precision 
                       = RTSeis::Precision::DOUBLE);
    /* @} */

    /*!
     * @brief Determines if the module is initialized.
     * @retval True indicates that the module is initialized.
     * @retval False indicates that the module is not initialized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Gets the length of the initial conditions array.
     * @result On successful exit this is the length of the initial
     *         conditions array.  On failure this is negative.
     */
    int getInitialConditionLength() const;
    /*!
     * @brief Applies the zero-phase IIR filter to the data.  Note,
     *        the class must be initialized prior to using this
     *        function.
     * @param[in] nz  The length of the initial conditions array.  This
     *                can be determined by first calling
     *                getInitialConditionLength().
     * @param[in] zi  The initial conditions.  This has dimension [nz].
     * @param[in] upPhase    The initial upsampling phase.  By default
     *                       this is 0.
     * @param[in] downPhase  The initial downsampler phase.  By default
     *                       this is 0.
     * @result 0 indicates success.
     */
    int setInitialConditions(const int nz,
                             const double zi[],
                             const int upPhase = 0,
                             const int downPhase = 0);
    /*!
     * @brief Sets the initial conditions.
     * @param[in] nz   The length of the initial conditions.  This
     *                 should equal getInitialConditionLength().
     * @param[in] zi   The initial conditions to set.  This has
     *                 has dimension [nwork_]. 
     * @result 0 indicates success.
     */
    int setInitialConditions(const int nz, const double zi[]);
    /*!
     * @brief Estimates the space required to store the output signal.
     * @param[in] n  The length of the input signal.
     * @result The array length required to store the output signal.
     *         If negative than error has occurred.
     */
    int estimateSpace(const int n) const;
    /*!
     * @brief Applies the multi-rate filter to the signal.
     * @param[in] n       Number of points in the signal.
     * @param[in] x       The signal to filter.  This has dimension [n].
     * @param[in] nywork  The workspace reserved to y.  An estimate 
     *                    can be obtained from estimateSpace().
     * @param[out] ny     The length of the output signal.
     * @param[out] y      The filtered signal.  This has dimension
     *                    [nywork] however only the first [ny] samples are set.
     * @result 0 indicates success.
     */
    int apply(const int n, const double x[],
              const int nywork, int *ny, double *y[]);
    /*! @copydoc apply */
    int apply(const int n, const float x[],
              const int nywork, int *ny, float *y[]);

    /*!
     * @brief Resets the initial conditions to those set in
     *        setInitialConditions().
     * @result 0 indicates success.
     */
    int resetInitialConditions();
    /*!
     * @brief Clears memory and resets the filter.  After applying
     *        this function the filter must be re-initialized prior
     *        to being applied to the data.
     */
    void clear() noexcept;
private:
    class MultiRateFIRImpl;
    std::unique_ptr<MultiRateFIRImpl> pFIR_;
}; // End multiratefir
}; // End filterimplementations
}; // End utilities
}; // End rtseis
#endif
