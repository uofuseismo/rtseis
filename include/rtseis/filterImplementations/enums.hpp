#ifndef RTSEIS_FILTERIMPLEMENTATIONS_ENUMS_HPP
#define RTSEIS_FILTERIMPLEMENTATIONS_ENUMS_HPP 1
namespace RTSeis::FilterImplementations
{
/// @brief Defines the implementation of the FIR filter.
/// @ingroup rtseis_filterImplemenations
enum class FIRImplementation
{
    DIRECT, /*!< Direct-form implementation.  This is 
                 advantageous for relatively short filters.  */
    FFT,    /*!< FFT overlap and add implementation.  This is
                 advantageous for relatively long filters
                 i.e., when \f$ \log_2 L < N \f$ where
                 \f$ L \f$ is the signal length \f$ N \f$ is
                 the number of taps. */
    AUTO    /*!< The implementation will decide
                 between DIRECT or FFT based. */
};

/// @brief Defines the IIR direct-form implementation.
/// @ingroup rtseis_filterImplemenations
enum IIRDFImplementation
{
    DF2_FAST, /*!< A fast IPP implementation.  For high order
                   filters where poles get close to the unit
                   circle this can generate results that are much
                   different than Matlab.  */
    DF2_SLOW  /*!< A slow implementation that is consistent
                   with Matlab. */
};

/// @brief Defines the detrend strategy.
/// @ingroup rtseis_filterImplemenations
enum DetrendType
{
    CONSTANT, /*!< Removes the mean from the time series. */
    LINEAR    /*!< Removes a best fitting line from the time series. */
};

}
#endif
