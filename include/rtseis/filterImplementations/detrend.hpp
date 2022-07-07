#ifndef RTSEIS_FILTERIMPLEMENTATIONS_DETREND_HPP
#define RTSEIS_FILTERIMPLEMENTATIONS_DETREND_HPP 1
#include <memory>
#include "rtseis/filterImplementations/enums.hpp"
namespace RTSeis::FilterImplementations
{
/// @class Detrend detrend.hpp "rtseis/filterImplementations/detrend.hpp"
/// @brief Removes the mean or trend from a signal.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_filterImplemenations
template<class T>
class Detrend
{
public:
    /// @name Constructors
    /// @{

    /// @brief Default constructor.
    Detrend();
    /// @brief Copy constructor.
    /// @param[in] detrend  Detrend class from which to initialize this class.
    Detrend(const Detrend &detrend);
    /// @brief Move constructor.
    /// @param[in,out] detrend  The detrend class from which to initialize this
    ///                         class.  On exit, detrend's behavior is
    ///                         undefined.
    Detrend(Detrend &&detrend) noexcept;
    /// @brief Constructs and initializes the detrend class.
    /// @param[in] type  The detrend type.
    explicit Detrend(DetrendType type);
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] detrend  The detrend class to copy.
    /// @result A deep copy of the detrend class.
    Detrend& operator=(const Detrend &detrend);
    /// @brief Move assignment operator.
    /// @param detrend   The detrend class to move to this.  On exit, detrend's
    ///                  behavior is undefined.
    /// @return The memory from detrend moved to this.
    Detrend& operator=(Detrend &&detrend) noexcept;
    /// @}
     

    /// @name Destructors
    /// @{

    /// @brief Destructor
    ~Detrend();
    /// @brief Resets the class.
    void clear() noexcept;
    /// @}

    /// @brief Initializes the class.
    /// @param[in] type   Defines the trend removal strategy.
    void initialize(DetrendType type);
    /// @result True indicates that the class is inititalized.
    [[nodiscard]] bool isInitialized() const noexcept;

    /// @brief Detrends the data.
    /// @param[in] nx  The number of samples in the signal.
    /// @param[in] x   The signal to demean or detrend.  This is an array of
    ///                dimension [nx].
    /// @param[out] y  The demeaned or detrended data.  This is an arra of
    ///                dimension [nx].
    /// @throws std::invalid_argument if nx is positive and x or y is NULL.
    /// @throws std::runtime_error if the class is not initialized.
    /// @note If the detrend type is linear and nx is 1 then only the mean will
    ///       be removed.
    void apply(int nx, const T x[], T *y[]);
private:
    class DetrendImpl;
    std::unique_ptr<DetrendImpl> pImpl;
};

/// @brief Function that removes the mean from the data.
/// @param[in] nx      The number of samples in the signal.
/// @param[in] x       The signal from which to remove the mean.  This is
///                    an array of dimension [nx].
/// @param[out] y      The demeaned version of x.  This is an array of
///                    dimension [nx].
/// @param[out] mean   The mean of the dataset.
/// @throws std::invalid_argument if nx is positive and x or y is NULL.
void removeMean(int nx, const double x[], double *y[], double *mean);
/// @copydoc removeMean
void removeMean(int nx, const float x[], float *y[], float *mean);

/// @brief Function that removes the linear trend from the data.
/// @param[in] nx          The number of samples in the signal.
/// @param[in] x           The signal from which to remove the trend.  This is
///                        an array of dimension [nx].
/// @param[out] y          The detrended version of x.  This is an array of
///                        dimension [nx].
/// @param[out] intercept  The intercept of the best-fitting line.
/// @param[out] slope      The slope of the best fitting line.
/// @throws std::invalid_argument if nx is positive and x or y is NULL.
/// @note If nx is less than 2 then the mean will be removed and the slope
///       set to 0.
void removeTrend(int nx, const double x[], double *y[],
                 double *intercept, double *slope);
/// @copydoc removeTrend
void removeTrend(int nx, const float x[], float *y[],
                 float *intercept, float *slope);

}
#endif
