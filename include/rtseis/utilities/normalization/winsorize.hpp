#ifndef RTSEIS_UTILITIES_NORMALIZATION_WINSORIZE_HPP
#define RTSEIS_UTILITIES_NORMALIZATION_WINSORIZE_HPP 1
#include <memory>
namespace RTSeis::Utilities::Normalization
{
/// @class Winsorize winsorize.hpp "include/rtseis/utilities/normalization/winsorize.hpp"
/// @brief Winsorizes the data.  For example, a 90 pct winsorization will set all data
///        below the 5th percentile the 5th percentile and all data above the 95th
///        percentile to the 95th percentile, i.e.,
///        \f[
///            y = \max \left \{x_5, \min \left \{x, x_95 \right \} \right \}
///        \f]
///        where \f$ x_5 \f$ is the 5th percentile of x and \f$ x_{95} \f$ is the
///        95'th percentile of x.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<class T = double>
class Winsorize
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor
    Winsorize();
    /// @brief Copy constructor
    /// @param[in] winsorize  The winsorization class from which to initialize
    ///                       this class.
    Winsorize(const Winsorize &winsorize);
    /// @brief Move constructor
    /// @param[in,out] winsorize  The winsorization class whose memory will
    ///                           be moved to this.  On exit winsorize's
    ///                           behavior is undefined.
    Winsorize(Winsorize &&winsorize) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] winsorize   The winsorization class to copy.
    /// @result A deep copy of winsorize.
    Winsorize& operator=(const Winsorize &winsorize);
    /// @brief Move assignment operator.
    /// @param[in,out] winsorize  The winsorization class whose memory will
    ///                           be moved to this.  On exit, winsorize's
    ///                           behavior will be undefined.
    /// @result The memory from winsorize moved to this.
    Winsorize& operator=(Winsorize &&winsorize) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Destructor
    ~Winsorize();
    /// @brief Resets the class.
    void clear() noexcept;
    /// @}

    /// @brief Initializes the class.
    /// @param[in] limits     The percentiles between which data is undisturbed.
    ///                       Specifically, limits.first is the lower percentile
    ///                       and limits.seconds is the upper percentile.  If
    ///                       there are n samples then n*limits.first/100 will
    ///                       be the smallest resulting Winsorized data and
    ///                       n*(1 - limits.second/100) will be the largest 
    ///                       Winsorized data.
    /// @param[in] inclusive  True indicates that the number of data being 
    ///                       masked should be truncated.  False indicates
    ///                       the number of data being msaked should be rounded.
    /// @param[in] lowMemory  If false then workspace arrays will be allocated
    ///                       during each application.  When working on 
    ///                       packetized data of with small, fixed lengths
    ///                       it is more efficient to set this to false.
    /// @throws std::invalid_argument if limits.first >= limits.second or
    ///         the limits are not between 0 and 100.
    void initialize(const std::pair<double, double> &limits = std::pair(5, 95),
                    bool inclusive = true,
                    bool lowMemory = false);
    /// @result True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @brief Winsorizes the data.
    /// @param[in] nx     The number of points in the time series.
    /// @param[in] x      The input time series to Winsorize.  This is an
    ///                   array of dimension [nx].
    /// @param[out] y     The Winsorized variant of x.  This is an array of
    ///                   dimension [nx].
    /// @throws std::invalid_argument if nx is positive and x or y is NULL.
    /// @throws std::runtime_error if the class is not initialized.
    /// @sa \c isInitialized()
    void apply(int npts, const T x[], T *y[]);
private:
    class WinsorizeImpl;
    std::unique_ptr<WinsorizeImpl> pImpl;
};
}
#endif

