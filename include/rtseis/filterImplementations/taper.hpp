#ifndef RTSEIS_FILTERIMPLEMENTATIONS_TAPER_HPP
#define RTSEIS_FILTERIMPLEMENTATIONS_TAPER_HPP 1
#include <memory>
namespace RTSeis::FilterImplementations
{
/// @brief The supported window types for tapering the ends of a signal.
/// @ingroup rtseis_filterImplemenations
/// @copyright Ben Baker distributed under the MIT license. 
enum class TaperWindowType
{
    Hamming,   /*!< Apply Hamming window to ends of signal. */
    Bartlett,  /*!< Apply Bartlett window to ends of signal. */
    Hann,      /*!< Apply Hann window to ends of signal. */
    Blackman,  /*!< Apply Blackman window to ends of signal. */
    Sine,      /*!< Apply a sine window to the ends of the signal. */
    Kaiser     /*!< Apply a Kaiser window the ends of the signal. */
};
/// @class Taper taper.hpp "rtseis/filterImplementations/taper.hpp"
/// @brief Tapers the ends of a waveform in a manner analogous to SAC.
/// @ingroup rtseis_filterImplemenations
/// @copyright Ben Baker distributed under the MIT license.
template<class T = double>
class Taper
{
public:
    /// @name Constructors
    /// @{ 
    /// @brief Default constructor.
    Taper();
    /// @brief Copy constructor.
    /// @param[in] taper  Taper class from which to initialize this class.
    Taper(const Taper &taper);
    /// @brief Move constructor.
    /// @param[in,out] taper  The taper class from which to initialize this
    ///                       class.  On exit, taper's behavior is undefined. 
    Taper(Taper &&taper) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] taper   The taper class to copy to this.
    /// @result A deep copy of the taper class.
    Taper& operator=(const Taper &taper);
    /// @brief Move assignment operator.
    /// @param[in,out] taper  The taper class whose memory will be moved
    ///                       to this.  On exit, taper's behavior is undefined.
    /// @result The memory from taper moved to this.
    Taper& operator=(Taper &&taper) noexcept;
    /// @}

    /// @name Initialization
    /// @{
    /// @param[in] percentage  The percentage of the signal to which the taper
    ///                        will be applied.  For example, if 5 percent then
    ///                        the initial 2.5 percent and final 2.5 percent of 
    ///                        the signal will be tapered with the given window.
    ///                        This must be in the range [0,100].
    /// @param[in] windowType  The taper window type.
    /// @param[in] beta        The beta for a Kaiser window.  For example: 
    ///                        0 is a rectangular window.
    ///                        5 is similar to a Hamming window.
    ///                        6 is similar to a Hann window.
    ///                        8.6 is similar to a Blackman window.
    /// @throws std::invalid_argument if percentage is out of range or
    ///         the window type is a Kaiser window and beta is negative or
    ///         too big.
    void initialize(double percentage,
                    TaperWindowType windowType = TaperWindowType::Hamming,
                    double beta = 8);
    /// @brief Determines if the class is initialized.
    /// @retval True indicates that the class is ready to be applied to data.
    [[nodiscard]] bool isInitialized() const;
    /// @}

    /// @name Application
    /// @{
    /// @brief Applies the taper to the data.
    /// @param[in] nx   Number of points in the signal.
    /// @param[in] x    The signal to taper.  This has dimension [nx].
    /// @param[out] y   The tapered signal.  This has dimension [nx].
    /// @throws std::invalid_argument if the parameters are invalid. 
    /// @throws std::runtime_error if \c isInitialized() is false.
    void apply(int nx, const T x[], T *y[]);
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Clears the memory and restores the defaults.
    void clear() noexcept;
    /// @brief Default destructor.
    ~Taper();
    /// @} 
private:
    class TaperImpl;
    std::unique_ptr<TaperImpl> pImpl;
}; // End Taper
}
#endif
