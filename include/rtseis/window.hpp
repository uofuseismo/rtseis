#ifndef RTSEIS_WINDOW_HPP
#define RTSEIS_WINDOW_HPP
#include <memory>
#include <rtseis/vector.hpp>
namespace RTSeis
{
enum class WindowType
{
    Bartlett, /*!< Barlett (triangle) window */
    Blackman, /*!< Standard Blackman */
    Hamming,  /*!< Hamming window */
    Hanning,  /*!< Hann(int) window */
    Kaiser,   /*!< Kaiser window.  This requires specification of \f$ \beta \f$. */
    Sine      /*!< Sine window. */
};
template<class T = double>
/// @class Window 
/// @brief Defines a window function to a signal.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Window
{
public:
    enum class Type
    {
        Bartlett, /*!< Barlett (triangle) window */
        Blackman, /*!< Standard Blackman */
        Hamming,  /*!< Hamming window */
        Hanning,  /*!< Hann(int) window */
        Kaiser,   /*!< Kaiser window.  This requires specification of \f$ \beta \f$. */
        Sine      /*!< Sine window. */
    };
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Window();
    /// @brief Copy constuctor.
    /// @param[in] window  The window from which to initialize this class.
    Window(const Window &window);
    /// @brief Move constructor.
    /// @param[in,out] window  The window from which to initialize this class.
    ///                        On exit, window's behavior is undefined.
    Window(Window &&window) noexcept;
    /// @}

    /// @brief Initializes the window function.
    /// @param[in] length  The number of samples in the window.
    /// @param[in] beta    For Kaiser window design only; this defines the
    ///                    trade-off between the main-lobe and side-lobe
    ///                    level.  As beta increases the window narrows 
    ///                    which in turn requires more samples for 
    ///                    numerical stability.
    /// @throws std::invalid_argument if the length is not positive or,
    ///         for Kaiser window design, beta is negative. 
    void initialize(int length, WindowType type, double beta = 0.5);
    /// @result True indicates the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @result The window type.
    /// @throw std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] WindowType getType() const;

    /// @result The window.
    /// @throws std::runtime_error \c isInitialized() is false.
    [[nodiscard]] Vector<T> getWindow() const;
    /// @result A reference to the underlying window. 
    /// @throws std::runtime_error \c isInitialized() is false.
    /// @note This exists for performance reasons.  You should prefer
    ///       \c getWindow().
    [[nodiscard]] const Vector<T> &getWindowReference() const;

    /// @name Destructors
    /// @{

    /// @brief Resets the class.
    void clear() noexcept;
    /// @brief Constructor.
    ~Window();
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] window  The window to copy to this.
    /// @result A deep copy of the input window.
    Window& operator=(const Window &window);
    /// @brief Move assignment.
    /// @param[in,out] window  The window whose memory will be moved to this.
    ///                        On exit, window's behavior is undefined.
    /// @result The memory from window moved to this.
    Window& operator=(Window &&window) noexcept;
    /// @}
private:
    class WindowImpl;
    std::unique_ptr<WindowImpl> pImpl;
};
}
#endif
