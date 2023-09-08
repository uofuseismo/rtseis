#ifndef RTSEIS_SYSTEM_SYSTEM_HPP
#define RTSEIS_SYSTEM_SYSTEM_HPP
#include <memory>
namespace RTSeis
{
template<class T> class Vector;
}
namespace RTSeis::System
{
template<class U, class T>
/// @brief A system, H, is an arbitrary unit that takes an input signal, x, of
///        type U and produces an output signal, y, of type T - i.e.,  
///        \f[
///            y(t) = H[x(t)]
///        \f]
///        Systems can be time invariant, linear-time invariant, non-linear etc.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class ISystem
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    ISystem();
    /// @brief Copy constructor.
    ISystem(const ISystem &system);
    /// @brief Move constructor.
    ISystem(ISystem &&system) noexcept;
    /// @}

    /// @name Operators
    /// @{

    ISystem& operator=(const ISystem &system);
    ISystem& operator=(ISystem &&system) noexcept;
    /// @}

    /// @brief Resets the class.
    virtual void clear() noexcept;
    /// @brief Destructor.
    virtual ~ISystem();

    /// @name Initialization
    /// @{

    /// @result True indicates the system is initialized and ready for
    ///         application.
    [[nodiscard]] virtual bool isInitialized() const noexcept = 0;

    /// @name Input Signal
    /// @{

    /// @brief Sets the input signal, x.
    /// @param[in] x   The input signal.
    virtual void setInput(const Vector<U> &x);
    /// @brief Sets the input signal, x.
    /// @param[in,out] x  The input signal to set.  On exit, x's behavior is
    ///                   undefined.
    virtual void setInput(Vector<U> &&x) noexcept;
    /// @result The input signal, x.
    virtual Vector<U> getInput() const;
    /// @result A reference to the input signal, x.
    virtual const Vector<U>& getInputReference() const noexcept;
    /// @}

    /// @name Apply System
    /// @{

    /// @brief Applies the system to the input.
    /// @throws std::runtime_error if \c isInitialized() is false.
    virtual void apply() = 0;
    /// @}

    /// @name Output
    /// @{

    /// @brief Sets the output signal.
    /// @param[in] y   The output signal.
    virtual void setOutput(const Vector<T> &y);
    /// @brief Sets the output signal with move semantics.
    /// @param[in,out] y  The output signal to set.  On exit, y's behavior is
    ///                   undefined.
    virtual void setOutput(Vector<T> &&y) noexcept;
    /// @result The output signal, y.
    virtual Vector<T> getOutput() const;
    /// @result A reference to the output signal, y.
    virtual const Vector<T>& getOutputReference() const noexcept;
    /// @}

    /// @brief Creates a clone of this class.
    //virtual ISystem* clone() const = 0;
private:
    class ISystemImpl;
    std::unique_ptr<ISystemImpl> pImpl;
};
}
#endif
