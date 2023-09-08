#ifndef RTSEIS_DETREND_HPP
#define RTSEIS_DETREND_HPP
#include <memory>
#include <rtseis/system/system.hpp>
namespace RTSeis
{
template<class T>
/// @class Detrend
/// @brief Removes the linear trend from a signal.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Detrend : public RTSeis::System::ISystem<T, T>
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Detrend();
    /// @brief Copy constuctor.
    Detrend(const Detrend &detrend);
    /// @brief Move constructor.
    Detrend(Detrend &&detrend) noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    Detrend& operator=(const Detrend &detrend);
    /// @brief Move assignment.
    Detrend& operator=(Detrend &&detrend) noexcept;
    /// @}    

    /// @result True indicates the class is initialized which, in this case,
    ///         is always true.
    [[nodiscard]] bool isInitialized() const noexcept override;

    /// @brief Detrends the signal.
    void apply() override;

    /// @result The slope of the best-fitting line.
    [[nodiscard]] double getSlope() const noexcept;
    /// @result The y intercept of the best-fitting line.
    [[nodiscard]] double getIntercept() const noexcept;
 
    /// @name
    /// @{

    /// @brief Resets class and releases memory.
    void clear() noexcept override;
    /// @brief Constructor.
    ~Detrend() override;
    /// @}
private:
    class DetrendImpl;
    std::unique_ptr<DetrendImpl> pImpl;
};
}
#endif
