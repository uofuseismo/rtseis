#ifndef RTSEIS_AMPLITUDE_TAUP_HPP
#define RTSEIS_AMPLITUDE_TAUP_HPP
#include <memory>
#include "rtseis/enums.hpp"
namespace RTSeis::Amplitude
{
class TauPParameters;
}
namespace RTSeis::Amplitude
{
/// @class TauP "tauP.hpp" "rtseis/amplitude/tauP.hpp"
/// @brief Computes the \f$ \tau_p \f$ which aims to find the predominant period
///        of an arrival.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<RTSeis::ProcessingMode E = RTSeis::ProcessingMode::POST,
         class T = double>
class TauP
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    TauP();
    /// @}

    /// @name Initialization
    /// @{

    /// @brief Initializes the TauP time-domain filter.
    /// @param[in] parameters  The TauP parameters.
    /// @throws std::invalid_argument if the sampling rate, input units,
    ///         or simple response is not set.
    void initialize(const TauPParameters &parameters);
    /// @result True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Destructor.
    ~TauP();
    /// @}
private:
    class TauPImpl;
    std::unique_ptr<TauPImpl> pImpl;
};
}
#endif
