#ifndef RTSEIS_POSTPROCESSING_SINGLECHANNEL_TAPER_HPP
#define RTSEIS_POSTPROCESSING_SINGLECHANNEL_TAPER_HPP 1
#include <memory>
#include "rtseis/enums.hpp"
namespace RTSeis::PostProcessing::SingleChannel
{
/// @defgroup rtseis_postprocessing_sc_taper Taper
/// @brief Utilities for tapering a signal.  This can be useful for mitigating
///        wraparound artificats prior to application of Fourier methods.
/// @ingroup rtseis_postprocessing_sc

/// @class TaperParameters taper.hpp "include/rtseis/processing/singleChannel/taper.hpp"
/// @brief Defines the parameters for tapering module.
/// @ingroup rtseis_postprocessing_sc_taper
/// @copyright Ben Baker distributed under the MIT license.
class TaperParameters
{
public:
    /// @brief Defines the supported tapers types that can be applied
    ///        to a signal.
    enum Type
    {
        HAMMING,   /*!< Apply Hamming window to ends of signal. */
        BARTLETT,  /*!< Apply Bartlett window to ends of signal. */
        HANN,      /*!< Apply Hann window to ends of signal. */
        BLACKMAN,  /*!< Apply Blackman window to ends of signal. */
        SINE       /*!< Apply a sine window to the ends of the signal. */
    };
/*
public:
    /// @name Constructors
    /// @{
    /// @brief Default constructor.
    /// @param[in] pct   The percentage of the signal to which the taper
    ///                  will be applied.  For example, if 5 percent then
    ///                  the initial 2.5 percent and final 2.5 percent of 
    ///                  the signal will be tapered with the given window.
    ///                  This must be in the range [0,100].
    /// @param[in] type  The window type.
    /// @throws std::invalid_argument if parameters are incorrect.
    TaperParameters(double pct = 5,
                    Type type = Type::HAMMING);
    /// @brief Copy constructor.
    /// @param[in] parms  Taper parameters class from which to initialize
    ///                   this class.
    TaperParameters(const TaperParameters &parms);
    /// @brief Move constructor.
    /// @param[in,out] parms  Taper parameters class to move to this class.
    ///                       On exit this class will no longer be usable.
    TaperParameters(TaperParameters &&parms) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Assignment operator.
    /// @param[in] parms  Taper parameters class to copy to this.
    /// @result A deep copy of the taper parameters class.
    TaperParameters& operator=(const TaperParameters &parms);
    /// @brief Move operator.
    ///@param[in,out] parms  Taper parameters whose memory will be moved 
    ///                      to this.  On exit, params's behavior is undefined.
    TaperParameters& operator=(TaperParameters &&parms) noexcept;
    /// @}
 
    /// @brief Sets the taper type.
    /// @param[in] type   The taper type.
    void setTaperType(Type type) noexcept;
    /// @brief Gets the taper type.
    /// @result The taper type to apply to the signal.
    [[nodiscard]] Type getTaperType() const;

    /// @brief Sets the percentage of the signal to taper.
    /// @param[in] pct   The percentage of the signal to which the taper
    ///                  will be applied.  For example, if 5 percent then
    ///                  the initial 2.5 percent and final 2.5 percent of 
    ///                  the signal will be tapered with the given window.
    ///                  This must be in the range [0,100].
    /// @throws std::invalid_argument if pct is out of range. 
    /// @note This differs from SAC which uses a fraction in the range
    ///       [0,0.5].  This and SAC are related by pct = 100*(frac/2).
    ///       So to taper 20 percent of the signal in SAC one would use
    ///       frac=0.2 which corresponds to pct=10.
    void setPercentage(double pct);
    /// @brief Gets the percentage of the signal to taper.
    /// @result The percentage of the signal to which the taper will
    ///         be applied.
    [[nodiscard]] double getPercentage() const;
    /// @result True indicates that this is a correctly initialized
    ///         parameter class.
    [[nodiscard]] bool isValid() const noexcept;

    /// @name Destructors
    /// @{
    /// @brief Clears variables in class and restores defaults.
    void clear() noexcept;
    /// @brief Default destructor.
    ~TaperParameters();
    /// @}
private:
    class TaperParametersImpl;
    std::unique_ptr<TaperParametersImpl> pImpl;
*/
}; // End Taper Parameters
}
#endif
