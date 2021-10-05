#ifndef RTSEIS_POSTPROCESSING_SC_TAPER
#define RTSEIS_POSTPROCESSING_SC_TAPER 1
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
}; // End Taper Parameters

/*!
 * @class Taper taper.hpp "include/rtseis/processing/singleChannel/taper.hpp"
 * @brief Tapers a waveform.
 * @ingroup rtseis_postprocessing_sc_taper
 * @copyright Ben Baker distributed under the MIT license.
 */
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
    /// @brief Constructs a taper command from the parameters.
    /// @param[in] parameters  The taper parameters.
    /// @throw std::invalid_argument if the parameters are invalid.
    explicit Taper(const TaperParameters &parameters);
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

    /// @name Destructors
    /// @{
    /// @brief Default destructor.
    ~Taper();
    /// @brief Clears the memory and restores the defaults.
    void clear() noexcept;
    /// @} 

    /// @brief Sets the taper parameters.
    /// @param[in] parameters  A correctly initialized taper parameters
    ///                        class.
    /// @throw std::invalid_argument if the parameters.isValid() is false.
    void setParameters(const TaperParameters &parameters);

    /// @brief Determines if the class is initialized.
    /// @retval True indicates that the class is ready to be applied to data.
    [[nodiscard]] bool isInitialized() const;
    /// @brief Applies the taper to the data.
    /// @param[in] nx   Number of points in the signal.
    /// @param[in] x    The signal to taper.  This has dimension [nx].
    /// @param[out] y   The tapered signal.  This has dimension [nx].
    /// @throws std::invalid_argument if the parameters are invalid. 
    /// @throws std::runtime_error if \c isInitialized() is false.
    void apply(int nx, const T x[], T *y[]);
private:
    class TaperImpl;
    std::unique_ptr<TaperImpl> pImpl;
}; // End Taper
} // End RTSeis
#endif
