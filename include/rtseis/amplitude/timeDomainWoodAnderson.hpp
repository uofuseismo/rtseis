#ifndef RTSEIS_AMPLITUDE_TIMEDOMAINWOODANDERSON_HPP
#define RTSEIS_AMPLITUDE_TIMEDOMAINWOODANDERSON_HPP 1
#include <memory>
#include "rtseis/enums.hpp"
namespace RTSeis::Amplitude
{
class TimeDomainWoodAndersonRealTimeParameters;
class TimeDomainWoodAndersonPostProcessingParameters;
/// @defgroup rtseis_amplitude Amplitudes
/// @brief These are the core real-time and post-processing 
///        filter implementations to aid in computing amplitudes.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis

/// @class TimeDomainWoodAnderson timeDomainWoodAnderson.hpp "rtseis/amplitude/timeDomainWoodAnderson.hpp"
/// @brief This is an implementation of the time-domain Wood-Anderson
///        filtering described by Kanamori et al.'s 1999 paper:
///        Continuous monitoring of ground-motion parameters.
///        The useful equations are Equation 3 and, for high-pass filtering
///        in advance, Equation 10.
/// @note This class is limited to instruments whose responses are proportional
///       to acceleration or velocity.  Moreover, only a handful of sampling
///       rates have optimized coefficients.  Hence, only a handful of sampling
///       rates are supported.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_amplitude
template<RTSeis::ProcessingMode E = RTSeis::ProcessingMode::POST,
         class T = double>
class TimeDomainWoodAnderson
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    TimeDomainWoodAnderson();
    /// @brief Copy constructor.
    /// @param[in] wa   The Wood Anderson class from which to initialize
    ///                 this class.
    TimeDomainWoodAnderson(const TimeDomainWoodAnderson &wa);
    /// @brief Move constructor.
    /// @bparam[in,out] wa  The Wood Anderson class from which to initialize
    ///                     this class.  On exit, wa's behavior is undefined.
    TimeDomainWoodAnderson(TimeDomainWoodAnderson &&wa) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] wa  The Wood-Anderson class to copy to this.
    /// @result A deep copy of the wa class.
    TimeDomainWoodAnderson& operator=(const TimeDomainWoodAnderson &wa);
    /// @brief Move assignment operator.
    /// @param[in,out] wa   The Wood-Anderson class to move to this.
    ///                     On exit, wa's behaveior is undefined.
    /// @result The memory from wa moved to this.
    TimeDomainWoodAnderson& operator=(TimeDomainWoodAnderson &&wa) noexcept;
    /// @}

    /// @name Initialization
    /// @{
    /// @brief Initializes the Wood-Anderson time-domain filter.
    /// @param[in] parameters   The Wood-Anderson filtering parameters.
    /// @throws std::invalid_argument if the sampling rate is not set.
    void initialize(const TimeDomainWoodAndersonParameters &parameters);
    /// @result The Wood-Anderson parameters.
    /// @result True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @result True indicates this is configured for processing signals
    ///         that are proportional to ground velocity.
    ///         False indicates this is configured for pcoessing signals
    ///         that are proportional to ground acceleration.
    [[nodiscard]] bool isVelocityFilter() const;
    /// @}

    /// @result The length of the initial conditions for the acceleration
    ///         filter.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] int getAccelerationInitialConditionsLength() const;
 
    /// @brief This sets the initial conditions for the differentiator to
    ///        first apply to the acceleration siganl prior to converting
    ///        to a Wood-Anderson response. 
    /// @param[in] nz   The number of coefficients in zi.  This should match
    ///                 \c getAccelerationInitialConditionsLength().
    /// @param[in] zi   The initial conditions.  This is an array of dimension 
    ///                 [nz].  
    /// @throws std::runtime_error if \c isInitialized() is false or if 
    ///         \c isVelocityFilter() is true.
    /// @throws std::invalid_argument if nz is the wrong size or zi is NULL.
    void setAccelerationInitialConditions(int nz, const double zi[]);
    void setWoodAndersonInitialConditions(int nz, const double zi[]); 
    /// @}

    /// @brief Applies the filter to the signal.
    /// @throws std::runtime_error if \c isInitialized() is false.
    void apply(int n, const T x[], T *y[]); 

    /// @name Destructors
    /// @{
    /// @brief Resets the class and releases all memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~TimeDomainWoodAnderson();
    /// @}
private:
    class TimeDomainWoodAndersonImpl;
    std::unique_ptr<TimeDomainWoodAndersonImpl> pImpl;
};
}
#endif
