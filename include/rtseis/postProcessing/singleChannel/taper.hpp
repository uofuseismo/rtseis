#ifndef RTSEIS_POSTPROCESSING_SC_TAPER
#define RTSEIS_POSTPROCESSING_SC_TAPER 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace PostProcessing
{
namespace SingleChannel
{

/*!
 * @defgroup rtseis_postprocessing_sc_taper Taper
 * @brief Utilities for tapering a signal.  This can be useful for mitigating
 *        wraparound artificats prior to application of Fourier methods.
 * @ingroup rtseis_postprocessing_sc
 */

/*!
 * @class TaperParameters taper.hpp "include/rtseis/processing/singleChannel/taper.hpp"
 * @brief Defines the parameters for tapering module.
 * @ingroup rtseis_postprocessing_sc_taper
 * @copyright Ben Baker distributed under the MIT license.
 */
class TaperParameters
{
public:
    /*!
     * @brief Defines the supported tapers types that can be applied
     *        to a signal.
     */
    enum Type
    {
        HAMMING,   /*!< Apply Hamming window to ends of signal. */
        BARTLETT,  /*!< Apply Bartlett window to ends of signal. */
        HANN,      /*!< Apply Hann window to ends of signal. */
        BLACKMAN,  /*!< Apply Blackman window to ends of signal. */
        SINE       /*!< Apply a sine window to the ends of the signal. */
    };

public:
    /*!
     * @brief Default constructor.
     * @param[in] pct   The percentage of the signal to which the taper
     *                  will be applied.  For example, if 5 percent then
     *                  the initial 2.5 percent and final 2.5 percent of 
     *                  the signal will be tapered with the given window.
     *                  This must be in the range [0,100].
     * @param[in] type  The window type.
     * @throws std::invalid_argument if parameters are incorrect.
     */ 
    TaperParameters(double pct = 5,
                    const Type type = Type::HAMMING);
    /*!
     * @brief Copy constructor.
     * @param[in] parms  Taper parameters class from which to initialize
     *                   this class.
     */
    TaperParameters(const TaperParameters &parms);
    /*!
     * @brief Move constructor.
     * @param[in,out] parms  Taper parameters class to move to this class.
     *                       On exit this class will no longer be usable.
     */
    TaperParameters(TaperParameters &&parms);
    /*!
     * @brief Assignment operator.
     * @param[in] parms  Taper class to copy.
     * @result A deep copy of the taper parameters class.
     */
    TaperParameters& operator=(const TaperParameters &parms);
    /*!
     * @brief Move operator.
     * @param[in,out] parms  Taper parameters to move.  On exit,
     *                       this class is no longer usable.
     */
    TaperParameters& operator=(TaperParameters &&parms);
    /*!
     * @brief Default destructor.
     */
    ~TaperParameters(void);
    /*!
     * @brief Sets the taper type.
     * @param[in] type   The taper type.
     */
    void setTaperType(const Type type);
    /*!
     * @brief Gets the taper type.
     * @result The taper type to apply to the signal.
     */
    Type getTaperType() const;
    /*!
     * @brief Sets the percentage of the signal to taper.
     * @param[in] pct   The percentage of the signal to which the taper
     *                  will be applied.  For example, if 5 percent then
     *                  the initial 2.5 percent and final 2.5 percent of 
     *                  the signal will be tapered with the given window.
     *                  This must be in the range [0,100].
     * @throws std::invalid_argument if pct is out of range. 
     * @note This differs from SAC which uses a fraction in the range
     *       [0,0.5].  This and SAC are related by pct = 100*(frac/2).
     *       So to taper 20 percent of the signal in SAC one would use
     *       frac=0.2 which corresponds to pct=10.
     */
    void setPercentage(double pct);
    /*!
     * @brief Gets the percentage of the signal to taper.
     * @result The percentage of the signal to which the taper will
     *         be applied.
     */
    double getPercentage() const;
    /*!
     * @brief Determines if the taper parameters are valid.
     * @retval True indicates that this is a correctly initialized
     *         parameter class.
     * @retval False indiates that this is an incorrectly initialized
     *         parameter class.
     */
    bool isValid() const;
    /*!
     * @brief Clears variables in class and restores defaults.
     */
    void clear();
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
template<class T>
class Taper
{
public:
    /*! @name Constructors
     * @{ 
     */
    /*!
     * @brief Default constructor.
     */
    Taper();
    /*!
     * @brief Copy constructor.
     * @param[in] taper   Taper class from which to initialize.
     */
    Taper(const Taper &taper);
    /*!
     * @brief Constructs a taper command from the parameters.
     * @param[in] parameters  The taper parameters.
     * @throw std::invalid_argument if the parameters are invalid.
     */
    explicit Taper(const TaperParameters &parameters);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] Taper   The taper class to copy.
     * @result A deep copy of the demean class.
     */
    Taper& operator=(const Taper &Taper);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~Taper();
    /*! 
     * @brief Clears the memory and restores the defaults.
     */
    void clear();
    /*! @} */

    /*!
     * @brief Sets the taper parameters.
     * @param[in] parameters  A correctly initialized taper parameters
     *                        class.
     * @throw std::invalid_argument if the parameters are invalid.
     */
    void setParameters(const TaperParameters &parameters);

    /*!
     * @brief Determines if the class is initialized.
     * @retval True indicates that the class is ready to be applied to data.
     * @retval False indicates that the class is not ready to be
     *         applied to data.
     */
    bool isInitialized() const;
    /*!
     * @brief Applies the taper to the data.
     * @param[in] nx   Number of points in the signal.
     * @param[in] x    The signal to taper.  This has dimension [nx].
     * @param[out] y   The tapered signal.  This has dimension [nx].
     * @throw std::invalid_argument if the parameters are invalid. 
     */
    void apply(int nx, const T x[], T y[]);
private:
    class TaperImpl;
    std::unique_ptr<TaperImpl> pImpl;
}; // End Taper


}; // End Single channel
}; // End Post-processing
}; // End RTSeis
#endif
