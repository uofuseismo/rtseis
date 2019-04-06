#ifndef RTSEIS_POSTPROCESSING_SC_ENVELOPE
#define RTSEIS_POSTPROCESSING_SC_ENVELOPE 1
#include <memory>
#include "rtseis/enums.h"

namespace RTSeis
{
namespace PostProcessing
{
namespace SingleChannel
{

/*!
 * @defgroup rtseis_postprocessing_sc_envelope Envelope
 * @brief Utilities for computing the envelope of a signal.
 * @ingroup rtseis_postprocessing_sc
 */

/*!
 * @class EnvelopeFIRParameters envelope.hpp "include/rtseis/processing/singleChannel/envelope.hpp"
 * @brief Defines the parameters for computing the envelope with an FIR
 *        filter.
 * @ingroup rtseis_postprocessing_sc_envelope
 * @copyright Ben Baker distributed under the MIT license.
 */
class EnvelopeFIRParameters
{
public:
    /*! @name Constructors
     * @{
     */ 
    /*!
     * @brief Default constructor.
     * @param[in] nfir  The number of elements in the FIR filter.  
     *                  This will use \c setFilterLength().
     * @param[in] beta  The \f$ \beta \f$ parameter used in designing the
     *                  window Kaiser-window FIR filter HIlbert transformer.
     * @param[in] precision  The precision of the module.
     * @throws std::invalid_arugment If nfir is invalid.
     */
    EnvelopeFIRParameters(const int nfir = 301,
                          const double beta = 8.0,
                          const RTSeis::Precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Copy constructor.
     * @param[in] envfir  The class from which to initialize this class.
     */
    EnvelopeFIRParameters(const EnvelopeFIRParameters &envfir);
    /*!
     * @brief Move constructor.
     * @param[in,out] envfir  The class to move to this class.  On exit envfir's
     *                        behavior will be undefined.
     */
    EnvelopeFIRParameters(EnvelopeFIRParameters &&envfir);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] envfir  Class to copy.
     * @result A copy of the envfir class.
     */
    EnvelopeFIRParameters &operator=(const EnvelopeFIRParameters &envfir);
    /*!
     * @brief Move assignement operator.
     * @param[in,out] envfir  Class to move.  On exit envfir's behavior will
     *                        be undefined.
     */
    EnvelopeFIRParameters &operator=(EnvelopeFIRParameters &&envfir);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~EnvelopeFIRParameters(void);
    /*!
     * @brief Resets the class.
     */
    void clear(void) noexcept;
    /*! @} */

    /*!
     * @brief Sets \f$ \beta \f$ used in the Kaiser window-based
     *        FIR Hilbert transformer design.
     * @param[in] beta  The \f$ \beta \f$ in the Kaiser window design.
     * @sa \c Utilities::WindowFunctions::kaiser()
     */
    void setBeta(const double beta) noexcept;
    /*!
     * @brief Gets the \f$ \beta \f$ in the FIR Kaiser window design.
     * @result The \f$ \beta \f$ in the window design.
     */
    double getBeta(void) const noexcept;

    /*!
     * @brief Sets the FIR filter length.
     * @param[in] nfir  The number of elements in the FIR filter.  This must
     *                  be positive.  Additionally, it is recommended that this
     *                  be an odd number to increase the filtering performance
     *                  and to simplify removal of the phase shift.
     * @throws std::invalid_argument if nfir is not positive.
     */
    void setFilterLength(const int nfir);
    /*!
     * @brief Gets the FIR filter length.
     * @result The number of taps in the FIR Hilbert transformer filter.
     */
    int getFilterLength(void) const noexcept;

    /*!
     * @brief Defines the precision of the computation.
     * @param[in] precision  The precision which the filter will be applied.
     */
    void setPrecision(const RTSeis::Precision precision) noexcept;
    /*!
     * @brief Gets the precision the filter will be applied.
     * @result The floating arithmetic precision in which the filter will
     *         be applied.
     */
    RTSeis::Precision getPrecision(void) const noexcept;

    /*!
     * @brief Utility routine to determine if the FIR-based envelope parameters
     *        are correctly set.
     * @result True indicates that the envelope parameters are valid.
     */
    bool isValid(void) const noexcept;
private:
    class EnvelopeFIRParametersImpl;
    std::unique_ptr<EnvelopeFIRParametersImpl> pImpl;
}; // End FIRParameters

/*!
 * @class EnvelopeDFTParameters envelope.hpp "include/rtseis/processing/singleChannel/envelope.hpp"
 * @brief Defines the parameters for computing the envelope with a DFT.
 * @ingroup rtseis_postprocessing_sc_envelope
 * @copyright Ben Baker distributed under the MIT license.
 */
class EnvelopeDFTParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     * @param[in] precision  The precision of the module.
     */
    EnvelopeDFTParameters(const RTSeis::Precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Copy constructor.
     * @param[in] envdft  The class from which to initialize this class.
     */
    EnvelopeDFTParameters(const EnvelopeDFTParameters &envdft);
    /*!
     * @brief Move constructor.
     * @param[in,out] envdft  The class to move to this class.  On exit envdft's
     *                        behavior will be undefined.
     */
    EnvelopeDFTParameters(EnvelopeDFTParameters &&envdft);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] envdft  Class to copy.
     * @result A copy of the envdft class.
     */
    EnvelopeDFTParameters &operator=(const EnvelopeDFTParameters &envdft);
    /*!
     * @brief Move assignement operator.
     * @param[in,out] envdft  Class to move.  On exit envdft's behavior will
     *                        be undefined.
     */
    EnvelopeDFTParameters &operator=(EnvelopeDFTParameters &&envdft);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~EnvelopeDFTParameters(void);
    /*!
     * @brief Resets the class.
     */
    void clear(void) noexcept;
    /*! @} */

    /*! 
     * @brief Defines the precision of the computation.
     * @param[in] precision  The precision which the filter will be applied.
     */
    void setPrecision(const RTSeis::Precision precision) noexcept;
    /*! 
     * @brief Gets the precision the filter will be applied.
     * @result The floating arithmetic precision in which the filter will
     *         be applied.
     */
    RTSeis::Precision getPrecision(void) const noexcept;

    /*!
     * @brief Utility routine to determine if the DFT-based envelope parameters
     *        are correctly set.
     * @result True indicates that the envelope parameters are valid.
     */
    bool isValid(void) const noexcept;
private:
    class EnvelopeDFTParametersImpl;
    std::unique_ptr<EnvelopeDFTParametersImpl> pImpl;
};

/*!
 * @class Envelope envelope.hpp "include/rtseis/processing/singleChannel/envelope.hpp"
 * @brief Computes the envelope of the signal.
 * @ingroup rtseis_postprocessing_sc_envelope
 * @copyright Ben Baker distributed under the MIT license.
 */
class Envelope
{

private:
    class EnvelopeImpl;
    std::unique_ptr<EnvelopeImpl> pImpl;
};



}; // End SingleChannel
}; // End PostProcessing
}; // End RTSeis
#endif
