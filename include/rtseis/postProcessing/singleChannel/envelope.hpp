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
    EnvelopeFIRParameters(EnvelopeFIRParameters &&envfir) noexcept;
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
    EnvelopeFIRParameters &operator=(EnvelopeFIRParameters &&envfir) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~EnvelopeFIRParameters();
    /*!
     * @brief Resets the class.
     */
    void clear() noexcept;
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
    double getBeta() const noexcept;

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
    int getFilterLength() const noexcept;

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
    RTSeis::Precision getPrecision() const noexcept;

    /*!
     * @brief Utility routine to determine if the FIR-based envelope parameters
     *        are correctly set.
     * @result True indicates that the envelope parameters are valid.
     */
    bool isValid() const noexcept;
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
    EnvelopeDFTParameters(EnvelopeDFTParameters &&envdft) noexcept;
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
    EnvelopeDFTParameters &operator=(EnvelopeDFTParameters &&envdft) noexcept;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~EnvelopeDFTParameters();
    /*!
     * @brief Resets the class.
     */
    void clear() noexcept;
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
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.  Note, the class is not yet ready to
     *        be applied to data.
     */
    Envelope();
    /*!
     * @brief Copy constructor.
     * @param[in] envelope  Envelope class from which to initialize.
     */
    Envelope(const Envelope &envelope);
    Envelope(Envelope &&envelope) = delete;
    /*!
     * @brief Initializes the FIR envelope computation.  This uses an FIR
     *        Hilbert transformer.
     * @param[in] firEnvelopeParameters  FIR envelope parameters.
     * @throws std::invalid_argument if the FIR envelope parameters are invalid.
     */
    explicit Envelope(const EnvelopeFIRParameters &firEnvelopeParameters);
    /*!
     * @brief Initializes the analytic envelope computation.  This uses
     *        the Hilbert transform which, in turn, relies on the discrete
     *        Fourier transform (DFT).
     * @param[in] dftEnvelopeParameters  The envelope parameters.
     * @throws std::invalid_argument if the analytic envelope parameters
     *         are invalid.
     */
    explicit Envelope(const EnvelopeDFTParameters &dftEnvelopeParameters);
    /*! @} */

    /*!
     * @brief Copy operator.
     * @param[in] envelope  Envelope class to copy.
     * @result A deep copy of the input envelope class.
     */
    Envelope& operator=(const Envelope &envelope);
    /*!
     * @brief Default destructor.
     */
    ~Envelope();
    /*!
     * @brief Resets the class.  Note, that the parameters must be set again
     *        prior to applying data.
     */
    void clear();

    /*!
     * @brief Sets the FIR-based envelope parameters.
     */
    void setParameters(const EnvelopeFIRParameters &firEnvelopeParameters);
    /*!
     * @brief Sets the DFT-based envelope parameters.
     */
    void setParameters(const EnvelopeDFTParameters &dftEnvelopeParameters);

    /*!
     * @brief Indicates whether or not the class is initialized.
     * @retval True indicates that the class is ready to be applied to data.
     * @retval False indicates that the class is not yet initialized.
     */
    bool isInitialized() const noexcept;

    /*!
     * @brief Computes the envelope of the signal.
     * @param[in] nx   The number of points in the input signal.
     * @param[in] x    The signal from which to compute the envelope.
     *                 This is an array of dimension [nx].
     * @param[out] y   The (upper) envelope of the data.  This is an array
     *                 whose dimension is [nx].
     * @throws std::runtime_error if the class is not yet initialized.
     */
    void apply(const int nx, const double x[], double y[]);
    /*! 
     * @brief Computes the envelope of the signal.
     * @param[in] nx   The number of points in the input signal.
     * @param[in] x    The signal from which to compute the envelope.
     *                 This is an array of dimension [nx].
     * @param[out] y   The (upper) envelope of the data.  This is an array
     *                 whose dimension is [nx].
     * @throws std::runtime_error if the class is not yet initialized.
     */
    void apply(const int nx, const float x[], float y[]);
    /*! 
     * @brief Computes the envelope of the signal.
     * @param[in] nx       The number of points in the input signal.
     * @param[in] x        The signal from which to compute the envelope.
     *                     This is an array of dimension [nx].
     * @param[out] ylower  The lower envelope of the data.  This is an array
     *                     whose dimension is [nx].
     * @param[out] yupper  The upper envelope of the data.  This is an array
     *                     whose dimension is [nx].
     * @throws std::runtime_error if the class is not yet initialized.
     */
    void apply(const int nx, const double x[],
               double ylower[], double yupper[]);
    /*! 
     * @brief Computes the envelope of the signal.
     * @param[in] nx       The number of points in the input signal.
     * @param[in] x        The signal from which to compute the envelope.
     *                     This is an array of dimension [nx].
     * @param[out] ylower  The lower envelope of the data.  This is an array
     *                     whose dimension is [nx].
     * @param[out] yupper  The upper envelope of the data.  This is an array
     *                     whose dimension is [nx].
     * @throws std::runtime_error if the class is not yet initialized.
     */
    void apply(const int nx, const float x[],
               float ylower[], float yupper[]);
private:
    class EnvelopeImpl;
    std::unique_ptr<EnvelopeImpl> pImpl;
};



}; // End SingleChannel
}; // End PostProcessing
}; // End RTSeis
#endif
