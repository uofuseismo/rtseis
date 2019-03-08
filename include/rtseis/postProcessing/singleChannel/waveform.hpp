#ifndef RTSEIS_POSTPROCESSING_SC_WAVEFORM
#define RTSEIS_POSTPROCESSING_SC_WAVEFORM 1
#include <memory>
#include <string>
#include <exception>
#ifndef RTSEIS_UTILS_FILTER_FIR_HPP
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#endif
#ifndef RTSEIS_UTILS_MATH_CONVOLVE_HPP
#include "rtseis/utilities/math/convolve.hpp"
#endif
#ifndef RTSEIS_POSTPROCESSING_SC_TAPER
#include "rtseis/postProcessing/singleChannel/taper.hpp"
#endif

namespace RTSeis
{
// Forward declare filter representations
namespace Utilities
{
namespace FilterRepresentations
{
class FIR;
class SOS;
class BA;
class ZPK;
}; // End filter representations
}; // End utilties
/*!
 * @defgroup rtseis_postprocessing Post-Processing
 * @brief These are higher level algorithms that expedite post-processing
 *        of seismic data. 
 */
namespace PostProcessing
{
/*!
 * @defgroup rtseis_postprocessing_sc Single-Channel Processing
 * @brief This section contains single-channel post-processing algorithms.
 * @ingroup rtseis_postprocessing
 */
namespace SingleChannel
{
/*!
 * @class Waveform Waveform "include/rtseis/processing/singleChannel/postProcessing.hpp"
 * @brief This class is used for single-channel post-processing applications.
 * @ingroup rtseis_postprocessing_sc
 */
class Waveform : public std::exception
{
public:
    /*!
     * @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     * @param[in] dt   Sets the sampling period to 1 second.
     * @throws std::invalid_argument if dt is not positive.
     */
    Waveform(const double dt = 1);
    /*!
     * @brief Constructs a waveform from time series data.  The sampling
     *        period will be unity.
     * @param[in] x   Signal from which to construct time series.
     * @param[in] dt  The sampling period.  By default this is set to 
     *                one sample per second. 
     * @throws std::invalid_argument if x is empty or the sampling period
     *         is not positive.
     */
    explicit Waveform(const std::vector<double> &x, const double dt = 1);
    /*!
     * @brief Constructs a waveform from time series data.
     * @param[in] dt   Sampling period in seconds.
     * @param[in] x    Signal from which to construct time series.
     * @throws std::invalid_argument if the sampling period is not positive
     *         or x is empty.
     */
    Waveform(const double dt, const std::vector<double> &x);
    /*! @} */

    /*!
     * @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */ 
    ~Waveform(void);
    /*! @} */

    /*!
     * @brief Sets a signal on waveform class.
     * @throws std::invalid_argument If x is empty.
     */
    void setData(const std::vector<double> &x);
    /*!
     * @brief Sets a waveform on the module.
     * @param[in] n   The number of points in the signal.
     * @param[in] x   The signal to set.  This is an array of dimension [n].
     * @throws std::invalid_argument if x is null or n is less than 1.
     */
    void setData(const size_t n, const double x[]);
    /*!
     * @brief Gets the processed waveform data.
     * @param[out] y  The processed waveform data.
     */
    void getData(std::vector<double> &y);
    /*!
     * @brief Gets the prcoessed waveform data.
     * @param[in] nwork  Max number of points allocated to y.
     * @param[out] y     The output time series.  This has dimension [nwork]
     *                   however only the first \c getOutputLength() samples
     *                   are accessed.
     * @throws std::invalid_argument if y is NULL or nwork is too small.
     */
    void getData(const size_t nwork, double y[]) const;
    /*!
     * @brief Gets the length of the output signal.
     * @result The length of the output signal, y.
     */
    size_t getOutputLength(void) const;
    /*!
     * @brief Computes the convolution \f$ x \ast s \f$ where the
     *        convolution sum is defined by
     *        \f$ y[k] = \sum_n x[n] s[n-k] \f$.
     * @param[in] s     The signal to convolve with x.
     * @param[in] mode  Defines the convolution output.
     * @param[in] implementation  Defines the implementation type.
     * @throws std::invalid_argument if s is empty or there is no data.
     */
    void convolve(const std::vector<double> &s,
         const Utilities::Math::Convolve::Mode mode = Utilities::Math::Convolve::Mode::FULL,
         const Utilities::Math::Convolve::Implementation implementation = Utilities::Math::Convolve::Implementation::AUTO);
    /*!
     * @brief Removes the mean from the data.
     * @throws std::invalid_argument if there is no data.
     *
     * @snippet testing/postProcessing/singleChannel.cpp ppSCDemeanExample
     */
    void demean(void);
    /*!
     * @brief Removes a best fitting line \f$ \hat{y} = a x + b \f$
     *        from the data by first computing \f$ a \f$ and \f$ b \f$
     *        then computing \f$ y - (a x + b) \f$.
     * @throws std::invalid_argument if there are less than 2 data points  in x.
     *
     * @snippet testing/postProcessing/singleChannel.cpp ppSCDetrendExample
     */
    void detrend(void);

    /*!
     * @name General Filtering
     * @{
     * @note It is the responsibility of hte user to ensure that the
     *       signal sampling rate and the sampling rate used in the digital
     *       filter design are compatible.
     */
    /*!
     * @brief Applies the digital FIR filter to the time series.
     * @param[in] fir  The digital FIR filter to apply to the signal.
     * @param[in] lremovePhase  If true, then this attempts to remove the 
     *                          the phase distortion by applying the
     *                          filter both forwards and backwards.  This is
     *                          done because FIR filters are not required to
     *                          have linear phase responses (i.e., a constant
     *                          group delay).
     * @param[in] implementation  Defines the implementation.
     * @throws std::invalid_argument if the filter is invalid.
     */
    void filter(const Utilities::FilterRepresentations::FIR &fir,
                const bool lremovePhase=false);
                //const Utilities::FilterRepresentations::FIRFilter::Implementation implementation);// = Utilities::FilterRepresentations::FIRFilter::Implementation::DIRECT);
    /*! 
     * @brief Applies the digital IIR filter to the time series.
     * @param[in] sos  The digital IIR filter stored in second order sections.
     * @param[in] lremovePhase  If true, then this removes the phase distortion
     *                          by applying the filter to both the time forward
     *                          and time reversed signal.
     * @note SOS filter application can be numerically more robust than
     *       its direct form counterpart and should therefore be preferred
     *       whenever possible.
     * @throws std::invalid_argument if the filter is invalid.
     */
    void filter(const Utilities::FilterRepresentations::SOS &sos,
                const bool lzeroPhase=false);
    /*! 
     * @brief Applies the digital IIR filter to the time series.
     * @param[in] ba   The digital IIR filter stored as feed-forward and
     *                 feed-back coefficients.
     * @param[in] lremovePhase  If true, then this removes the phase distortion
     *                          by applying the filter to both the time forward
     *                          and time reversed signal.
     * @note For higher-order filters this can be less numerically stable than
     *        an SOS implementation.
     * @throws std::invalid_argument if the filter is invalid.
     */
    void filter(const Utilities::FilterRepresentations::BA &ba,
                const bool lzeroPhase=false);
    /*! @} */

    /*!
     * @brief Tapers the ends of a signal.
     * @param[in] pct  The percentage of the signal to which the taper will
     *                 be applied.  For example, 5 percent indicates that
     *                 the first 2.5 and final 2.5 percent of the signal
     *                 will be tapered.  This must be in the range (0,100). 
     * @param[in] window  Defines the window function used to generate
     *                    the taper.
     * @note The SAC convention would require a fraction in the range (0,0.5).
     *       Hence, to convert from SAC to RTSeis one would do something like
     *       pct = 100*(2*fraction).
     * @throw std::invalid_argument if the parameters are invalid.
     */
    void taper(const double pct = 5,
               const TaperParameters::Type window = TaperParameters::Type::HAMMING);
private:
    class DataImpl;
    std::unique_ptr<DataImpl> pImpl;
}; // end waveform
}; // End SingleChannel
}; // End PostProcessing
}; // End RTSeis

# endif
