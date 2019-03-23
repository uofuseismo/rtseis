#ifndef RTSEIS_POSTPROCESSING_SC_WAVEFORM
#define RTSEIS_POSTPROCESSING_SC_WAVEFORM 1
#include <memory>
#include <vector>
#include <string>
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
 * @brief Defines the IIR filter implementation.
 * @ingroup rtseis_postprocessing_sc
 */
enum class IIRFilterImplementation
{
    SOS,    /*!< Apply the filter as a cascade of second order sections.
                 This is numerically more stable than a direct form
                 implementation. */
    DIRECT  /*!< Direct form IIR implementation.  The design and filter
                 application is slightly faster than using second order
                 sections filtering. */
};
/*!
 * @brief Defines the analog prototype from which the IIR filters are designed.
 * @ingroup rtseis_postprocessing_sc
 */
enum class IIRPrototype
{
    BUTTERWORTH, /*!< Butterworth filter design. */
    BESSEL,      /*!< Bessel filter design. */
    CHEBYSHEV1,  /*!< Chebyshev I filter design. */
    CHEBYSHEV2   /*!< Chebyshev II filter design. */
};
/*!
 * @brief Defines the window used in the FIR filter design.
 * @ingroup rtseis_postprocessing_sc 
 */
enum class FIRWindow
{
    HAMMING,     /*!< Hamming window. */
    BARTLETT,    /*!< Bartlett (triangle) window. */
    HANN,        /*!< Hann window. */
    BLACKMAN_OPT /*!< Optimal Blackman window. */
}; 
/*!
 * @brief Defines the filter passband.
 * @ingroup rtseis_postprocessing_sc
 */
enum class Bandtype
{
    LOWPASS,  /*!< Lowpass filter. */
    HIGHPASS, /*!< Highpass filter. */
    BANDPASS, /*!< Bandpass filter. */
    BANDSTOP  /*!< Bandstop (notch) filter. */ 
};

/*!
 * @brief Defines the nature of the convolution or correlation
 *        and the consequence with respect to edge effects.
 * @ingroup rtseis_postprocessing_sc
 */
enum class ConvolutionMode
{
    FULL,   /*!< A full discrete convolution or correlation of
                 inputs which will have length \f$ m + n - 1 \f$.
                 Because the signals do not overlap completely
                 at the convolution edges boundary effects can be
                 seen. */
    VALID,  /*!< The output consists only of those elements that
                 do not rely on zero-padding.  The return 
                 convolution or correlation will have length
                 \f$ \max(m,n) - \min(m,n) + 1 \f$.  This will
                 only be computed where the input signals completely
                 overlap so that there will not be edge effects. */
    SAME,  /*!< The output is the same size as the first input
                 and centered with respect to the FULL output.
                 The resulting convolution or correlation will
                 have length \f$ \max(m, n) \f$. */
};
/*!
 * @brief Defines the convolution or correlation implementation.
 * @ingroup rtseis_postprocessing_sc
 */
enum class ConvolutionImplementation
{
    AUTO,   /*!< Let the implementation decide. */
    DIRECT, /*!< Time domain implementation. */
    FFT     /*!< Frequency domain implementaiton. */
};
/*!
 * @class Waveform Waveform "include/rtseis/processing/singleChannel/postProcessing.hpp"
 * @brief This class is to be used for single-channel post-processing
 *        applications.
 * @ingroup rtseis_postprocessing_sc
 */
class Waveform
{
public:
    /*!
     * @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.  The sampling period defaults to unity.
     */
    Waveform(void); //const double dt = 1);
    /*!
     * @brief Constructs a waveform from time series data.  The sampling
     *        period will be unity.
     * @param[in] x   Signal from which to construct time series.
     * @throws std::invalid_argument if x is empty or the sampling period
     *         is not positive.
     */
    //explicit Waveform(const std::vector<double> &x); //, const double dt = 1);
    /*!
     * @brief Constructs a waveform from time series data.
     * @param[in] dt   Sampling period in seconds.  This must be positive.
     * @param[in] x    Signal from which to construct time series.
     * @throws std::invalid_argument if the sampling period is not positive
     *         or x is empty.
     */
    //Waveform(const double dt, const std::vector<double> &x);
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

    /*! @name Convolution and Correlation
     * @{
     */
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
         const ConvolutionMode mode = ConvolutionMode::FULL,
         const ConvolutionImplementation implementation = ConvolutionImplementation::AUTO);
    /*! 
     * @brief Computes the correlation \f$ x \star s \f$ where the
     *        correlation sum is defined by
     *        \f$ y[k] = \sum_n x[n] s[n+k] \f$.
     * @param[in] s     The signal to correlate with x.
     * @param[in] mode  Defines the correlation output.
     * @param[in] implementation  Defines the implementation type.
     * @throws std::invalid_argument if s is empty or there is no data.
     */
    void correlate(const std::vector<double> &s,
         const ConvolutionMode mode = ConvolutionMode::FULL,
         const ConvolutionImplementation implementation = ConvolutionImplementation::AUTO);
    /*! 
     * @brief Computes the autocorrelation \f$ x \star x \f$.
     * @param[in] mode  Defines the correlation output.
     * @param[in] implementation  Defines the implementation type.
     * @throws std::invalid_argument if there is no data.
     */
    void autocorrelate(const ConvolutionMode mode = ConvolutionMode::FULL,
         const ConvolutionImplementation implementation = ConvolutionImplementation::AUTO);
    /*! @} */

    /*! @name Demeaning and Detrending
     * @{
     */
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
    /*! @} */

    /*! @name Finite Impulse Response Filtering
     * @{
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
     * @note It is the responsibility of the user to ensure that the
     *       signal sampling rate and the sampling rate used in the digital
     *       filter design are compatible.
     */
    void firFilter(const Utilities::FilterRepresentations::FIR &fir,
                   const bool lremovePhase=false);
    /*! 
     * @brief Lowpass filters a signal using an FIR filter.
     * @param[in] ntaps   The number of filter taps.  This must be at least 4.
     *                    Moreover, if the phase shift is to be removed then
     *                    ntaps must be odd.  If it is not odd then ntaps will
     *                    bet set to ntaps + 1.
     * @param[in] fc      The critical frequency in Hz.  This must be between
     *                    0 and and the Nyquist frequency.  The latter can be
     *                    obtained from \c getNyquistFrequency().
     * @param[in] window  The window using the FIR filter design. 
     * @param[in] lremovePhase  If true then the phase shift is to be removed.
     *                          The filters designed have linear group delays
     *                          hence, the first ntaps/2 data points will not
     *                          be saved to the output signal.
     * @throws std::invalid_argument if the number of taps or critical frequency
     *         is invalid.
     */
    void firLowpassFilter(const int ntaps, const double fc, 
                          const FIRWindow window,
                          const bool lremovePhase=false);
    /*! 
     * @brief Highpass filters a signal using an FIR filter.
     * @param[in] ntaps   The number of filter taps.
     * @param[in] fc      The critical frequency in Hz.
     * @param[in] window  The window using the FIR filter design. 
     * @param[in] lremovePhase  If true then the phase shift is to be removed.
     * @throws std::invalid_argument if the number of taps or critical frequency
     *         is invalid.
     * @sa firLowpassFilter()
     */
    void firHighpassFilter(const int ntaps, const double fc,
                           const FIRWindow window,
                           const bool lremovePhase=false);
    /*! 
     * @brief Bandpass filters a signal using an FIR filter.
     * @param[in] ntaps   The number of filter taps.
     * @param[in] fc      The critical frequencies in Hz.  fc.first is the low
     *                    corner and fc.second is the high corner.
     *                    Both critical frequencies must be greater than 0 and
     *                    less than the Nyquist frequency.  Additionally,
     *                    fc.second must be greater than fc.first.  
     * @param[in] window  The window using the FIR filter design. 
     * @param[in] lremovePhase  If true then the phase shift is to be removed.
     * @throws std::invalid_argument if the number of taps or critical frequency
     *         is invalid.
     * @sa firLowpassFilter()
     */
    void firBandpassFilter(const int ntaps, const std::pair<double,double> fc,
                           const FIRWindow window,
                           const bool lremovePhase=false);
    /*! 
     * @brief Bandstop (notch) filters a signal using an FIR filter.
     * @param[in] ntaps   The number of filter taps.
     * @param[in] fc      The critical frequencies in Hz.  fc.first is the low
     *                    corner and fc.second is the high corner.
     *                    Both critical frequencies must be greater than 0 and
     *                    less than the Nyquist frequency.  Additionally,
     *                    fc.second must be greater than fc.first.  
     * @param[in] window  The window using the FIR filter design.
     * @param[in] lremovePhase  If true then the phase shift is to be removed.
     * @throws std::invalid_argument if the number of taps or critical frequency
     *         is invalid.
     * @sa firBandpassfilter
     */
    void firBandstopFilter(const int ntaps, const std::pair<double,double> fc,
                           const FIRWindow window,
                           const bool lremovePhase=false);
    /*! @} */

    /*! @name Second Order Section (Biquadratic) Filtering
     * @{
     * @note SOS filter application can be numerically more robust than
     *       its direct form counterpart and should therefore be preferred
     *       whenever possible.
     */
    /*! 
     * @brief Applies the digital IIR filter represented as cascaded second
     *        order sections to the time series.
     * @param[in] sos  The digital IIR filter stored in second order sections.
     * @param[in] lremovePhase  If true, then this removes the phase distortion
     *                          by applying the filter to both the time forward
     *                          and time reversed signal.
     * @throws std::invalid_argument if the filter is invalid.
     * @note It is the responsibility of the user to ensure that the
     *       signal sampling rate and the sampling rate used in the digital
     *       filter design are compatible.
     */
    void sosFilter(const Utilities::FilterRepresentations::SOS &sos,
                   const bool lzeroPhase=false);
    /*! 
     * @brief Lowpass filters a signal using an IIR filter specified as a series
     *        of second order sections.
     * @param[in] order   The filter order which equals the number of npoles.
     *                    This must be positive.
     * @param[in] fc      The critical frequency in Hz.  This must be between
     *                    0 and and the Nyquist frequency.  The latter can be
     *                    obtained from \c getNyquistFrequency().
     * @param[in] prototype  The analog prototype from which to design the
     *                       the digital lowpass filter.
     * @param[in] ripple  For Chebyshev I filters this controls the maximum
     *                    attenuation in the passband and is measured in dB.
     *                    For Chebyshev II filters this controls the minimum
     *                    attenuation in the stopband and is measured in dB.
     *                    For both Chebyshev I and II filters this must be
     *                    positive.  For Butterworth and Bessel filters this
     *                    parameter is ignored.
     * @param[in] lzeroPhase  If true then the phase shift is removed by
     *                        filtering the signal in both directions.  This
     *                        effectively squares the magnitude of the filter
     *                        response while conveniently annihilating
     *                        the nonlinear phase response.
     * @throws std::invalid_argument if the order isn't positive or the ripple
     *         isn't positive for Chebyshev I or Chebyshev II filter design.
     *
     * @snippet testing/postProcessing/singleChannel.cpp ppSCSOSLowpassExample
     *
     */
    void sosLowpassFilter(const int order, const double fc, 
                          const IIRPrototype prototype,
                          const double ripple,
                          const bool lzeroPhase=false);
    /*!
     * @brief Highpass filters a signal using an IIR filter.
     * @param[in] order   The filter order which equals the number of npoles.
     * @param[in] fc      The critical frequency in Hz.
     * @param[in] prototype  The analog prototype from which to design the
     *                       the digital lowpass filter.
     * @param[in] ripple  Controls the ripple size in Chebyshev I and 
     *                    Chebyshev II design.
     * @param[in] lzeroPhase  If true then the phase shift is removed.
     * @throws std::invalid_argument if the order isn't positive or the ripple
     *         isn't positive for Chebyshev I or Chebyshev II filter design.
     * @sa sosLowpassFilter()
     */
    void sosHighpassFilter(const int order, const double fc,
                           const IIRPrototype prototype,
                           const double ripple,
                           const bool lzeroPhase=false);
    /*! 
     * @brief Bandpass filters a signal using an IIR filter.
     * @param[in] order   The filter order which equals the number of npoles.
     * @param[in] fc      The critical frequencies in Hz.  fc.first is the low
     *                    corner and fc.second is the high corner.
     *                    Both critical frequencies must be greater than 0 and
     *                    less than the Nyquist frequency.  Additionally,
     *                    fc.second must be greater than fc.first.
     * @param[in] prototype  The analog prototype from which to design the
     *                       the digital lowpass filter.
     * @param[in] ripple  Controls the ripple size in Chebyshev I and 
     *                    Chebyshev II design.
     * @param[in] lzeroPhase  If true then the phase shift is removed.
     * @param[in] implementation  Defines the IIR implementation. 
     * @throws std::invalid_argument if the order isn't positive or the ripple
     *         isn't positive for Chebyshev I or Chebyshev II filter design.
     * @sa sosLowpassFilter()
     */
    void sosBandpassFilter(const int order, const std::pair<double,double> fc,
                           const IIRPrototype prototype,
                           const double ripple,
                           const bool lzeroPhase=false);
    /*! 
     * @brief Bandstop filters a signal using an IIR filter.
     * @param[in] order   The filter order which equals the number of npoles.
     * @param[in] fc      The critical frequencies in Hz.
     * @param[in] prototype  The analog prototype from which to design the
     *                       the digital lowpass filter.
     * @param[in] ripple  Controls the ripple size in Chebyshev I and 
     *                    Chebyshev II design.
     * @param[in] lzeroPhase  If true then the phase shift is removed.
     * @param[in] implementation  Defines the IIR implementation. 
     * @throws std::invalid_argument if the order isn't positive or the ripple
     *         isn't positive for Chebyshev I or Chebyshev II filter design.
     * @sa sosLowpassFilter()
     * @sa sosBandpassFilter()
     */
    void sosBandstopFilter(const int order, const std::pair<double,double> fc,
                           const IIRPrototype prototype,
                           const double ripple,
                           const bool lzeroPhase=false);
    /*! @} */

    /*! @name Infinite Impulse Response Filtering
     * @{
     */
    /*! 
     * @brief Applies the digital IIR filter using a direct form implementation
     *        to the time series.
     * @param[in] ba   The digital IIR filter stored as feed-forward and
     *                 feed-back coefficients.
     * @param[in] lremovePhase  If true, then this removes the phase distortion
     *                          by applying the filter to both the time forward
     *                          and time reversed signal.
     * @note For higher-order filters this can be less numerically stable than
     *        an SOS implementation.
     * @throws std::invalid_argument if the filter is invalid.
     */
    void iirFilter(const Utilities::FilterRepresentations::BA &ba,
                   const bool lzeroPhase=false);
    /*!
     * @brief Lowpass filters a signal using an IIR direct form filter.
     * @param[in] order   The filter order which equals the number of npoles.
     *                    This must be positive.
     * @param[in] fc      The critical frequency in Hz.  This must be between
     *                    0 and and the Nyquist frequency.  The latter can be
     *                    obtained from \c getNyquistFrequency().
     * @param[in] prototype  The analog prototype from which to design the
     *                       the digital lowpass filter.
     * @param[in] ripple  Controls the ripple size in Chebyshev I and 
     *                    Chebyshev II design.
     * @param[in] lzeroPhase  If true then the phase shift is removed.
     * @throws std::invalid_argument if the order isn't positive or the ripple
     *         isn't positive for Chebyshev I or Chebyshev II filter design.
     * @sa sosLowpassFilter
     */
    void iirLowpassFilter(const int order, const double fc,
                          const IIRPrototype prototype,
                          const double ripple,
                          const bool lzeroPhase=false);
    /*!
     * @brief Highpass filters a signal using an IIR filter.
     * @param[in] order   The filter order which equals the number of npoles.
     * @param[in] fc      The critical frequency in Hz.
     * @param[in] prototype  The analog prototype from which to design the
     *                       the digital lowpass filter.
     * @param[in] ripple  Controls the ripple size in Chebyshev I and 
     *                    Chebyshev II design.
     * @param[in] lzeroPhase  If true then the phase shift is removed.
     * @throws std::invalid_argument if the order isn't positive or the ripple
     *         isn't positive for Chebyshev I or Chebyshev II filter design.
     * @sa sosLowpassFilter()
     */
    void iirHighpassFilter(const int order, const double fc,
                           const IIRPrototype prototype,
                           const double ripple,
                           const bool lzeroPhase=false);
    /*! 
     * @brief Bandpass filters a signal using an IIR filter.
     * @param[in] order   The filter order which equals the number of npoles.
     * @param[in] fc      The critical frequencies in Hz.
     * @param[in] prototype  The analog prototype from which to design the
     *                       the digital lowpass filter.
     * @param[in] ripple  Controls the ripple size in Chebyshev I and 
     *                    Chebyshev II design.
     * @param[in] lzeroPhase  If true then the phase shift is removed.
     * @param[in] implementation  Defines the IIR implementation. 
     * @throws std::invalid_argument if the order isn't positive or the ripple
     *         isn't positive for Chebyshev I or Chebyshev II filter design.
     * @sa iirLowpassFilter()
     * @sa iirBandpassFilter()
     */
    void iirBandpassFilter(const int order, const std::pair<double,double> fc,
                           const IIRPrototype prototype,
                           const double ripple,
                           const bool lzeroPhase=false);
    /*! 
     * @brief Bandstop filters a signal using an IIR filter.
     * @param[in] order   The filter order which equals the number of npoles.
     * @param[in] fc      The critical frequencies in Hz.
     * @param[in] prototype  The analog prototype from which to design the
     *                       the digital lowpass filter.
     * @param[in] ripple  Controls the ripple size in Chebyshev I and 
     *                    Chebyshev II design.
     * @param[in] lzeroPhase  If true then the phase shift is removed.
     * @param[in] implementation  Defines the IIR implementation. 
     * @throws std::invalid_argument if the order isn't positive or the ripple
     *         isn't positive for Chebyshev I or Chebyshev II filter design.
     * @sa sosLowpassFilter()
     * @sa sosBandpassFilter()
     */
    void iirBandstopFilter(const int order, const std::pair<double,double> fc,
                           const IIRPrototype prototype,
                           const double ripple,
                           const bool lzeroPhase=false);
    /*! @} */

    /*! @name Tapering and Cutting
     * @{
     */
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
    /*! @} */

    /*! @name Utilities
     * @{
     */ 
    /*!
     * @brief Sets the sampling period.
     * @param[in] dt  The signal sampling period in seconds.
     *                This must be positive.
     * @throws std::invalid_argument if dt is not positive.
     */
    void setSamplingPeriod(const double dt);
    /*!
     * @brief Gets the sampling period.
     * @result The signal sampling period in seconds.
     */
    double getSamplingPeriod(void) const noexcept;
    /*!
     * @brief Gets the Nyquist frequency of the signal.
     * @result The Nyquist freuqency in Hz.
     */
    double getNyquistFrequency(void) const noexcept;
    /*! @} */
private:
    class WaveformImpl;
    std::unique_ptr<WaveformImpl> pImpl;
}; // end waveform
}; // End SingleChannel
}; // End PostProcessing
}; // End RTSeis

# endif
