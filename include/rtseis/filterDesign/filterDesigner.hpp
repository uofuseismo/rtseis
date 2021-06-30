#ifndef RTSEIS_FILTERDESIGN_FILTERDESIGNER_HPP
#define RTSEIS_FILTERDESIGN_FILTERDESIGNER_HPP
#include <memory>
#include "rtseis/filterDesign/enums.hpp"

// Forward declarations
namespace RTSeis::FilterRepresentations
{
class BA; 
class SOS;
class ZPK;
class FIR;
}
namespace RTSeis::FilterDesign
{
/// @class FilterDesigner filterDesigner.hpp "rtseis/filterDesign/filterDesigner.hpp"
/// @brief A class for filter design.  If designing many filters then using
///        this class may be advantageous as it will save previous filter
///        designs.
/// @copyright Ben Baker distributed under the MIT license.
/// @ingroup rtseis_filterdesign_filterDesigner
class FilterDesigner
{
public:
    /// @name Constructors
    /// @{
    /// @brief Default constructor.
    FilterDesigner(); 
    /// @brief Copy constructor.
    /// @param[in] design  Class from which to initialize this class.
    FilterDesigner(const FilterDesigner &design); 
    /// @brief Move constructor.
    /// @param[in,out] design  The class from which to initialize this class.
    ///                        On exit, design's behavior is undefined.
    FilterDesigner(FilterDesigner &&design) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] design  The class to copy to this.
    /// @result A deep copy of the filter designer class.
    FilterDesigner& operator=(const FilterDesigner &design);
    /// @brief Move assignment operator.
    /// @param[in,out] design  The class whose memory will be moved to this.
    ///                        On exit, design's behavior is undefined.
    /// @result The memory from design moved to this.
    FilterDesigner& operator=(FilterDesigner &&design) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Default destructor.
    ~FilterDesigner();
    /// @brief Erases all existing filter designs and releases all memory.
    void clear() noexcept;
    /// @}

    /// @name FIR Window-Based Filter Design
    /// @{
    /// @brief Designs an FIR lowpass filter.
    /// @param[in] order   The filter order.  The number of taps is order+1.
    ///                    This must be at least 4.
    /// @param[in] r       The normalized cutoff frequency where 1 is the
    ///                    Nyquist frequency.
    /// @param[in] window  The FIR window design.
    /// @param[out] fir    The lowpass filter that corresponds to the 
    ///                    given design parameters.
    /// @throws std::invalid_argument If any of the arguments are invalid.
    void designLowpassFIRFilter(int order,
                                double r,
                                FIRWindow window,
                                FilterRepresentations::FIR &fir) const;
    /// @brief Designs an FIR highpass filter.
    /// @param[in] order   The filter order.  The number of taps is order+1.
    ///                    This must be at least 4.
    /// @param[in] r       The normalized cutoff frequency where 1 is the
    ///                    Nyquist frequency.
    /// @param[in] window  The FIR window design.
    /// @param[out] fir    The highpass filter that corresponds to the 
    ///                    given design parameters.
    /// @throws std::invalid_argument If any of the arguments are invalid.
    void designHighpassFIRFilter(int order,
                                 double r,
                                 FIRWindow window,
                                 FilterRepresentations::FIR &fir) const;
    /// @brief Designs an FIR bandpass filter.
    /// @param[in] order   The filter order.  The number of taps is order+1.
    ///                    This must be at least 4.
    /// @param[in] r       Normalized low cutoff frequency and high cutoff
    ///                    frequencies where 1 is the Nyquist.  Here,
    ///                    r.first is the low cutoff and r.second is the
    ///                    high cutoff.
    /// @param[in] window  The FIR window design.
    /// @param[out] fir    If the bandpass filter that corresponds to the 
    ///                    given design parameters.
    /// @throws std::invalid_argument If any of the arguments are invalid.
    void designBandpassFIRFilter(int order,
                                 const std::pair<double,double> &r,
                                 FIRWindow window,
                                 FilterRepresentations::FIR &fir) const;
    /// @brief Designs an FIR bandstop filter.
    /// @param[in] order   The filter order.  The number of taps is order+1.
    ///                    This must be at least 4.
    /// @param[in] r       Normalized low cutoff frequency and high cutoff
    ///                    frequencies where 1 is the Nyquist.  Here,
    ///                    r.first is the low cutoff and r.second is the
    ///                    high cutoff.
    /// @param[in] window  The FIR window design.
    /// @param[out] fir    If the bandstop filter that corresponds to the 
    ///                    given design parameters.
    /// @throws std::invalid_argument If any of the arguments are invalid.
    void designBandstopFIRFilter(int order,
                                 const std::pair<double,double> &r,
                                 FIRWindow window,
                                 FilterRepresentations::FIR &fir) const;
    /// @}

    /// @name IIR Analog Prototype-Based Design
    /// @{
    /// @brief Designs an IIR lowpass filter.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequency where 1 is the
    ///                   Nyquist frequency.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] zpk    The lowpass filter design.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designLowpassIIRFilter(int n, double r,
                                IIRPrototype ftype,
                                double ripple,
                                FilterRepresentations::ZPK &zpk,
                                IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an IIR lowpass filter.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequency where 1 is the
    ///                   Nyquist frequency.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] ba    The lowpass filter design.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designLowpassIIRFilter(int n, double r,
                                IIRPrototype ftype,
                                double ripple,
                                FilterRepresentations::BA &ba,
                                IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an IIR lowpass filter stored as second-order sections.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequency where 1 is the
    ///                   Nyquist frequency.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] sos   The lowpass filter design.
    /// @param[in] pairing   Defines the pole pairing policy.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designLowpassIIRFilter(int n, double r,
                                IIRPrototype ftype,
                                double ripple,
                                FilterRepresentations::SOS &sos,
                                SOSPairing pairing = SOSPairing::NEAREST,
                                IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);

    /// @brief Designs an  IIR highpass filter.
    /// @param[in] n       The order of the filter.
    /// @param[in] r       The normalized cutoff frequency where 1 is the
    ///                    Nyquist frequency.
    /// @param[in] ftype   The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] zpk    The highpass filter design.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designHighpassIIRFilter(int n, double r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::ZPK &zpk,
                                 IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an  IIR highpass filter.
    /// @param[in] n       The order of the filter.
    /// @param[in] r       The normalized cutoff frequency where 1 is the
    ///                    Nyquist frequency.
    /// @param[in] ftype   The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] ba     The highpass filter design.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designHighpassIIRFilter(int n, double r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::BA &ba,
                                 IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an IIR highpass filter stored as second-order sections.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequency where 1 is the
    ///                   Nyquist frequency.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] sos   The highpass filter design.
    /// @param[in] pairing   Defines the pole pairing policy.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designHighpassIIRFilter(int n, double r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::SOS &sos,
                                 SOSPairing pairing = SOSPairing::NEAREST,
                                 const IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an IIR bandpass filter.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequencies where 1 is the
    ///                   Nyquist frequency.  Here, r.first is the low corner
    ///                   and r.second is the high corner, and it is required
    ///                   that r.second > r.first.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] zpk    The bandpass filter.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designBandpassIIRFilter(int n, const std::pair<double, double> &r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::ZPK &zpk,
                                 IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an IIR bandpass filter.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequencies where 1 is the
    ///                   Nyquist frequency.  Here, r.first is the low corner
    ///                   and r.second is the high corner, and it is required
    ///                   that r.second > r.first.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] ba     The bandpass filter.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designBandpassIIRFilter(int n, const std::pair<double, double> &r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::BA &ba,
                                 IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an IIR bandpass filter stored as second order sections.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequencies where 1 is the
    ///                   Nyquist frequency.  Here, r.first is the low corner
    ///                   and r.second is the high corner, and it is required
    ///                   that r.second > r.first.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] sos   The bandpass filter design.
    /// @param[in] pairing   Defines the pole pairing policy.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designBandpassIIRFilter(int n, const std::pair<double, double> &r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::SOS &sos,
                                 SOSPairing pairing = SOSPairing::NEAREST,
                                 IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);

    /// @brief Designs an IIR bandstop filter.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequencies where 1 is the
    ///                   Nyquist frequency.  Here, r.first is the low corner
    ///                   and r.second is the high corner, and it is required
    ///                   that r.second > r.first.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] zpk    The bandstop filter.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designBandstopIIRFilter(int n, const std::pair<double, double> &r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::ZPK &zpk,
                                 IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an IIR bandstop filter.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequencies where 1 is the
    ///                   Nyquist frequency.  Here, r.first is the low corner
    ///                   and r.second is the high corner, and it is required
    ///                   that r.second > r.first.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] ba     The bandstop filter.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designBandstopIIRFilter(int n, const std::pair<double, double> &r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::BA &ba,
                                 IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @brief Designs an IIR bandstop filter stored as second order sections.
    /// @param[in] n      The order of the filter.
    /// @param[in] r      The normalized cutoff frequencies where 1 is the
    ///                   Nyquist frequency.  Here, r.first is the low corner
    ///                   and r.second is the high corner, and it is required
    ///                   that r.second > r.first.
    /// @param[in] ftype  The filter prototype.
    /// @param[in] ripple  For Chebyshev I filters this is the maximum ripple
    ///                    in the passband specified in dB.
    ///                    For Chebyshev II filters this is the maximum ripple
    ///                    in the stopband specified in dB.
    ///                    For Butterworth and Bessel filters this is ignored.
    /// @param[out] sos   The bandpass filter design.
    /// @param[in] pairing   Defines the pole pairing policy.
    /// @param[in] ldigital  Specifies whether the filter is digital or analog.
    /// @throws std::invalid_argument if any parameters are incorrect.
    void designBandstopIIRFilter(int n, const std::pair<double, double> &r,
                                 IIRPrototype ftype,
                                 double ripple,
                                 FilterRepresentations::SOS &sos,
                                 SOSPairing pairing = SOSPairing::NEAREST,
                                 IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL);
    /// @}
private:
    class FilterDesignerImpl;
    mutable std::unique_ptr<FilterDesignerImpl> pImpl;
};
}
#endif
