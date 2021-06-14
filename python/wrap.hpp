#ifndef PYRTSEIS_WRAP_HPP
#define PYRTSEIS_WRAP_HPP 1
#include <memory>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "rtseis/postProcessing/singleChannel/waveform.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include "rtseis/filterRepresentations/sos.hpp"

namespace PBPostProcessing
{
class Waveform
{
public:
    /// Constructor
    Waveform();
    /// Destructor
    ~Waveform();
    /// Convolution
    void convolve(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &s,
                  const std::string &smode);
    /// Demean
    void demean();
    /// Detrend
    void detrend();
    /// Downsampler
    void downsample(const int nq);
    /// Envelope
    void firEnvelope(const int nfir);
    void envelope();
    /// Generic FIR filter data
    void firFilter(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &taps);
    /// Generic IIR filter data using biquad second-order sections
    /// IIR lowpass filter using second order sections
    void sosLowpassFilter(double fc, int order,
                          RTSeis::PostProcessing::SingleChannel::IIRPrototype prototype,
                          double ripple,
                          bool zeroPhase);
    /// IIR highpass filter using second order sections
    void sosHighpassFilter(double fc, int order,
                           RTSeis::PostProcessing::SingleChannel::IIRPrototype prototype,
                           double ripple,
                           bool zeroPhase); 
    /// IIR bandpass filter using second order sections
    void sosBandpassFilter(const std::pair<double,double> fc, int order,
                           RTSeis::PostProcessing::SingleChannel::IIRPrototype prototype,
                           double ripple,
                           bool zeroPhase);
    /// IIR bandstop filter using second order sections
    void sosBandstopFilter(const std::pair<double,double> fc, int order,
                           RTSeis::PostProcessing::SingleChannel::IIRPrototype prototype,
                           double ripple,
                           bool zeroPhase);
    /// Normalization
    void normalizeMinMax(const std::pair<double, double> &targetRange);
    void normalizeSignBit();
    void normalizeZScore();
    /// Taper
    void taper(double pct, const std::string &taperName); 
    /// Set/get sampling period
    void setSamplingPeriod(double dt);
    double getSamplingPeriod() const; 
    /// Set data
    void setData(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &x);
    /// Interpolate
    void interpolate(double newSamplingPeriod,
                     RTSeis::PostProcessing::SingleChannel::InterpolationMethod method);
    /// Get data
    pybind11::array_t<double> getData();
    /// Checks if class is initialized
    bool isInitialized() const;
private:
    std::unique_ptr<RTSeis::PostProcessing::SingleChannel::Waveform<double>> waveform_;
};
};

namespace PBFilterRepresentations
{
class FIR
{
public:
    FIR();
    void setTaps(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> taps);
    pybind11::array_t<double> getTaps() const;
private:
    std::unique_ptr<RTSeis::FilterRepresentations::FIR> fir_;
};

class SOS
{
public:
    SOS();
    int getNumberOfSections() const;
    pybind11::array_t<double> getNumeratorCoefficients() const;
    pybind11::array_t<double> getDenominatorCoefficients() const;
private:
    std::unique_ptr<RTSeis::FilterRepresentations::SOS> sos_;
};

class BA 
{
public:
    BA();
    void setNumeratorCoefficients(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &b);
    void setDenominatorCoefficients(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &x);
    pybind11::array_t<double> getNumeratorCoefficients() const;
    pybind11::array_t<double> getDenominatorCoefficients() const;
private:
    std::unique_ptr<RTSeis::FilterRepresentations::BA> ba_;
};

}; /// End representations

#endif
