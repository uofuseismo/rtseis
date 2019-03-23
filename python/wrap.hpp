#ifndef RTSEIS_WRAP_HPP
#define RTSEIS_WRAP_HPP 1
#include <memory>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "include/rtseis/postProcessing/singleChannel/waveform.hpp"
#include "include/rtseis/utilities/filterRepresentations/ba.hpp"
#include "include/rtseis/utilities/filterRepresentations/fir.hpp"
#include "include/rtseis/utilities/filterRepresentations/sos.hpp"

namespace PBPostProcessing
{
class Waveform
{
public:
    /// Constructor
    Waveform(void);
    /// Destructor
    ~Waveform(void);
    /// Convolution
    void convolve(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &s,
                  const std::string &smode);
    /// Demean
    void demean(void);
    /// Detrend
    void detrend(void);
    /// Generic FIR filter data
    void firFilter(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &taps);
    /// IIR lowpass filter using second-order-sections
    void sosLowpassFilter(const double fc, const int order,
                          const std::string &prototype,
                          const double ripple,
                          const bool zeroPhase);
    /// Taper
    void taper(const double pct, const std::string &taperName); 
    /// Set data
    void setData(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &x);
    /// Get data
    pybind11::array_t<double> getData(void);
    /// Checks if class is initialized
    bool isInitialized(void) const;
private:
    std::unique_ptr<RTSeis::PostProcessing::SingleChannel::Waveform> waveform_;
};
};

namespace PBFilterRepresentations
{
class FIR
{
public:
    FIR(void);
    void setTaps(pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> taps);
    pybind11::array_t<double> getTaps(void);
private:
    RTSeis::Utilities::FilterRepresentations::FIR fir;
};

};

#endif
