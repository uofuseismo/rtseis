//#include "wrap.hpp"
#include "modules.hpp"
#include "rtseis/version.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

PYBIND11_MODULE(pyrtseis, modules)
{
    modules.attr("__version__") = RTSEIS_VERSION;
    modules.attr("__name__") = "pyrtseis"; 
    modules.attr("__doc__") = "The Python interface to RTSeis.";
    //------------------------------------------------------------------------//
    //                         PostProcessing Group                           //
    //------------------------------------------------------------------------// 
    py::module m = modules.def_submodule("PostProcessing");
    init_pp_waveform(m);

    py::module mTransforms = modules.def_submodule("Transforms");
    PTransforms::init_transforms(mTransforms);
/*

    py::class_<PBPostProcessing::Waveform> singleChannelWaveform(m, "Waveform");
    //singleChannelWaveform
    singleChannelWaveform.def(py::init<>());
    singleChannelWaveform.doc() = "Single channel waveform post-processing";
    singleChannelWaveform.def("set_data", &PBPostProcessing::Waveform::setData,
                              "Sets the signal to process on the class");
    singleChannelWaveform.def("get_data", &PBPostProcessing::Waveform::getData,
                              "Gets the filtered data as a NumPy array");
    singleChannelWaveform.def("convolve",  &PBPostProcessing::Waveform::convolve,
                              "Convolves the time series with the input signal",
                              py::arg("s"),
                              py::arg("smode") = "full");
    singleChannelWaveform.def("demean",  &PBPostProcessing::Waveform::demean,
                              "Removes the mean from the time series");
    singleChannelWaveform.def("detrend", &PBPostProcessing::Waveform::detrend,
                              "Removes the trend from the time series");
    singleChannelWaveform.def("fir_filter", &PBPostProcessing::Waveform::firFilter,
                              py::arg("taps"));
    singleChannelWaveform.def("sos_lowpass_filter", &PBPostProcessing::Waveform::sosLowpassFilter,
                              "Lowpass filters a signal using a biquadratic (second-order-section) filter",
                              py::arg("fc"),
                              py::arg("order") = 2,
                              py::arg("prototype") = "butterworth",
                              py::arg("ripple") = 5,
                              py::arg("zero_phase") = false);
    singleChannelWaveform.def("sos_highpass_filter", &PBPostProcessing::Waveform::sosHighpassFilter,
                              "Highpass filters a signal using a biquadratic (second-order-section) filter",
                              py::arg("fc"),
                              py::arg("order") = 2,
                              py::arg("prototype") = "butterworth",
                              py::arg("ripple") = 5,
                              py::arg("zero_phase") = false);
    singleChannelWaveform.def("sos_bandpass_filter", &PBPostProcessing::Waveform::sosBandpassFilter,
                              "Bandpass filters a signal using a biquadratic (second-order-section) filter",
                              py::arg("fc"),
                              py::arg("order") = 2,
                              py::arg("prototype") = "butterworth",
                              py::arg("ripple") = 5,
                              py::arg("zero_phase") = false);
    singleChannelWaveform.def("sos_bandstop_filter", &PBPostProcessing::Waveform::sosBandstopFilter,
                              "Lowpass filters a signal using a biquadratic (second-order-section) filter",
                              py::arg("fc"),
                              py::arg("order") = 2,
                              py::arg("prototype") = "butterworth",
                              py::arg("ripple") = 5,
                              py::arg("zero_phase") = false);
    singleChannelWaveform.def("taper",   &PBPostProcessing::Waveform::taper,
                              "Tapers the ends of a signal",
                              py::arg("pct") = 5,
                              py::arg("type") = "hamming");
    singleChannelWaveform.def("is_initialized", &PBPostProcessing::Waveform::isInitialized,
                              "Checks if the class is initialized");

    //m.doc() = "pybind11 example plugin"; // optional module docstring

    //m.def("add", &add, "A function which adds two numbers");
*/

    //------------------------------------------------------------------------//
    //                       Filter Representations Group                     //
    //------------------------------------------------------------------------// 
/*
    py::module mfr = modules.def_submodule("FilterRepresentations");
    mfr.doc() = "Utilities for representing filters";

    py::class_<PBFilterRepresentations::FIR> fir(mfr, "FIR");
    fir.def(py::init<>());
    fir.def_property("taps",
                      &PBFilterRepresentations::FIR::getTaps,
                      &PBFilterRepresentations::FIR::setTaps); 
*/
/*
    fir.def("setTaps", &PBFilterRepresentations::FIR::setTaps,
            "Sets the filter taps");
    fir.def("getTaps", &PBFilterRepresentations::FIR::getTaps,
            "Gets the filter taps");
*/
}
