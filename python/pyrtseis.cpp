#include "wrap.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

PYBIND11_MODULE(libpyrtseis, modules)
{
    //------------------------------------------------------------------------//
    //                         PostProcessing Group                           //
    //------------------------------------------------------------------------// 
    py::module m = modules.def_submodule("PostProcessing");
    m.doc() = "Utilities for post-processing waveforms";

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
    singleChannelWaveform.def("taper",   &PBPostProcessing::Waveform::taper,
                              "Tapers the ends of a signal",
                              py::arg("pct") = 5,
                              py::arg("type") = "hamming");
    singleChannelWaveform.def("is_initialized", &PBPostProcessing::Waveform::isInitialized,
                              "Checks if the class is initialized");

    //m.doc() = "pybind11 example plugin"; // optional module docstring

    //m.def("add", &add, "A function which adds two numbers");

    //------------------------------------------------------------------------//
    //                       Filter Representations Group                     //
    //------------------------------------------------------------------------// 
    py::module mfr = modules.def_submodule("FilterRepresentations");
    mfr.doc() = "Utilities for representing filters";

    py::class_<PBFilterRepresentations::FIR> fir(mfr, "FIR");
    fir.def(py::init<>());
    fir.def("setTaps", &PBFilterRepresentations::FIR::setTaps,
            "Sets the filter taps");
    fir.def("getTaps", &PBFilterRepresentations::FIR::getTaps,
            "Gets the filter taps");
}
