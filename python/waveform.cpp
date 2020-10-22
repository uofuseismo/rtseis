#include <vector>
#include <memory>
#include <exception>
#include "wrap.hpp"
#include "modules.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterRepresentations/sos.hpp"

namespace py = pybind11;
using namespace PBPostProcessing;

static RTSeis::PostProcessing::SingleChannel::IIRPrototype
stringToAnalogPrototype(const std::string &prototype);

Waveform::Waveform(void) :
    waveform_(new RTSeis::PostProcessing::SingleChannel::Waveform())
{
   return;
}

Waveform::~Waveform(void){return;}

/// Convolve
void Waveform::convolve(
    py::array_t<double, py::array::c_style | py::array::forcecast> &s,
    const std::string &smode)
{
    std::vector<double> svec(s.size());
    std::memcpy(svec.data(), s.data(), s.size()*sizeof(double));
    RTSeis::PostProcessing::SingleChannel::ConvolutionMode mode
         = RTSeis::PostProcessing::SingleChannel::ConvolutionMode::FULL;
    if (smode == "full")
    {
        mode = RTSeis::PostProcessing::SingleChannel::ConvolutionMode::FULL;
    }
    else if (smode == "valid")
    {
        mode = RTSeis::PostProcessing::SingleChannel::ConvolutionMode::VALID;
    }
    else if (smode == "same")
    {
        mode = RTSeis::PostProcessing::SingleChannel::ConvolutionMode::SAME;
    }
    else
    {
        throw std::invalid_argument("Invalid mode");
    }
    // Perform convolution
    try
    {
        waveform_->convolve(svec, mode);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s\n", ia.what());
        throw std::invalid_argument("Convolution failed");
    }
}

/// Remove mean
void Waveform::demean(void)
{
    try
    {
        waveform_->demean();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s\n", ia.what());
        throw std::invalid_argument("Demean failed");
    }
    return;
}

/// Remove trend
void Waveform::detrend(void)
{
    try
    {
        waveform_->detrend();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s\n", ia.what());
        throw std::invalid_argument("Detrend failed");
    }
    return;
}

/// Downsample
void Waveform::downsample(const int nq)
{
    if (nq < 1)
    {
        throw std::invalid_argument("nq = " + std::to_string(nq)
                                   + "must be positive");
    }
    try
    {
        waveform_->downsample(nq);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        throw std::invalid_argument("Downsample failed");
    }
    catch (const std::runtime_error &rt)
    {
        fprintf(stderr, "%s", rt.what());
        throw std::invalid_argument("Internal in downsample error");
    }
    return;
}

void Waveform::firEnvelope(const int nfir)
{
    if (nfir < 1)
    {
        throw std::invalid_argument("Number of FIR coeffs must be positive");
    }
    waveform_->firEnvelope(nfir);
}

void Waveform::envelope()
{
    waveform_->envelope();
}

/// FIR filter
void Waveform::firFilter(py::array_t<double, py::array::c_style | py::array::forcecast> &taps)
{
    if (taps.size() < 1)
    {
        throw std::invalid_argument("No filter coefficients!");
    }
    std::vector<double> tapsVec(taps.size());
    std::memcpy(tapsVec.data(), taps.data(), taps.size()*sizeof(double));
    try
    {
        RTSeis::Utilities::FilterRepresentations::FIR fir(tapsVec);
        waveform_->firFilter(fir); 
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        throw std::invalid_argument("FIR filtering failed"); 
    }
    return;
}

void Waveform::sosLowpassFilter(const double fc, const int order,
                                const std::string &prototype,
                                const double ripple,
                                const bool zeroPhase)
{
    RTSeis::PostProcessing::SingleChannel::IIRPrototype ptype; 
    ptype = stringToAnalogPrototype(prototype);
    try
    {
        waveform_->sosLowpassFilter(order, fc, ptype, ripple, zeroPhase);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        throw std::invalid_argument("SOS filtering failed");
    } 
    return;    
}

void Waveform::sosHighpassFilter(const double fc, const int order,
                                 const std::string &prototype,
                                 const double ripple,
                                 const bool zeroPhase)
{
    RTSeis::PostProcessing::SingleChannel::IIRPrototype ptype;
    ptype = stringToAnalogPrototype(prototype);
    try
    {
        waveform_->sosHighpassFilter(order, fc, ptype, ripple, zeroPhase);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        throw std::invalid_argument("SOS filtering failed");
    }
    return;
}

void Waveform::sosBandpassFilter(const std::pair<double,double> fc,
                                 const int order,
                                 const std::string &prototype,
                                 const double ripple,
                                 const bool zeroPhase)
{
    RTSeis::PostProcessing::SingleChannel::IIRPrototype ptype;
    ptype = stringToAnalogPrototype(prototype);
    try
    {
        waveform_->sosBandpassFilter(order, fc, ptype, ripple, zeroPhase);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        throw std::invalid_argument("SOS filtering failed");
    }
    return;
}

void Waveform::sosBandstopFilter(const std::pair<double,double> fc, 
                                 const int order,
                                 const std::string &prototype,
                                 const double ripple,
                                 const bool zeroPhase)
{
    RTSeis::PostProcessing::SingleChannel::IIRPrototype ptype;
    ptype = stringToAnalogPrototype(prototype);
    try
    {
        waveform_->sosBandstopFilter(order, fc, ptype, ripple, zeroPhase);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        throw std::invalid_argument("SOS filtering failed");
    }
}

/// Normalization
void Waveform::normalizeMinMax(const std::pair<double, double> &targetRange)
{
    try
    {
        waveform_->normalizeMinMax(targetRange);
    }
    catch (const std::exception &e)
    {
        fprintf(stderr, "%s", e.what());
        throw std::runtime_error("minMax normalization failed");
    } 
}

void Waveform::normalizeSignBit()
{
    try
    {
        waveform_->normalizeSignBit();
    }
    catch (const std::exception &e)
    {
        fprintf(stderr, "%s", e.what());
        throw std::runtime_error("signBit normalization failed");
    }
}

void Waveform::normalizeZScore()
{
    try
    {
        waveform_->normalizeZScore();
    }
    catch (const std::exception &e)
    {
        fprintf(stderr, "%s", e.what());
        throw std::runtime_error("zScore normalization failed");
    }
}

/// Taper
void Waveform::taper(const double pct, const std::string &taperName)
{
    if (pct < 0 || pct > 100)
    {
        throw std::invalid_argument("Invalid percentage");
    }
    // Classify the wavefor type
    RTSeis::PostProcessing::SingleChannel::TaperParameters::Type window
        = RTSeis::PostProcessing::SingleChannel::TaperParameters::Type::HAMMING;
    if (taperName == "hamming")
    {
        window = RTSeis::PostProcessing::SingleChannel::TaperParameters::Type::HAMMING;
    }
    else if (taperName == "hann" || taperName == "hanning")
    {
        window = RTSeis::PostProcessing::SingleChannel::TaperParameters::Type::HANN;
    }
    else if (taperName == "blackman")
    {
        window = RTSeis::PostProcessing::SingleChannel::TaperParameters::Type::BLACKMAN;
    }
    else if (taperName == "bartlett" || taperName == "triangle")
    {
        window = RTSeis::PostProcessing::SingleChannel::TaperParameters::Type::BARTLETT;
    }
    else if (taperName == "sine")
    {
        window = RTSeis::PostProcessing::SingleChannel::TaperParameters::Type::SINE;
    }
    else
    {
        throw std::invalid_argument("Unknown taper: " + taperName);
    }
    // Taper
    try
    {
        waveform_->taper(pct, window);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s\n", ia.what());
        throw std::invalid_argument("Taper failed");
    }
    return;
}

/// Sets the data
void Waveform::setData(py::array_t<double, py::array::c_style | py::array::forcecast> &x)
{
/*
    // This is really slow but simple
    std::vector<double> xvec(x.size());
    std::memcpy(xvec.data(), x.data(), x.size()*sizeof(double)); 
    waveform_->setData(xvec);
*/
    // Use pointers
    py::buffer_info xbuf = x.request();
    size_t len = xbuf.size;
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr)
    {
        throw std::runtime_error("x is null");
    }
    //waveform_->setDataPointer(len, xptr);
    waveform_->setData(len, xptr);
    return;
}

/// Gets the data
py::array_t<double> Waveform::getData()
{
    size_t ny = waveform_->getOutputLength();
    auto y = py::array_t<double, py::array::c_style> (ny);
    py::buffer_info ybuf = y.request();
    auto yptr = static_cast<double *> (ybuf.ptr);
    waveform_->getData(ny, &yptr);
    return y;
/*
    std::vector<double> y;
    waveform_->getData(y);
    return py::array(y.size(), y.data());
*/
}

void Waveform::setSamplingPeriod(const double dt)
{
    if (dt <= 0)
    {
        throw std::invalid_argument("Sampling period = "
                                  + std::to_string(dt)
                                  + " must be positive");
    }
    waveform_->setSamplingPeriod(dt);
}

double Waveform::getSamplingPeriod() const
{
    return waveform_->getSamplingPeriod();
}

/// Interpolate a signal
void Waveform::interpolate(
    const double newSamplingPeriod,
    const RTSeis::PostProcessing::SingleChannel::InterpolationMethod method) 
{
    if (newSamplingPeriod <= 0)
    {
        throw std::invalid_argument("New sampling period = "
                                  + std::to_string(newSamplingPeriod)
                                  + " must be positive");
    }
    waveform_->interpolate(newSamplingPeriod, method);
}

/// Checks if initialized
bool Waveform::isInitialized() const
{
    return true;//waveform_->isInitialized();
}

///==========================================================================///
///                          Private Functions                               ///
///==========================================================================///

RTSeis::PostProcessing::SingleChannel::IIRPrototype
stringToAnalogPrototype(const std::string &prototype)
{
    RTSeis::PostProcessing::SingleChannel::IIRPrototype ptype;
    if (prototype == "bessel")
    {
        ptype = RTSeis::PostProcessing::SingleChannel::IIRPrototype::BESSEL;
    }
    else if (prototype == "butterworth" || prototype == "butter")
    {
        ptype = RTSeis::PostProcessing::SingleChannel::IIRPrototype::BUTTERWORTH;
    }
    else if (prototype == "chebyshev1" || prototype == "cheby1")
    {
        ptype = RTSeis::PostProcessing::SingleChannel::IIRPrototype::CHEBYSHEV1;
    }
    else if (prototype == "chebyshev2" || prototype == "cheby2")
    {
        ptype = RTSeis::PostProcessing::SingleChannel::IIRPrototype::CHEBYSHEV2;
    }
    else
    {
        throw std::invalid_argument("Unknown prototype " + prototype);
    }
    return ptype;

}

///==========================================================================///
///                          Python Bindings                                 ///
///==========================================================================///

void init_pp_waveform(py::module &m)
{
    m.doc() = "Utilities for post-processing waveforms";

    py::class_<PBPostProcessing::Waveform> singleChannelWaveform(m, "Waveform");

    singleChannelWaveform.def(py::init<>());
    singleChannelWaveform.doc() = "Single channel waveform post-processing";
    singleChannelWaveform.def("set_data", &PBPostProcessing::Waveform::setData,
                              "Sets the signal to process on the class");
    singleChannelWaveform.def("get_data", &PBPostProcessing::Waveform::getData,
                              "Gets the filtered data as a NumPy array");
    singleChannelWaveform.def("set_sampling_period", 
                              &PBPostProcessing::Waveform::setSamplingPeriod,
                              "Sets the sampling period (seconds)");
    singleChannelWaveform.def("get_sampling_period",
                              &PBPostProcessing::Waveform::getSamplingPeriod,
                              "Gets the sampling period (seconds).");
    singleChannelWaveform.def("convolve",  &PBPostProcessing::Waveform::convolve,
                              "Convolves the time series with the input signal",
                              py::arg("s"),
                              py::arg("smode") = "full");
    singleChannelWaveform.def("demean",  &PBPostProcessing::Waveform::demean,
                              "Removes the mean from the time series");
    singleChannelWaveform.def("detrend", &PBPostProcessing::Waveform::detrend,
                              "Removes the trend from the time series");
    singleChannelWaveform.def("downsample", &PBPostProcessing::Waveform::downsample,
                              "Downsamples sample signal by integer factor",
                              py::arg("nq"));
    singleChannelWaveform.def("fir_filter", &PBPostProcessing::Waveform::firFilter,
                              py::arg("taps"));

    singleChannelWaveform.def("envelope", &PBPostProcessing::Waveform::envelope,
                              "Computes the envelope using the Fourier transform");
    singleChannelWaveform.def("fir_envelope", &PBPostProcessing::Waveform::firEnvelope,
                              "Computes the envelope using a FIR-based Hilbert transformer",
                              py::arg("nfir") = 301);
    singleChannelWaveform.def("interpolate", &PBPostProcessing::Waveform::interpolate,
                              "Interpolates a signal",
                              py::arg("new_sampling_period"),
                              py::arg("method") = RTSeis::PostProcessing::SingleChannel::InterpolationMethod::DFT);

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
    singleChannelWaveform.def("normalize_min_max", &PBPostProcessing::Waveform::normalizeMinMax,
                              "Min-max normalization of a signal",
                              py::arg("target_range") = std::make_pair<double, double> (0.0, 1.0));
    singleChannelWaveform.def("normalize_sign_bit", &PBPostProcessing::Waveform::normalizeSignBit,
                              "Sign-bit normalization of a signal");
    singleChannelWaveform.def("normalize_z_score", &PBPostProcessing::Waveform::normalizeZScore,
                              "z score normalization of a signal");
    singleChannelWaveform.def("taper",   &PBPostProcessing::Waveform::taper,
                              "Tapers the ends of a signal",
                              py::arg("pct") = 5,
                              py::arg("type") = "hamming");
    singleChannelWaveform.def("is_initialized", &PBPostProcessing::Waveform::isInitialized,
                              "Checks if the class is initialized");

    // Add some enums
    pybind11::enum_<RTSeis::PostProcessing::SingleChannel::InterpolationMethod> (m, "InterpolationType")
        .value("dft", RTSeis::PostProcessing::SingleChannel::InterpolationMethod::DFT,
               "Interpolates using Fourier interpolation - i.e., zero stuffing in the frequency domain.")
        .value("weighted_average_slopes", RTSeis::PostProcessing::SingleChannel::InterpolationMethod::WEIGHTED_AVERAGE_SLOPES,
               "Interpolates using the weighted-average slopes method of Wiggins.  This is the algorithm used in SAC.");
}
