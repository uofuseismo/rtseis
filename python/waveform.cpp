#include <vector>
#include <memory>
#include <exception>
#include "wrap.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include "include/rtseis/utilities/filterRepresentations/ba.hpp"
#include "include/rtseis/utilities/filterRepresentations/fir.hpp"
#include "include/rtseis/utilities/filterRepresentations/sos.hpp"

namespace py = pybind11;
using namespace PBPostProcessing;

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
        waveform_->filter(fir); 
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        throw std::invalid_argument("FIR filtering failed"); 
    }
    return;
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
    int len = xbuf.size;
    const double *xptr = (double *) (xbuf.ptr);
    if (xptr == nullptr)
    {
        throw std::runtime_error("x is null");
    }
    waveform_->setData(len, xptr);
    return;
}

/// Gets the data
py::array_t<double> Waveform::getData(void)
{
    size_t ny = waveform_->getOutputLength();
    auto y = py::array_t<double, py::array::c_style> (ny);
    py::buffer_info ybuf = y.request();
    double *yptr = (double *) (ybuf.ptr);
    waveform_->getData(ny, yptr);
    return y;
/*
    std::vector<double> y;
    waveform_->getData(y);
    return py::array(y.size(), y.data());
*/
}
/// Checks if initialized
bool Waveform::isInitialized(void) const
{
    return true;//waveform_->isInitialized();
}
