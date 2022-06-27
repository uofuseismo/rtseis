#include <iostream>
#include <vector>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include "rtseis/transforms/continuousWavelet.hpp"
#include "rtseis/transforms/wavelets/morlet.hpp"
#include "modules.hpp"

using namespace PTransforms;
namespace RTransforms = RTSeis::Transforms;

namespace
{
class Morlet
{
public:
    Morlet() :
        mWavelet(std::make_unique<RTransforms::Wavelets::Morlet> ())
    {
    }
    Morlet(const Morlet &m)
    {
        *this = m;
    }
    Morlet(Morlet &&m) noexcept
    {
        *this = std::move(m);
    }
    ~Morlet() = default;
    Morlet& operator=(const Morlet &m)
    {
        if (&m == this){return *this;}
        mWavelet = std::make_unique<RTransforms::Wavelets::Morlet> (*m.mWavelet);
        return *this;
    }
    Morlet& operator=(Morlet &&m) noexcept
    {
        if (&m == this){return *this;}
        mWavelet = std::move(m.mWavelet);
        return *this;
    }
    void enableNormalization() noexcept
    {
        mWavelet->enableNormalization();
    }
    void disableNormalization() noexcept
    {
        mWavelet->disableNormalization();
    }
    [[nodiscard]] bool normalize() const noexcept
    {
        return mWavelet->normalize();
    }
    void setParameter(const double omega0)
    {
        mWavelet->setParameter(omega0);
    }
    [[nodiscard]] double getParameter() const noexcept
    {
        return mWavelet->getParameter();
    }
    [[nodiscard]] std::vector<std::complex<double>>
        evaluate(const int n, const double scale)
    {
        std::vector<std::complex<double>> wavelet(std::max(n, 0));
        if (n < 1){return wavelet;}
        auto wPtr = wavelet.data();
        mWavelet->evaluate(n, scale, &wPtr);
        return wavelet;
    }
    [[nodiscard]] RTransforms::Wavelets::Morlet getNativeClass() const
    {
        return *mWavelet;
    }
    std::unique_ptr<RTransforms::Wavelets::Morlet> mWavelet; 
};
///--------------------------------------------------------------------------///
class ContinuousWaveletTransform
{
public:
    /// Constructor
    ContinuousWaveletTransform() :
        mCWT(std::make_unique<RTransforms::ContinuousWavelet<double>> ())
    {
    }
    /// Destructor
    ~ContinuousWaveletTransform() = default;
    /// Initialize
    void initialize(int nSamples,
                    const std::vector<double> &scales,
                    const ::Morlet &morlet, 
                    double samplingRate = 100)
    {
        auto nScales = static_cast<int> (scales.size());
        auto wavelet = morlet.getNativeClass();
        mCWT->initialize(nSamples, nScales, scales.data(),
                         wavelet, samplingRate);
    }
    /// Initialized?
    [[nodiscard]] bool isInitialized() const noexcept
    {
        return mCWT->isInitialized();
    }
    [[nodiscard]] int getNumberOfScales() const
    {
        return mCWT->getNumberOfScales();
    }
    [[nodiscard]] int getNumberOfSamples() const
    {
        return mCWT->getNumberOfSamples();
    }
    /// Transform
    void transform(const std::vector<double> &x)
    {
         mCWT->transform(x.size(), x.data());
    }
    /// Have transform?
    [[nodiscard]] bool haveTransform() const noexcept
    {
         return mCWT->haveTransform();
    }
    std::vector<double> getAmplitude() const
    {
         auto nSamples = getNumberOfSamples();
         auto nScales = getNumberOfScales();
         auto nWork = nSamples*nScales;
         std::vector<double> amplitude(nWork, -100);
         auto ptr = amplitude.data(); 
         mCWT->getAmplitudeTransform(nSamples, nScales, &ptr);
         return amplitude;
    } 
    std::vector<double> getPhase() const
    {   
         auto nSamples = getNumberOfSamples();
         auto nScales = getNumberOfScales();
         auto nWork = nSamples*nScales;
         std::vector<double> phase(nWork, 0);
         auto ptr = phase.data(); 
         mCWT->getPhaseTransform(nSamples, nScales, &ptr);
         return phase;
    }                                     
    std::vector<std::complex<double>> getTransform() const
    {
         auto nSamples = getNumberOfSamples();
         auto nScales = getNumberOfScales();
         auto nWork = nSamples*nScales;
         std::vector<std::complex<double>> tform(nWork);
         auto ptr = tform.data(); 
         mCWT->getTransform(nSamples, nScales, &ptr);
         return tform;
    }
    /// Gets the amplitude transform
    std::unique_ptr<RTransforms::ContinuousWavelet<double>> mCWT;

    ContinuousWaveletTransform(const ContinuousWaveletTransform &) = delete;
    ContinuousWaveletTransform(ContinuousWaveletTransform &&) noexcept = delete;
    ContinuousWaveletTransform& operator=(const ContinuousWaveletTransform &) = delete;
    ContinuousWaveletTransform& operator=(ContinuousWaveletTransform &&) noexcept = delete;
}; 
}

void PTransforms::initialize(pybind11::module &m)
{
    pybind11::class_<::Morlet> morlet(m, "MorletWavelet");
    morlet.def(pybind11::init<> ());
    morlet.doc() = R""""(
Defines the Morlet wavelet.

Options :
    wavelet_parameter : The wavelet parameter.  This represents a tradeoff
                        between time and frequency resolution.  The default is 6. 
)"""";
    morlet.def("enable_normalization",  &::Morlet::enableNormalization, "This will normalize by the 1/sqrt(s) where s is the scale.");
    morlet.def("disable_normalization", &::Morlet::disableNormalization, "This will not normalize.");
    morlet.def("normalize", &::Morlet::normalize, "True indicates the wavelet will be normalized.");
    morlet.def("evaluate", &::Morlet::evaluate, "Evaluates the Morlet wavelet at the given given number of scales and scale.");
    morlet.def_property("wavelet_parameter",
                        &::Morlet::getParameter,
                        &::Morlet::setParameter);

    pybind11::class_<::ContinuousWaveletTransform> cwt(m, "ContinuousWaveletTransform");
    cwt.def(pybind11::init<> ());
    cwt.doc() = R""""(
Computes the continuous wavelet transform of a signal.

)"""";
    cwt.def("initialize", &::ContinuousWaveletTransform::initialize,
            "Initializes the CWT calculator for a signal size of nSamples, the scales which are proportional to (2*pi*omega0*samplingRate/frequency) for a Morlet wavelet and omega0 = morlet.wavelet_parameter, the (Morlet) wavelet parameters, and the signal sampling rate in Hz."); 
    cwt.def("transform", &::ContinuousWaveletTransform::transform,
            "Transforms a signal.  The class must be initialized and the signal length must match cwt.get_number_of_samples.");
    cwt.def("get_transform", &::ContinuousWaveletTransform::getTransform,
            "After calling transform this will return the complex-valued transform.  This can be reshaped by r.reshape((n_scales, n_samples), order='C')"); 
    cwt.def("get_amplitude_transform", &::ContinuousWaveletTransform::getAmplitude,
            "After calling transform this will return the amplitude transform.");
    cwt.def("get_phase_transform", &::ContinuousWaveletTransform::getPhase,
            "After calling transform this will return the phase transform in radians.");
    cwt.def_property_readonly("initialized",
                              &::ContinuousWaveletTransform::isInitialized); 
    cwt.def_property_readonly("number_of_samples",
                              &::ContinuousWaveletTransform::getNumberOfSamples);
    cwt.def_property_readonly("number_of_scales",
                              &::ContinuousWaveletTransform::getNumberOfScales);
    cwt.def_property_readonly("have_transform",
                              &::ContinuousWaveletTransform::haveTransform);

}
