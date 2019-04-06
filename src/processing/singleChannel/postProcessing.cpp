#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <exception>
#ifdef DEBUG
#include <cassert>
#endif
#include <memory>
#include <algorithm>
#include <ipps.h>
#include <ippcore.h>
#ifdef __INTEL_COMPILER
#include <pstl/execution>
#include <pstl/algorithm>
#endif
#include <vector>
#include "rtseis/private/throw.hpp"
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/postProcessing/singleChannel/waveform.hpp"
#include "rtseis/postProcessing/singleChannel/detrend.hpp"
#include "rtseis/postProcessing/singleChannel/demean.hpp"
#include "rtseis/postProcessing/singleChannel/taper.hpp"
#include "rtseis/utilities/design/enums.hpp"
#include "rtseis/utilities/design/filterDesigner.hpp"
#include "rtseis/utilities/math/convolve.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/utilities/filterRepresentations/sos.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#include "rtseis/utilities/filterImplementations/iirFilter.hpp"
#include "rtseis/utilities/filterImplementations/iiriirFilter.hpp"
#include "rtseis/utilities/filterImplementations/sosFilter.hpp"

using namespace RTSeis;
using namespace PostProcessing::SingleChannel;

static inline 
double computeNyquistFrequencyFromSamplingPeriod(const double dt);
static inline std::pair<double,double>
computeNormalizedFrequencyFromSamplingPeriod(const std::pair<double,double> fc,
                                            const double dt);
static inline
double computeNormalizedFrequencyFromSamplingPeriod(const double fc,
                                                    const double dt);
static inline Utilities::Math::Convolve::Mode
classifyConvolveMode(const ConvolutionMode mode);
static inline Utilities::Math::Convolve::Implementation
classifyConvolveImplementation(const ConvolutionImplementation implementation);
static inline Utilities::FilterDesign::IIRPrototype
classifyIIRPrototype(const IIRPrototype prototype);
static inline Utilities::FilterDesign::FIRWindow
classifyFIRWindow(const FIRWindow windowIn);

#define FIR_REMOVE_PHASE(pImpl, fir) \
{ \
    const std::vector<double> taps = fir.getFilterTaps(); \
    const int ntaps = static_cast<int> (taps.size()); \
    int nhalf = ntaps/2; \
    int npad = len + nhalf; \
    double *xtemp = ippsMalloc_64f(npad); \
    double *ytemp = ippsMalloc_64f(npad); \
    const double *x = pImpl->getInputDataPointer(); \
    ippsCopy_64f(x, xtemp, len); \
    ippsZero_64f(&xtemp[len], npad - len); \
    RTSeis::Utilities::FilterImplementations::FIRFilter firFilter; \
    firFilter.initialize(ntaps, taps.data(), \
                   ProcessingMode::POST_PROCESSING, \
                   Precision::DOUBLE, \
                   Utilities::FilterImplementations::FIRImplementation::DIRECT); \
    firFilter.apply(npad, xtemp, ytemp); \
    ippsFree(xtemp); \
    pImpl->resizeOutputData(len); \
    double *yout = pImpl->getOutputDataPointer(); \
    ippsCopy_64f(&ytemp[nhalf], yout, len); \
    ippsFree(ytemp); \
};

/*
static inline void reverse(std::vector<double> &x)
{
#ifdef __INTEL_COMPLIER
    std::reverse(pstl::execution::unseq, x.begin(), x.end()); 
#else
    int len = static_cast<int> (x.size());
    ippsFlip_64f_I(x.data(), len);
#endif
    return;
}
*/


class Waveform::WaveformImpl
{
public:
    /// Default constructor
    WaveformImpl(void)
    {
        return;
    }
/*
    /// Copy constructor
    WaveformImpl(const WaveformImpl &waveform)
    {
        *this = waveform;
        return;
    }
    /// Move constructor
    WaveformImpl(WaveformImpl &&waveform)
    {
        *this = std::move(waveform);
        return;
    }
*/
    /// Default destructor
    ~WaveformImpl(void)
    {
        clear();
        return;
    }
    /// Resets the module
    void clear(void) noexcept
    {
        filterDesigner.clear();
        if (x_){ippsFree(x_);}
        if (y_){ippsFree(y_);}
        x_ = nullptr;
        y_ = nullptr;
        xptr_ = nullptr; //.release(); // = nullptr;
        maxx_ = 0;
        nx_ = 0;
        maxy_ = 0;
        ny_ = 0;
        return;
    }
    /// Gets the number of input sapmles
    int getNumberOfInputSamples(void) const noexcept
    {
        return nx_;
    }
    /// Gets the number of output samples
    int getNumberOfOutputSamples(void) const noexcept
    {
        return ny_;
    }
    /// Sets a pointer to the input data
    void setInputDataPointer(const int nx,
                             const double *x) noexcept
    {
        xptr_ = nullptr; //.release();
        nx_ = nx;
        xptr_ = x; //std::move(x);
    }
    /// Releases the input pointer to the data
    void releaseInputDataPointer(void) noexcept
    {
        xptr_ = nullptr; //.release();
        nx_ = 0;
        return;
    }
    /// Sets the input time series
    void setData(const size_t n, const double x[]) 
    {
        xptr_ = nullptr; //.release();
        nx_ = static_cast<int> (n); 
        ny_ = 0; // Can't have processed data before new data
        if (nx_ > maxx_)
        {
            if (x_){ippsFree(x_);}
            x_ = ippsMalloc_64f(nx_); 
            maxx_ = nx_; 
        }
        if (nx_ == 0){return;} // Nothing to copy
#ifdef __INTEL_COMPILER
        std::copy(pstl::execution::unseq, x, x+n, x_); 
#else
        ippsCopy_64f(x, x_, nx_);
#endif
        return;
    }
    /// Resizes the output
    void resizeOutputData(const int ny)
    {
        ny_ = ny;
        if (ny_ > maxy_)
        {
            if (y_){ippsFree(y_);}
            y_ = ippsMalloc_64f(ny_);
            maxy_ = ny_;
        }
        return;
    }
    /// Returns a pointer to the input data
    const double *getInputDataPointer(void) const
    {
        if (xptr_)
        {
            return xptr_;
        }
        else
        {
            return x_;
        }
    }
    /// Gets a pointer to the output data
    double *getOutputDataPointer(void)
    {
        return y_;
    }

    int getLengthOfInputSignal(void) const noexcept
    {
        return nx_; //static_cast<int> (x_.size());
    }
//private:
    Utilities::FilterDesign::FilterDesigner filterDesigner;
    /// A pointer to the input data
    const double *xptr_ = nullptr;
    /// The input data
    double *x_ = nullptr;
    /// The output data
    double *y_ = nullptr;
    /// Sampling period
    double dt = 1;
    /// Max space reserved to x
    int maxx_ = 0;
    /// Max space reserved to y
    int maxy_ = 0;
    /// Number of elements in x
    int nx_ = 0;
    /// Number of elements in y
    int ny_ = 0;
};

Waveform::Waveform(void) :
    pImpl(std::make_unique<WaveformImpl>())
{
    return;
}

Waveform::~Waveform(void)
{
    return;
}

void Waveform::setData(const std::vector<double> &x)
{
    size_t n = x.size();
    if (n < 1)
    {
        RTSEIS_ERRMSG("%s", "x has zero length");
        throw std::invalid_argument("x has zero length");
    }
    setData(n, x.data());
    return;
}

void Waveform::setDataPointer(const size_t n,
                              const double *x)
{
    if (n < 1 || x == nullptr)
    {
        if (n < 1){RTSEIS_THROW_IA("%s", "x has zero length");}
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "Invalid arguments");
    }
    pImpl->setInputDataPointer(static_cast<int> (n), x);
    return; 
}

void Waveform::releaseDataPointer(void) noexcept
{
    pImpl->releaseInputDataPointer();
    return;
}

void Waveform::setData(const size_t n, const double x[])
{
    if (n < 1 || x == nullptr)
    {
        if (n < 1){RTSEIS_THROW_IA("%s", "x has zero length");}
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "Invalid arguments");
    }
    pImpl->setData(n, x);
    return;
}

void Waveform::getData(std::vector<double> &y)
{
    int ny = pImpl->getNumberOfOutputSamples();
    y.resize(ny);
    if (ny > 0)
    {
        const double *yout = pImpl->getOutputDataPointer();
        ippsCopy_64f(yout, y.data(), ny);
    }
    return;
}

void Waveform::getData(const size_t nwork, double y[]) const
{
    int leny = pImpl->getNumberOfOutputSamples();
    if (nwork < static_cast<size_t> (leny))
    {
        throw std::invalid_argument("nwork = " + std::to_string(nwork)
                                  + " must be at least = "
                                   + std::to_string(leny));
    }
    if (leny == 0){return;}
    if (y == nullptr)
    {
        throw std::invalid_argument("y is NULL");
    }
    const double *yout = pImpl->getOutputDataPointer();
    ippsCopy_64f(yout, y, leny);
    return;
}

//----------------------------------------------------------------------------//
//                                 Utilities                                  //
//----------------------------------------------------------------------------//

/// TODO delete this function
size_t Waveform::getOutputLength(void) const
{
    return pImpl->getNumberOfOutputSamples(); //pImpl->ny_;
}

void Waveform::setSamplingPeriod(const double dt)
{
    if (dt <= 0)
    {
        RTSEIS_THROW_IA("Sampling period = %lf must be positive", dt);
    }
    pImpl->dt = dt;
    return;
}

double Waveform::getSamplingPeriod(void) const noexcept
{
    return pImpl->dt;
}

double Waveform::getNyquistFrequency(void) const noexcept
{
    double fnyq = 1.0/(2.0*pImpl->dt);
    return fnyq;
} 

//----------------------------------------------------------------------------//
//                     Convolution/Correlation/AutoCorrelation                //
//----------------------------------------------------------------------------//

void Waveform::convolve(
    const std::vector<double> &s,
    const ConvolutionMode mode,
    const ConvolutionImplementation implementation)
{
    int nx = pImpl->getLengthOfInputSignal();
    int ny = static_cast<int> (s.size());
    if (nx < 1){RTSEIS_THROW_IA("%s", "No data is set on the module");}
    if (ny < 1){RTSEIS_THROW_IA("%s", "No data points in s");}
    // Classify the convolution mode
    Utilities::Math::Convolve::Mode convcorMode;
    convcorMode = classifyConvolveMode(mode); // throws
    // Classify the convolution implementation
    Utilities::Math::Convolve::Implementation convcorImpl;
    convcorImpl = classifyConvolveImplementation(implementation); // throws
    // Perform the convolution
    try
    {
        int lenc = 0; 
        lenc = Utilities::Math::Convolve::computeConvolutionLength(nx, ny,
                                                                   convcorMode);
        pImpl->resizeOutputData(lenc); 
        int nyout;
        const double *x = pImpl->getInputDataPointer();
        double *yout = pImpl->getOutputDataPointer();
        Utilities::Math::Convolve::convolve(nx, x,
                                            ny, s.data(),
                                            lenc, &nyout, yout,
                                            convcorMode, convcorImpl);
#ifdef DEBUG
        assert(lenc == nyout);
#endif
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to compute convolution: %s", ia.what());
        pImpl->resizeOutputData(0);
    } 
    return;
}

void Waveform::correlate(
    const std::vector<double> &s, 
    const ConvolutionMode mode,
    const ConvolutionImplementation implementation)
{
    int nx = pImpl->getLengthOfInputSignal();
    int ny = static_cast<int> (s.size());
    if (nx < 1){RTSEIS_THROW_IA("%s", "No data is set on the module");}
    if (ny < 1){RTSEIS_THROW_IA("%s", "No data points in s");}
    // Classify the convolution mode
    Utilities::Math::Convolve::Mode convcorMode;
    convcorMode = classifyConvolveMode(mode); // throws
    // Classify the convolution implementation
    Utilities::Math::Convolve::Implementation convcorImpl;
    convcorImpl = classifyConvolveImplementation(implementation); // throws
    // Perform the correlation
    try
    {
        int lenc = 0;
        lenc = Utilities::Math::Convolve::computeConvolutionLength(nx, ny,
                                                                   convcorMode);
        pImpl->resizeOutputData(lenc); 
        int nyout;
        const double *x = pImpl->getInputDataPointer();
        double *yout = pImpl->getOutputDataPointer();
        Utilities::Math::Convolve::correlate(nx, x,
                                             ny, s.data(),
                                             lenc, &nyout, yout,
                                             convcorMode, convcorImpl);
#ifdef DEBUG
        assert(lenc == nyout);
#endif
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to compute correlation: %s", ia.what());
        pImpl->resizeOutputData(0);
    }
    return;
}

void Waveform::autocorrelate(
    const ConvolutionMode mode,
    const ConvolutionImplementation implementation)
{
    int nx = pImpl->getLengthOfInputSignal();
    if (nx < 1){RTSEIS_THROW_IA("%s", "No data is set on the module");}
    // Classify the convolution mode
    Utilities::Math::Convolve::Mode convcorMode;
    convcorMode = classifyConvolveMode(mode); // throws
    // Classify the convolution implementation
    Utilities::Math::Convolve::Implementation convcorImpl;
    convcorImpl = classifyConvolveImplementation(implementation); // throws
    // Perform the correlation
    try
    {
        int lenc = 0; 
        lenc = Utilities::Math::Convolve::computeConvolutionLength(nx, nx,
                                                                   convcorMode);
        pImpl->resizeOutputData(lenc);
        int nyout;
        const double *x = pImpl->getInputDataPointer();
        double *yout = pImpl->getOutputDataPointer();
        Utilities::Math::Convolve::autocorrelate(nx, x,
                                                 lenc, &nyout, yout,
                                                 convcorMode, convcorImpl);
#ifdef DEBUG
        assert(lenc == nyout);
#endif
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to compute autocorrelation: %s", ia.what());
        pImpl->resizeOutputData(0);
    }
    return;
}

//----------------------------------------------------------------------------//
//                            Demeaning/detrending                            //
//----------------------------------------------------------------------------//

void Waveform::demean(void)
{
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Demean the data
    try
    {
        DemeanParameters parms(RTSeis::Precision::DOUBLE);
        Demean demean(parms);
        const double *x = pImpl->getInputDataPointer();
        pImpl->resizeOutputData(len);
        double  *y = pImpl->getOutputDataPointer();
        demean.apply(len, x, y);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        throw std::runtime_error("Algorithmic failure");
    }
    return; 
}

void Waveform::detrend(void)
{
    int len = pImpl->getLengthOfInputSignal();
    if (len < 2)
    {
        RTSEIS_WARNMSG("%s", "At least 2 data points required to detrend");
        return;
    }
    // Detrend the data
    DetrendParameters parms(RTSeis::Precision::DOUBLE);
    Detrend detrend(parms);
    const double *x = pImpl->getInputDataPointer();
    pImpl->resizeOutputData(len);
    double  *y = pImpl->getOutputDataPointer();
    detrend.apply(len, x, y); 
    return;
}
//----------------------------------------------------------------------------//
//                           Band-specific Filters                            //
//----------------------------------------------------------------------------//

void Waveform::iirLowpassFilter(const int order, const double fc,
                                const IIRPrototype prototype,
                                const double ripple,
                                const bool lzeroPhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {   
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    Utilities::FilterRepresentations::BA ba;
    pImpl->filterDesigner.designLowpassIIRFilter(
                        order, r, ptype, ripple, ba,
                        Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    iirFilter(ba, lzeroPhase);
    return;
}

void Waveform::sosLowpassFilter(const int order, const double fc, 
                                const IIRPrototype prototype,
                                const double ripple,
                                const bool lzeroPhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    Utilities::FilterRepresentations::SOS sos;
    pImpl->filterDesigner.designLowpassIIRFilter(
                     order, r, ptype, ripple, sos,
                     Utilities::FilterDesign::SOSPairing::NEAREST,
                     Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    sosFilter(sos, lzeroPhase);
    return;
}

void Waveform::firLowpassFilter(const int ntapsIn, const double fc,
                                const FIRWindow windowIn,
                                const bool lremovePhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    int ntaps = ntapsIn;
    if (lremovePhase && ntaps%2 == 0)
    {
        RTSEIS_WARNMSG("%s", "Adding a filter tap");
        ntaps = ntaps + 1;
    }
    if (ntaps < 5){RTSEIS_THROW_IA("ntaps = %d  must be at least 5", ntaps);}
    int order = ntaps - 1;
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::FIRWindow window;
    window = classifyFIRWindow(windowIn);
    Utilities::FilterRepresentations::FIR fir;
    pImpl->filterDesigner.designLowpassFIRFilter(order, r, window, fir);
    if (!lremovePhase)
    {
        firFilter(fir);
    }
    else
    {
        FIR_REMOVE_PHASE(pImpl, fir);
    }
    return;
}

void Waveform::iirHighpassFilter(const int order, const double fc, 
                                 const IIRPrototype prototype,
                                 const double ripple,
                                 const bool lzeroPhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {   
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }   
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    Utilities::FilterRepresentations::BA ba; 
    pImpl->filterDesigner.designHighpassIIRFilter(
                     order, r, ptype, ripple, ba, 
                     Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    iirFilter(ba, lzeroPhase);
    return;
}

void Waveform::sosHighpassFilter(const int order, const double fc, 
                                 const IIRPrototype prototype,
                                 const double ripple,
                                 const bool lzeroPhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {   
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }   
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    Utilities::FilterRepresentations::SOS sos;
    pImpl->filterDesigner.designHighpassIIRFilter(
                    order, r, ptype, ripple, sos,
                    Utilities::FilterDesign::SOSPairing::NEAREST,
                    Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    sosFilter(sos, lzeroPhase);
    return;
}

void Waveform::firHighpassFilter(const int ntapsIn, const double fc, 
                                 const FIRWindow windowIn,
                                 const bool lremovePhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    int ntaps = ntapsIn;
    if (lremovePhase && ntaps%2 == 0)
    {
        RTSEIS_WARNMSG("%s", "Adding a filter tap");
        ntaps = ntaps + 1;
    }
    if (ntaps < 5){RTSEIS_THROW_IA("ntaps = %d  must be at least 5", ntaps);}
    int order = ntaps - 1;
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::FIRWindow window;
    window = classifyFIRWindow(windowIn);
    Utilities::FilterRepresentations::FIR fir;
    pImpl->filterDesigner.designHighpassFIRFilter(order, r, window, fir);
    // Standard FIR filtering
    if (!lremovePhase)
    {
        firFilter(fir);
    }
    else
    {
        FIR_REMOVE_PHASE(pImpl, fir);
    }
    return;
}


void Waveform::iirBandpassFilter(const int order,
                                 const std::pair<double,double> fc, 
                                 const IIRPrototype prototype,
                                 const double ripple,
                                 const bool lzeroPhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {   
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    Utilities::FilterRepresentations::BA ba; 
    pImpl->filterDesigner.designBandpassIIRFilter(
                    order, r, ptype, ripple, ba, 
                     Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    iirFilter(ba, lzeroPhase);
    return;
}

void Waveform::sosBandpassFilter(const int order,
                                 const std::pair<double,double> fc, 
                                 const IIRPrototype prototype,
                                 const double ripple,
                                 const bool lzeroPhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {   
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }   
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    Utilities::FilterRepresentations::SOS sos;
    pImpl->filterDesigner.designBandpassIIRFilter(
                    order, r, ptype, ripple, sos,
                    Utilities::FilterDesign::SOSPairing::NEAREST,
                    Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    sosFilter(sos, lzeroPhase);
    return;
}

void Waveform::firBandpassFilter(const int ntapsIn,
                                 const std::pair<double,double> fc, 
                                 const FIRWindow windowIn,
                                 const bool lremovePhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    int ntaps = ntapsIn;
    if (lremovePhase && ntaps%2 == 0)
    {
        RTSEIS_WARNMSG("%s", "Adding a filter tap");
        ntaps = ntaps + 1;
    }
    if (ntaps < 5){RTSEIS_THROW_IA("ntaps = %d  must be at least 5", ntaps);}
    int order = ntaps - 1;
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::FIRWindow window;
    window = classifyFIRWindow(windowIn);
    Utilities::FilterRepresentations::FIR fir;
    pImpl->filterDesigner.designBandpassFIRFilter(order, r, window, fir);
    if (!lremovePhase)
    {
        firFilter(fir);
    }
    else
    {
        FIR_REMOVE_PHASE(pImpl, fir);
    }
    return;
}

void Waveform::iirBandstopFilter(const int order,
                                 const std::pair<double,double> fc, 
                                 const IIRPrototype prototype,
                                 const double ripple,
                                 const bool lzeroPhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    Utilities::FilterRepresentations::BA ba;
    pImpl->filterDesigner.designBandstopIIRFilter(
                    order, r, ptype, ripple, ba,
                    Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    iirFilter(ba, lzeroPhase);
    return;
}

void Waveform::sosBandstopFilter(const int order,
                                 const std::pair<double,double> fc,
                                 const IIRPrototype prototype,
                                 const double ripple,
                                 const bool lzeroPhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    Utilities::FilterRepresentations::SOS sos;
    pImpl->filterDesigner.designBandstopIIRFilter(
                    order, r, ptype, ripple, sos,
                    Utilities::FilterDesign::SOSPairing::NEAREST,
                    Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    sosFilter(sos, lzeroPhase);
    return;
}

void Waveform::firBandstopFilter(const int ntapsIn,
                                 const std::pair<double,double> fc,
                                 const FIRWindow windowIn,
                                 const bool lremovePhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    int ntaps = ntapsIn;
    if (lremovePhase && ntaps%2 == 0)
    {
        RTSEIS_WARNMSG("%s", "Adding a filter tap");
        ntaps = ntaps + 1;
    }
    if (ntaps < 5){RTSEIS_THROW_IA("ntaps = %d  must be at least 5", ntaps);}
    int order = ntaps - 1;
    // Compute normalized frequencies
    std::pair<double,double> r
         = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt);
    Utilities::FilterDesign::FIRWindow window;
    window = classifyFIRWindow(windowIn);
    Utilities::FilterRepresentations::FIR fir;
    pImpl->filterDesigner.designBandstopFIRFilter(order, r, window, fir);
    if (!lremovePhase)
    {
        firFilter(fir);
    }
    else
    {
        FIR_REMOVE_PHASE(pImpl, fir);
    }
    return;
}

//----------------------------------------------------------------------------//
//                               General Filtering                            //
//----------------------------------------------------------------------------//

void Waveform::firFilter(const Utilities::FilterRepresentations::FIR &fir,
                         const bool lremovePhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Initialize the FIR filter
    const std::vector<double> taps = fir.getFilterTaps();
    const int nb = static_cast<int> (taps.size());
    if (nb < 1)
    {
        RTSEIS_THROW_IA("%s", "No filter taps");
        return;
    }
    // Initialize filter
    RTSeis::Utilities::FilterImplementations::FIRFilter firFilter;
    firFilter.initialize(nb, taps.data(),
                   ProcessingMode::POST_PROCESSING,
                   Precision::DOUBLE,
                   Utilities::FilterImplementations::FIRImplementation::DIRECT);
    pImpl->resizeOutputData(len);
    // Standard FIR filtering 
    const double *x = pImpl->getInputDataPointer();
    double *yout = pImpl->getOutputDataPointer();
    if (!lremovePhase)
    {
        firFilter.apply(len, x, yout);
    }
    else
    {
        double *ywork = ippsMalloc_64f(len);
        firFilter.apply(len, x,    ywork); // Filter forwards
        ippsFlip_64f(ywork, yout,  len);   // Reverse y
        firFilter.apply(len, yout, ywork); // Filter y backwards
        ippsFlip_64f(ywork, yout,  len);   // Reverse it
        ippsFree(ywork);
    }
    return;
}

void Waveform::iirFilter(const Utilities::FilterRepresentations::BA &ba,
                         const bool lremovePhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Initialize the IIR filter
    const std::vector<double> b = ba.getNumeratorCoefficients();
    const std::vector<double> a = ba.getDenominatorCoefficients();
    const int nb = static_cast<int> (b.size());
    const int na = static_cast<int> (a.size());
    if (nb < 1 || na < 1)
    {
        if (na < 1){RTSEIS_THROW_IA("%s", "No denominator coefficients");}
        if (nb < 1){RTSEIS_THROW_IA("%s", "No numerator coefficients");}
        RTSEIS_THROW_IA("%s", "No filter coefficients");
        return;
    }
    // Initialize filter
    if (!lremovePhase)
    {
        RTSeis::Utilities::FilterImplementations::IIRFilter iirFilter;
        iirFilter.initialize(nb, b.data(),
               na, a.data(),
               ProcessingMode::POST_PROCESSING,
               Precision::DOUBLE,
               Utilities::FilterImplementations::IIRDFImplementation::DF2_FAST);
        pImpl->resizeOutputData(len);
        const double *x = pImpl->getInputDataPointer();
        double *yout = pImpl->getOutputDataPointer();
        iirFilter.apply(len, x, yout);
    }
    else
    {
        RTSeis::Utilities::FilterImplementations::IIRIIRFilter iiriirFilter;
        iiriirFilter.initialize(nb, b.data(),
                                na, a.data(),
                                Precision::DOUBLE);
        pImpl->resizeOutputData(len);
        const double *x = pImpl->getInputDataPointer();
        double *yout = pImpl->getOutputDataPointer();
        iiriirFilter.apply(len, x, yout);
    }
    return;
}

void Waveform::sosFilter(const Utilities::FilterRepresentations::SOS &sos,
                         const bool lremovePhase)
{
    // Check that there's data
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Initialize the FIR filter
    const int ns = sos.getNumberOfSections();
    if (ns < 1)
    {
        RTSEIS_THROW_IA("%s", "No sections in fitler");
    }
    const std::vector<double> bs = sos.getNumeratorCoefficients();
    const std::vector<double> as = sos.getDenominatorCoefficients();
    // Initialize filter
    RTSeis::Utilities::FilterImplementations::SOSFilter sosFilter;
    sosFilter.initialize(ns, bs.data(), as.data(),
                         ProcessingMode::POST_PROCESSING,
                         Precision::DOUBLE);
    pImpl->resizeOutputData(len);
    // Get handles on pointers
    const double *x = pImpl->getInputDataPointer();
    double *yout = pImpl->getOutputDataPointer();
    // Zero-phase filtering needs workspace so that x isn't annihalated
    if (lremovePhase)
    {
        double *ywork = ippsMalloc_64f(len);
        sosFilter.apply(len, x,    ywork); // Filter forwards
        ippsFlip_64f(ywork, yout,  len);   // Reverse y
        sosFilter.apply(len, yout, ywork); // Filter y backwards
        ippsFlip_64f(ywork, yout,  len);   // Reverse it
        ippsFree(ywork);
    }
    else
    {
        sosFilter.apply(len, x, yout);
    }
    return;
}

//----------------------------------------------------------------------------//
//                                   Tapering                                 //
//----------------------------------------------------------------------------//

void Waveform::taper(const double pct,
                     const TaperParameters::Type window)
{
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Taper the data
    TaperParameters parms(pct, window, RTSeis::Precision::DOUBLE);
    Taper taper(parms);
    const double *x = pImpl->getInputDataPointer();
    pImpl->resizeOutputData(len);
    double  *y = pImpl->getOutputDataPointer();
    taper.apply(len, x, y);
    return;
}

//----------------------------------------------------------------------------//
//                              Private Functions                             //
//----------------------------------------------------------------------------//
std::pair<double,double>
computeNormalizedFrequencyFromSamplingPeriod(const std::pair<double,double> fc,
                                            const double dt)
{
    std::pair<double,double> r(0, 0);
    double fnyq = computeNyquistFrequencyFromSamplingPeriod(dt);
    if (fc.first < 0)
    {
        RTSEIS_THROW_IA("fc.first=%lf must be positive", fc.first);
    }
    if (fc.first >= fc.second) 
    {
        RTSEIS_THROW_IA("fc.first=%lf must be less than fc.second=%lf",
                        fc.first, fc.second);
    }
    if (fc.second > fnyq)
    {
        RTSEIS_THROW_IA("fc.seoncd=%lf must be in range [%lf,%lf]",
                        fc.second, fc.first, fnyq);
    }
    r.first  = fc.first/fnyq;
    r.second = fc.second/fnyq;
    return r;
}

double computeNormalizedFrequencyFromSamplingPeriod(const double fc,
                                                    const double dt)
{
    double r = 0;
    double fnyq = computeNyquistFrequencyFromSamplingPeriod(dt); 
    if (fc < 0 || fc > fnyq)
    {
        RTSEIS_THROW_IA("fc=%lf must be in range [0,%lf]", fc, fnyq);
    }
    r = fc/fnyq;
    return r;
}
double computeNyquistFrequencyFromSamplingPeriod(const double dt)
{
#ifdef DEBUG
    assert(dt > 0);
#endif
    return 1.0/(2.0*dt);
}
Utilities::Math::Convolve::Mode
classifyConvolveMode(const ConvolutionMode mode)
{
    Utilities::Math::Convolve::Mode convcorMode;
    if (mode == ConvolutionMode::FULL)
    {
        convcorMode = Utilities::Math::Convolve::Mode::FULL;
    }
    else if (mode == ConvolutionMode::SAME)
    {
        convcorMode = Utilities::Math::Convolve::Mode::SAME;
    }
    else if (mode == ConvolutionMode::VALID)
    {
        convcorMode = Utilities::Math::Convolve::Mode::VALID;
    }
    else
    {
        RTSEIS_THROW_IA("Unsupported convolution mode=%d",
                        static_cast<int> (mode));
    }
    return convcorMode;
}

Utilities::Math::Convolve::Implementation 
classifyConvolveImplementation(const ConvolutionImplementation implementation)
{
    // Classify the convolution implementaiton
    Utilities::Math::Convolve::Implementation convcorImpl;
    if (implementation == ConvolutionImplementation::AUTO)
    {   
        convcorImpl = Utilities::Math::Convolve::Implementation::AUTO;
    }   
    else if (implementation == ConvolutionImplementation::DIRECT)
    {   
        convcorImpl = Utilities::Math::Convolve::Implementation::DIRECT;
    }   
    else if (implementation == ConvolutionImplementation::FFT)
    {   
        convcorImpl = Utilities::Math::Convolve::Implementation::FFT;
    }   
    else
    {   
        RTSEIS_THROW_IA("Unsupported convolution implementation=%d",
                        static_cast<int> (implementation));
    }
    return convcorImpl;
}

Utilities::FilterDesign::IIRPrototype
classifyIIRPrototype(const IIRPrototype prototype)
{
    Utilities::FilterDesign::IIRPrototype ptype;
    if (prototype == IIRPrototype::BESSEL)
    {
        ptype = Utilities::FilterDesign::IIRPrototype::BESSEL;
    }
    else if (prototype == IIRPrototype::BUTTERWORTH)
    {
        ptype = Utilities::FilterDesign::IIRPrototype::BUTTERWORTH;
    }
    else if (prototype == IIRPrototype::CHEBYSHEV1)
    {
        ptype = Utilities::FilterDesign::IIRPrototype::CHEBYSHEV1;
    }
    else if (prototype == IIRPrototype::CHEBYSHEV2)
    {
        ptype = Utilities::FilterDesign::IIRPrototype::CHEBYSHEV2;
    }
    else
    {
        RTSEIS_THROW_IA("Unsupported prototype=%d",
                        static_cast<int> (prototype));
    }
    return ptype;
}

Utilities::FilterDesign::FIRWindow
classifyFIRWindow(const FIRWindow windowIn)
{
    Utilities::FilterDesign::FIRWindow window;
    if (windowIn == FIRWindow::HAMMING)
    {
        window = Utilities::FilterDesign::FIRWindow::HAMMING;
    }
    else if (windowIn == FIRWindow::HANN)
    {
        window = Utilities::FilterDesign::FIRWindow::HANN;
    }
    else if (windowIn == FIRWindow::BLACKMAN_OPT)
    {
        window = Utilities::FilterDesign::FIRWindow::BLACKMAN_OPT;
    }
    else if (windowIn == FIRWindow::BARTLETT)
    {
        window = Utilities::FilterDesign::FIRWindow::BARTLETT;
    }
    else
    {
        RTSEIS_THROW_IA("Unsupported window=%d",
                        static_cast<int> (windowIn));
    }
    return window;
}

