#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <exception>
#ifndef NDEBUG
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
#include "private/throw.hpp"
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/postProcessing/singleChannel/waveform.hpp"
#include "rtseis/postProcessing/singleChannel/detrend.hpp"
#include "rtseis/postProcessing/singleChannel/demean.hpp"
#include "rtseis/filterDesign/enums.hpp"
#include "rtseis/filterDesign/filterDesigner.hpp"
#include "rtseis/utilities/math/convolve.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "rtseis/filterImplementations/decimate.hpp"
#include "rtseis/filterImplementations/detrend.hpp"
#include "rtseis/filterImplementations/downsample.hpp"
#include "rtseis/filterImplementations/firFilter.hpp"
#include "rtseis/filterImplementations/iirFilter.hpp"
#include "rtseis/filterImplementations/iiriirFilter.hpp"
#include "rtseis/filterImplementations/sosFilter.hpp"
#include "rtseis/filterImplementations/taper.hpp"
#include "rtseis/utilities/interpolation/interpolate.hpp"
#include "rtseis/utilities/interpolation/weightedAverageSlopes.hpp"
#include "rtseis/utilities/normalization/minMax.hpp"
#include "rtseis/utilities/normalization/signBit.hpp"
#include "rtseis/utilities/normalization/zscore.hpp"
#include "rtseis/transforms/firEnvelope.hpp"
#include "rtseis/transforms/envelope.hpp"

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
static inline FilterDesign::IIRPrototype
classifyIIRPrototype(const IIRPrototype prototype);
static inline FilterDesign::FIRWindow
classifyFIRWindow(const FIRWindow windowIn);

#define FIR_REMOVE_PHASE(pImpl, fir) \
{ \
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();} \
    int len = pImpl->getLengthOfInputSignal(); \
    if (len < 1) \
    { \
        RTSEIS_WARNMSG("%s", "No data is set on the module"); \
        return; \
    } \
    const std::vector<double> taps = fir.getFilterTaps(); \
    const int nt = static_cast<int> (taps.size()); \
    int nhalf = nt/2; \
    int npad = len + nhalf; \
    double *xtemp = ippsMalloc_64f(npad); \
    double *ytemp = ippsMalloc_64f(npad); \
    const double *x = pImpl->getInputDataPointer(); \
    ippsCopy_64f(x, xtemp, len); \
    ippsZero_64f(&xtemp[len], npad - len); \
    RTSeis::FilterImplementations::FIRFilter<RTSeis::ProcessingMode::POST, double> firFilter; \
    firFilter.initialize(nt, taps.data(), \
                   RTSeis::FilterImplementations::FIRImplementation::DIRECT); \
    firFilter.apply(npad, xtemp, &ytemp); \
    ippsFree(xtemp); \
    pImpl->resizeOutputData(len); \
    double *yout = pImpl->getOutputDataPointer(); \
    ippsCopy_64f(&ytemp[nhalf], yout, len); \
    ippsFree(ytemp); \
    pImpl->lfirstFilter_ = false; \
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


template<>
class Waveform<double>::WaveformImpl
{
public:
    /// Default constructor
    WaveformImpl() = default;
    /// Copy c'tor
    WaveformImpl(const WaveformImpl &waveform)
    {
        *this = waveform;
    }
    /// Copy assignment
    WaveformImpl& operator=(const WaveformImpl &waveform)
    {
        filterDesigner = waveform.filterDesigner;
        xptr_ = waveform.xptr_;
        dt0_ = waveform.dt0_;
        dt_ = waveform.dt_;
        maxx_ = waveform.maxx_;
        maxy_ = waveform.maxy_;
        nx_ = waveform.nx_;
        ny_ = waveform.ny_;
        lfirstFilter_ = true; 
        if (maxx_ > 0)
        {
            x_ = ippsMalloc_64f(maxx_);
            if (waveform.x_ != nullptr){ippsCopy_64f(waveform.x_, x_, maxx_);}
        }
        if (maxy_ > 0)
        {
            y_ = ippsMalloc_64f(maxy_);
            if (waveform.y_ != nullptr){ippsCopy_64f(waveform.y_, y_, maxy_);}
        }
        return *this;
    }
    /// Default destructor
    ~WaveformImpl()
    {
        clear();
    }
    /// Resets the module
    void clear() noexcept
    {
        filterDesigner.clear();
        if (x_){ippsFree(x_);}
        if (y_){ippsFree(y_);}
        x_ = nullptr;
        y_ = nullptr;
        xptr_ = nullptr; //.release(); // = nullptr;
        dt0_ = 1;
        dt_ = 1;
        maxx_ = 0;
        nx_ = 0;
        maxy_ = 0;
        ny_ = 0;
        lfirstFilter_ = true;
    }
    /// Gets the number of input sapmles
    int getNumberOfInputSamples() const noexcept
    {
        return nx_;
    }
    /// Gets the number of output samples
    int getNumberOfOutputSamples() const noexcept
    {
        return ny_;
    }
    /// Sets a pointer to the input data
    void setInputDataPointer(const int nx,
                             const double *x,
                             const bool lfirst = true) noexcept
    {
        xptr_ = nullptr; //.release();
        nx_ = nx;
        xptr_ = x; //std::move(x);
        lfirstFilter_ = lfirst;
    }
    /// Releases the input pointer to the data
    void releaseInputDataPointer() noexcept
    {
        xptr_ = nullptr; //.release();
        nx_ = 0;
    }
    /// Overwrite input data with filtered data
    void overwriteInputWithOutput()
    {
        xptr_ = nullptr; // Release
        setData(static_cast<size_t> (ny_), y_, false); // Overwrite
        ny_ = 0; // Try to avoid unnecessary allocation
    }
    /// Restores the sampling period
    void restoreSamplingPeriod() noexcept
    {
        dt_ = dt0_;
    }
    /// Sets the input time series
    void setData(const size_t n, const double x[],
                 const bool lfirst = true)
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
        lfirstFilter_ = lfirst;
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
    }
    /// Returns a pointer to the input data
    const double *getInputDataPointer() const
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
    double *getOutputDataPointer()
    {
        return y_;
    }

    int getLengthOfInputSignal() const noexcept
    {
        return nx_; //static_cast<int> (x_.size());
    }
//private:
    FilterDesign::FilterDesigner filterDesigner;
    /// A pointer to the input data
    const double *xptr_ = nullptr;
    /// The input data
    double *x_ = nullptr;
    /// The output data
    double *y_ = nullptr;
    /// Input sampling period
    double dt0_ = 1;
    /// Sampling period
    double dt_ = 1;
    /// Max space reserved to x
    int maxx_ = 0;
    /// Max space reserved to y
    int maxy_ = 0;
    /// Number of elements in x
    int nx_ = 0;
    /// Number of elements in y
    int ny_ = 0;
    /// Flag indicating this is the first filtering operation on the input data
    bool lfirstFilter_ = true; 
};

/// C'tor
template<class T>
Waveform<T>::Waveform() :
    pImpl(std::make_unique<WaveformImpl>())
{
}

/// Copy c'tor
template<class T>
Waveform<T>::Waveform(const Waveform<T> &waveform)
{
    *this = waveform;
}

/// Move c'tor
template<class T>
Waveform<T>::Waveform(Waveform<T> &&waveform) noexcept
{
    *this = std::move(waveform);
}

/// Destructor
template<class T>
Waveform<T>::~Waveform() = default;

/// Move assignment
template<class T>
Waveform<T>& Waveform<T>::operator=(Waveform &&waveform) noexcept
{
    if (&waveform == this){return *this;}
    pImpl = std::move(waveform.pImpl);
    return *this;
}

/// Copy assignment
template<class T>
Waveform<T>& Waveform<T>::operator=(const Waveform &waveform)
{
    if (&waveform == this){return *this;}
    pImpl = std::make_unique<WaveformImpl> (*waveform.pImpl);
    return *this;
}

template<class T>
void Waveform<T>::setData(const std::vector<T> &x)
{
    size_t n = x.size();
    if (n < 1)
    {
        RTSEIS_ERRMSG("%s", "x has zero length");
        throw std::invalid_argument("x has zero length");
    }
    setData(n, x.data());
}

template<class T>
void Waveform<T>::setDataPointer(const size_t n,
                                 const T *x)
{
    if (n < 1 || x == nullptr)
    {
        if (n < 1){RTSEIS_THROW_IA("%s", "x has zero length");}
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "Invalid arguments");
    }
    pImpl->restoreSamplingPeriod();
    pImpl->setInputDataPointer(static_cast<int> (n), x, true);
}

template<class T>
void Waveform<T>::releaseDataPointer() noexcept
{
    pImpl->releaseInputDataPointer();
}

template<class T>
void Waveform<T>::setData(const size_t n, const T x[])
{
    if (n < 1 || x == nullptr)
    {
        if (n < 1){RTSEIS_THROW_IA("%s", "x has zero length");}
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "Invalid arguments");
    }
    pImpl->restoreSamplingPeriod();
    pImpl->setData(n, x, true);
}

template<>
std::vector<double> Waveform<double>::getData() const
{
    std::vector<double> y;
    int ny = pImpl->getNumberOfOutputSamples();
    y.resize(ny);
    if (ny > 0)
    {
        const double *yout = pImpl->getOutputDataPointer();
        ippsCopy_64f(yout, y.data(), ny);
    }
    return y;
}

template<class T>
void Waveform<T>::getData(const size_t nwork, T *yIn[]) const
{
    int leny = pImpl->getNumberOfOutputSamples();
    if (nwork < static_cast<size_t> (leny))
    {
        throw std::invalid_argument("nwork = " + std::to_string(nwork)
                                  + " must be at least = "
                                   + std::to_string(leny));
    }
    if (leny == 0){return;}
    T *y = *yIn;
    if (y == nullptr)
    {
        throw std::invalid_argument("y is NULL");
    }
    const T *yout = pImpl->getOutputDataPointer();
    ippsCopy_64f(yout, y, leny);
}

//----------------------------------------------------------------------------//
//                                 Utilities                                  //
//----------------------------------------------------------------------------//

/// TODO delete this function
template<class T>
size_t Waveform<T>::getOutputLength() const
{
    return pImpl->getNumberOfOutputSamples(); //pImpl->ny_;
}

template<class T>
void Waveform<T>::setSamplingPeriod(const double dt)
{
    if (dt <= 0)
    {
        RTSEIS_THROW_IA("Sampling period = %lf must be positive", dt);
    }
    pImpl->dt0_ = dt;
    pImpl->dt_ = dt;
}

template<class T>
double Waveform<T>::getSamplingPeriod() const noexcept
{
    return pImpl->dt_;
}

template<class T>
double Waveform<T>::getNyquistFrequency() const noexcept
{
    double fnyq = 1.0/(2.0*pImpl->dt_);
    return fnyq;
} 

//----------------------------------------------------------------------------//
//                     Convolution/Correlation/AutoCorrelation                //
//----------------------------------------------------------------------------//

template<class T>
void Waveform<T>::convolve(
    const std::vector<T> &s,
    const ConvolutionMode mode,
    const ConvolutionImplementation implementation)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
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
        const T *x = pImpl->getInputDataPointer();
        T *yout = pImpl->getOutputDataPointer();
        Utilities::Math::Convolve::convolve(nx, x,
                                            ny, s.data(),
                                            lenc, &nyout, &yout,
                                            convcorMode, convcorImpl);
#ifndef NDEBUG
        assert(lenc == nyout);
#endif
        pImpl->lfirstFilter_ = false;
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to compute convolution: %s", ia.what());
        pImpl->resizeOutputData(0);
    }
}

template<class T>
void Waveform<T>::correlate(
    const std::vector<T> &s, 
    const ConvolutionMode mode,
    const ConvolutionImplementation implementation)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
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
        const T *x = pImpl->getInputDataPointer();
        T *yout = pImpl->getOutputDataPointer();
        Utilities::Math::Convolve::correlate(nx, x,
                                             ny, s.data(),
                                             lenc, &nyout, &yout,
                                             convcorMode, convcorImpl);
#ifndef NDEBUG
        assert(lenc == nyout);
#endif
        pImpl->lfirstFilter_ = false;
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to compute correlation: %s", ia.what());
        pImpl->resizeOutputData(0);
    }
}

template<class T>
void Waveform<T>::autocorrelate(
    const ConvolutionMode mode,
    const ConvolutionImplementation implementation)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
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
        const T *x = pImpl->getInputDataPointer();
        T *yout = pImpl->getOutputDataPointer();
        Utilities::Math::Convolve::autocorrelate(nx, x,
                                                 lenc, &nyout, &yout,
                                                 convcorMode, convcorImpl);
#ifndef NDEBUG
        assert(lenc == nyout);
#endif
        pImpl->lfirstFilter_ = false;
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to compute autocorrelation: %s", ia.what());
        pImpl->resizeOutputData(0);
    }
}

//----------------------------------------------------------------------------//
//                            Demeaning/detrending                            //
//----------------------------------------------------------------------------//

template<class T>
void Waveform<T>::demean()
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_THROW_RTE("%s", "No data is set on the module");
    }
    // Demean the data
    constexpr auto type = RTSeis::FilterImplementations::DetrendType::CONSTANT;
    try
    {
        RTSeis::FilterImplementations::Detrend<T> demean;
        demean.initialize(type);
        const T *x = pImpl->getInputDataPointer();
        pImpl->resizeOutputData(len);
        T *y = pImpl->getOutputDataPointer();
        demean.apply(len, x, &y);
        pImpl->lfirstFilter_ = false;
/* 
        DemeanParameters parms(RTSeis::Precision::DOUBLE);
        Demean demean(parms);
        const double *x = pImpl->getInputDataPointer();
        pImpl->resizeOutputData(len);
        double  *y = pImpl->getOutputDataPointer();
        demean.apply(len, x, y);
        pImpl->lfirstFilter_ = false;
*/
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("%s", e.what());
        throw std::runtime_error("Algorithmic failure");
    }
}

template<class T>
void Waveform<T>::detrend()
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_THROW_RTE("%s", "No data iset set on the module");
    }
    // Detrend the data
    constexpr auto type = RTSeis::FilterImplementations::DetrendType::LINEAR;
    try  
    {    
        RTSeis::FilterImplementations::Detrend<T> detrend;
        detrend.initialize(type);
        const T *x = pImpl->getInputDataPointer();
        pImpl->resizeOutputData(len);
        T *y = pImpl->getOutputDataPointer();
        detrend.apply(len, x, &y); 
        pImpl->lfirstFilter_ = false;
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("%s", e.what());
        throw std::runtime_error("Algorithmic failure");
    }
/*
    DetrendParameters parms(RTSeis::Precision::DOUBLE);
    Detrend detrend(parms);
    const double *x = pImpl->getInputDataPointer();
    pImpl->resizeOutputData(len);
    double  *y = pImpl->getOutputDataPointer();
    detrend.apply(len, x, y); 
    pImpl->lfirstFilter_ = false;
*/
}

//----------------------------------------------------------------------------//
//                        Downsampling and Decimation                         //
//----------------------------------------------------------------------------//

template<class T>
void Waveform<T>::downsample(const int nq)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    if (nq < 1)
    {
        RTSEIS_THROW_IA("Downsampling factor = %d must be at least 1", nq); 
    }
    // Initialize the downsampler
    RTSeis::FilterImplementations::Downsample
        <RTSeis::ProcessingMode::POST_PROCESSING, T> downsample;
    try
    {
        downsample.initialize(nq);
        // Space estimate
        int leny = downsample.estimateSpace(len);
        pImpl->resizeOutputData(leny);
        T *y = pImpl->getOutputDataPointer(); // Handle on output
        const T *x = pImpl->getInputDataPointer(); // Handle on input
        int nyout;
        downsample.apply(len, x, leny, &nyout, &y);  // Finally downsample
#ifndef NDEBUG
        assert(nyout == leny);
#endif
        pImpl->dt_ = pImpl->dt_*static_cast<T> (nq);
        pImpl->lfirstFilter_ = false;
    }
    catch (const std::runtime_error &ra)
    {
        RTSEIS_ERRMSG("Downsampling failed: %s", ra.what());
    }
}

template<class T>
void Waveform<T>::decimate(const int nq, const int filterLength)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    if (nq < 2)
    {
        RTSEIS_THROW_IA("Downsampling factor = %d must be at least 2", nq);
    }
    if (filterLength < 5)
    {
        RTSEIS_THROW_IA("filterLength = %d must be at least 5", filterLength);
    }
    // Handle odd length so I can remove the phase shift from the FIR filter
    int nfir = filterLength;
    if (nfir%2 == 0){nfir = nfir + 1;}
    RTSeis::FilterImplementations::Decimate<RTSeis::ProcessingMode::POST, T> decimate;
    try
    {
        constexpr bool lremovePhaseShift = true;
        decimate.initialize(nq, nfir, lremovePhaseShift);
        // Space estimate
        int leny = decimate.estimateSpace(len);
        pImpl->resizeOutputData(leny);
        T *y = pImpl->getOutputDataPointer(); // Handle on output
        const T *x = pImpl->getInputDataPointer(); // Handle on input
        int nyout;
        decimate.apply(len, x, leny, &nyout, &y);  // Finally downsample
#ifndef NDEBUG
        assert(nyout == leny);
#endif
        pImpl->dt_ = pImpl->dt_*static_cast<T> (nq);
        pImpl->lfirstFilter_ = false;
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("Decimation failed: %s", e.what());
    }
}

template<class T>
void Waveform<T>::interpolate(const double newSamplingPeriod,
                              const InterpolationMethod method)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data set in module");
        return;
    }
    if (newSamplingPeriod <= 0)
    {
        RTSEIS_THROW_IA("New sampling period = %lf must be positive",
                        newSamplingPeriod);
    }
    const T *x = pImpl->getInputDataPointer(); // Handle on input
    if (method == InterpolationMethod::DFT) 
    {
        // From Matlab documentation: if x as size m with sampling rate dx,
        // then the new sampling rate from interpft dy = dx*m/n.
        // Solving for n we obtain: n = m*(dx/dy) where n > m.
        auto npnew
            = static_cast<int> (len*(pImpl->dt_/newSamplingPeriod) + 0.5);
        pImpl->resizeOutputData(npnew);
        T *y = pImpl->getOutputDataPointer(); // Handle on output
        RTSeis::Utilities::Interpolation::interpft(len, x, npnew, &y);
    }
    else if (method == InterpolationMethod::WEIGHTED_AVERAGE_SLOPES)
    {
        // The spline interpolators are more stringent.  In the Fourier
        // domain extrapolation amounts to wrap around (periodicity).
        // Here, the underlying code will throw if it as asked to extrapolate.
        // So let's approximate then refine the output length.
        std::pair<T, T> xInterval(0, (len - 1)*pImpl->dt_);
        auto npnew
            = static_cast<int> (len*(pImpl->dt_/newSamplingPeriod) + 0.5);
        auto tmax = (npnew - 1)*newSamplingPeriod;
        while (tmax > xInterval.second && npnew > 1)
        {
            npnew = npnew - 1;
            tmax = (npnew - 1)*newSamplingPeriod;
        }
        std::pair<T, T> xIntervalNew(0, (npnew - 1)*newSamplingPeriod);
        // Get pointers
        pImpl->resizeOutputData(npnew);
        T *y = pImpl->getOutputDataPointer(); // Handle on output 
        // Now interpolate 
        RTSeis::Utilities::Interpolation::WeightedAverageSlopes<T> was;
        was.initialize(len, xInterval, x);
        was.interpolate(npnew, xIntervalNew, &y);
    }
    else
    {
        RTSEIS_THROW_IA("%s", "Invalid interpolation type");
    }
    // Update sampling period
    pImpl->dt_ = newSamplingPeriod;
    pImpl->lfirstFilter_ = false;
}


//----------------------------------------------------------------------------//
//                                   Envelope                                 //
//----------------------------------------------------------------------------//

/// FIR-based
template<class T>
void Waveform<T>::firEnvelope(const int nfir)
{
    if (nfir < 1)
    {
        RTSEIS_THROW_IA("Number of FIR coefficients = %d must be positive",
                        nfir);
    }
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int nx = pImpl->getLengthOfInputSignal();
    if (nx < 1){RTSEIS_THROW_IA("%s", "No data is set on the module");}
    pImpl->resizeOutputData(nx);
    T *y = pImpl->getOutputDataPointer(); // Handle on output
    const T *x = pImpl->getInputDataPointer(); // Handle on input
    Transforms::FIREnvelope<RTSeis::ProcessingMode::POST, T>envelope;
    envelope.initialize(nfir);
    envelope.transform(nx, x, &y);
    pImpl->lfirstFilter_ = false;
}

/// FFT-based
template<class T>
void Waveform<T>::envelope()
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int nx = pImpl->getLengthOfInputSignal();
    if (nx < 1){RTSEIS_THROW_IA("%s", "No data is set on the module");}
    Transforms::Envelope<T> envelope;
    envelope.initialize(nx);
    pImpl->resizeOutputData(nx);
    T *y = pImpl->getOutputDataPointer(); // Handle on output
    const T *x = pImpl->getInputDataPointer(); // Handle on input
    envelope.transform(nx, x, &y);
    pImpl->lfirstFilter_ = false;
}
//----------------------------------------------------------------------------//
//                           Band-specific Filters                            //
//----------------------------------------------------------------------------//

template<class T>
void Waveform<T>::iirLowpassFilter(const int order, const double fc,
                                   const IIRPrototype prototype,
                                   const double ripple,
                                   const bool lzeroPhase)
{
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    RTSeis::FilterRepresentations::BA ba;
    pImpl->filterDesigner.designLowpassIIRFilter(
                        order, r, ptype, ripple, ba,
                        FilterDesign::IIRFilterDomain::DIGITAL);
    iirFilter(ba, lzeroPhase);
}

template<class T>
void Waveform<T>::sosLowpassFilter(const int order, const double fc, 
                                   const IIRPrototype prototype,
                                   const double ripple,
                                   const bool lzeroPhase)
{
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    RTSeis::FilterRepresentations::SOS sos;
    pImpl->filterDesigner.designLowpassIIRFilter(
                     order, r, ptype, ripple, sos,
                     FilterDesign::SOSPairing::NEAREST,
                     FilterDesign::IIRFilterDomain::DIGITAL);
    sosFilter(sos, lzeroPhase);
}

template<class T>
void Waveform<T>::firLowpassFilter(const int ntapsIn, const double fc,
                                   const FIRWindow windowIn,
                                   const bool lremovePhase)
{
    int ntaps = ntapsIn;
    if (lremovePhase && ntaps%2 == 0)
    {
        std::cerr << "Adding a filter tap" << std::endl;
        ntaps = ntaps + 1;
    }
    if (ntaps < 5){RTSEIS_THROW_IA("ntaps = %d  must be at least 5", ntaps);}
    int order = ntaps - 1;
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::FIRWindow window;
    window = classifyFIRWindow(windowIn);
    RTSeis::FilterRepresentations::FIR fir;
    pImpl->filterDesigner.designLowpassFIRFilter(order, r, window, fir);
    if (!lremovePhase)
    {
        firFilter(fir);
    }
    else
    {
        FIR_REMOVE_PHASE(pImpl, fir);
    }
}

template<class T>
void Waveform<T>::iirHighpassFilter(const int order, const double fc, 
                                    const IIRPrototype prototype,
                                    const double ripple,
                                    const bool lzeroPhase)
{
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    RTSeis::FilterRepresentations::BA ba; 
    pImpl->filterDesigner.designHighpassIIRFilter(
                     order, r, ptype, ripple, ba, 
                     FilterDesign::IIRFilterDomain::DIGITAL);
    iirFilter(ba, lzeroPhase);
}

template<class T>
void Waveform<T>::sosHighpassFilter(const int order, const double fc, 
                                    const IIRPrototype prototype,
                                    const double ripple,
                                    const bool lzeroPhase)
{
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    RTSeis::FilterRepresentations::SOS sos;
    pImpl->filterDesigner.designHighpassIIRFilter(
                    order, r, ptype, ripple, sos,
                    FilterDesign::SOSPairing::NEAREST,
                    FilterDesign::IIRFilterDomain::DIGITAL);
    sosFilter(sos, lzeroPhase);
}

template<class T>
void Waveform<T>::firHighpassFilter(const int ntapsIn, const double fc, 
                                    const FIRWindow windowIn,
                                    const bool lremovePhase)
{
    int ntaps = ntapsIn;
    if (lremovePhase && ntaps%2 == 0)
    {
        std::cerr << "Adding a filter tap" << std::endl;
        ntaps = ntaps + 1;
    }
    if (ntaps < 5){RTSEIS_THROW_IA("ntaps = %d  must be at least 5", ntaps);}
    int order = ntaps - 1;
    // Compute normalized frequencies
    double r = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::FIRWindow window;
    window = classifyFIRWindow(windowIn);
    RTSeis::FilterRepresentations::FIR fir;
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
}

template<class T>
void Waveform<T>::iirBandpassFilter(const int order,
                                    const std::pair<double,double> fc, 
                                    const IIRPrototype prototype,
                                    const double ripple,
                                    const bool lzeroPhase)
{
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    RTSeis::FilterRepresentations::BA ba; 
    pImpl->filterDesigner.designBandpassIIRFilter(
                    order, r, ptype, ripple, ba, 
                    FilterDesign::IIRFilterDomain::DIGITAL);
    iirFilter(ba, lzeroPhase);
}

template<class T>
void Waveform<T>::sosBandpassFilter(const int order,
                                    const std::pair<double,double> fc, 
                                    const IIRPrototype prototype,
                                    const double ripple,
                                    const bool lzeroPhase)
{
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    RTSeis::FilterRepresentations::SOS sos;
    pImpl->filterDesigner.designBandpassIIRFilter(
                    order, r, ptype, ripple, sos,
                    FilterDesign::SOSPairing::NEAREST,
                    FilterDesign::IIRFilterDomain::DIGITAL);
    sosFilter(sos, lzeroPhase);
}

template<class T>
void Waveform<T>::firBandpassFilter(const int ntapsIn,
                                    const std::pair<double,double> fc, 
                                    const FIRWindow windowIn,
                                    const bool lremovePhase)
{
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
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::FIRWindow window;
    window = classifyFIRWindow(windowIn);
    RTSeis::FilterRepresentations::FIR fir;
    pImpl->filterDesigner.designBandpassFIRFilter(order, r, window, fir);
    if (!lremovePhase)
    {
        firFilter(fir);
    }
    else
    {
        FIR_REMOVE_PHASE(pImpl, fir);
    }
}

template<class T>
void Waveform<T>::iirBandstopFilter(const int order,
                                    const std::pair<double,double> fc, 
                                    const IIRPrototype prototype,
                                    const double ripple,
                                    const bool lzeroPhase)
{
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    RTSeis::FilterRepresentations::BA ba;
    pImpl->filterDesigner.designBandstopIIRFilter(
                    order, r, ptype, ripple, ba,
                    FilterDesign::IIRFilterDomain::DIGITAL);
    iirFilter(ba, lzeroPhase);
}

template<class T>
void Waveform<T>::sosBandstopFilter(const int order,
                                    const std::pair<double,double> fc,
                                    const IIRPrototype prototype,
                                    const double ripple,
                                    const bool lzeroPhase)
{
    // Compute normalized frequencies
    std::pair<double,double> r
        = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::IIRPrototype ptype;
    ptype = classifyIIRPrototype(prototype);
    RTSeis::FilterRepresentations::SOS sos;
    pImpl->filterDesigner.designBandstopIIRFilter(
                    order, r, ptype, ripple, sos,
                    FilterDesign::SOSPairing::NEAREST,
                    FilterDesign::IIRFilterDomain::DIGITAL);
    sosFilter(sos, lzeroPhase);
}

template<class T>
void Waveform<T>::firBandstopFilter(const int ntapsIn,
                                    const std::pair<double,double> fc,
                                    const FIRWindow windowIn,
                                    const bool lremovePhase)
{
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
         = computeNormalizedFrequencyFromSamplingPeriod(fc, pImpl->dt_);
    FilterDesign::FIRWindow window;
    window = classifyFIRWindow(windowIn);
    RTSeis::FilterRepresentations::FIR fir;
    pImpl->filterDesigner.designBandstopFIRFilter(order, r, window, fir);
    if (!lremovePhase)
    {
        firFilter(fir);
    }
    else
    {
        FIR_REMOVE_PHASE(pImpl, fir);
    }
}

//----------------------------------------------------------------------------//
//                               General Filtering                            //
//----------------------------------------------------------------------------//

template<>
void Waveform<double>::firFilter(
    const RTSeis::FilterRepresentations::FIR &fir,
    const bool lremovePhase)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
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
    }
    // Initialize filter
    RTSeis::FilterImplementations::FIRFilter<RTSeis::ProcessingMode::POST, double> firFilter;
    firFilter.initialize(nb, taps.data(),
                   RTSeis::FilterImplementations::FIRImplementation::DIRECT);
    pImpl->resizeOutputData(len);
    // Standard FIR filtering 
    const double *x = pImpl->getInputDataPointer();
    double *yout = pImpl->getOutputDataPointer();
    if (!lremovePhase)
    {
        firFilter.apply(len, x, &yout);
    }
    else
    {
        double *ywork = ippsMalloc_64f(len);
        firFilter.apply(len, x,    &ywork); // Filter forwards
        ippsFlip_64f(ywork, yout,  len);    // Reverse y
        firFilter.apply(len, yout, &ywork); // Filter y backwards
        ippsFlip_64f(ywork, yout,  len);    // Reverse it
        ippsFree(ywork);
    }
    pImpl->lfirstFilter_ = false;
}

template<class T>
void Waveform<T>::iirFilter(const RTSeis::FilterRepresentations::BA &ba,
                            const bool lremovePhase)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
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
    }
    // Initialize filter
    if (!lremovePhase)
    {
        RTSeis::FilterImplementations::IIRFilter
            <ProcessingMode::POST, T> iirFilter;
        iirFilter.initialize(nb, b.data(),
                             na, a.data(),
               RTSeis::FilterImplementations::IIRDFImplementation::DF2_FAST);
        pImpl->resizeOutputData(len);
        const T *x = pImpl->getInputDataPointer();
        T *yout = pImpl->getOutputDataPointer();
        iirFilter.apply(len, x, &yout);
    }
    else
    {
        RTSeis::FilterImplementations::IIRIIRFilter<T> iiriirFilter;
        iiriirFilter.initialize(nb, b.data(),
                                na, a.data());
        pImpl->resizeOutputData(len);
        const T *x = pImpl->getInputDataPointer();
        T *yout = pImpl->getOutputDataPointer();
        iiriirFilter.apply(len, x, &yout);
    }
    pImpl->lfirstFilter_ = false;
}

template<>
void Waveform<double>::sosFilter(
    const RTSeis::FilterRepresentations::SOS &sos,
    const bool lremovePhase)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
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
    RTSeis::FilterImplementations::SOSFilter
        <RTSeis::ProcessingMode::POST, double> sosFilter;
    sosFilter.initialize(ns, bs.data(), as.data());
    pImpl->resizeOutputData(len);
    // Get handles on pointers
    const double *x = pImpl->getInputDataPointer();
    double *yout = pImpl->getOutputDataPointer();
    // Zero-phase filtering needs workspace so that x isn't annihalated
    if (lremovePhase)
    {
        double *ywork = ippsMalloc_64f(len);
        sosFilter.apply(len, x,    &ywork); // Filter forwards
        ippsFlip_64f(ywork, yout,  len);    // Reverse y
        sosFilter.apply(len, yout, &ywork); // Filter y backwards
        ippsFlip_64f(ywork, yout,  len);    // Reverse it
        ippsFree(ywork);
    }
    else
    {
        sosFilter.apply(len, x, &yout);
    }
    pImpl->lfirstFilter_ = false;
}

//----------------------------------------------------------------------------//
//                                Normalization                               //
//----------------------------------------------------------------------------//

template<class T>
void Waveform<T>::normalizeMinMax(const std::pair<double, double> targetRange)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        std::cerr << "No data is set on the module" << std::endl;
        return;
    }
    // Normalize the data
    RTSeis::Utilities::Normalization::MinMax<T> minMax;
    const T *x = pImpl->getInputDataPointer();
    minMax.initialize(len, x, targetRange); // Throws
    pImpl->resizeOutputData(len);
    T *y = pImpl->getOutputDataPointer();
    minMax.apply(len, x, &y);
    pImpl->lfirstFilter_ = false;
}

template<class T>
void Waveform<T>::normalizeSignBit()
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        std::cerr << "No data is set on the module" << std::endl;
        return;
    }
    // Normalize the data
    RTSeis::Utilities::Normalization::SignBit signBit;
    signBit.initialize();
    const T *x = pImpl->getInputDataPointer();
    pImpl->resizeOutputData(len);
    T *y = pImpl->getOutputDataPointer();
    signBit.apply(len, x, &y);
    pImpl->lfirstFilter_ = false;
}

template<class T>
void Waveform<T>::normalizeZScore()
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        std::cerr << "No data is set on the module" << std::endl;
        return;
    }
    // Normalize the data
    RTSeis::Utilities::Normalization::ZScore zscore;
    const T *x = pImpl->getInputDataPointer();
    pImpl->resizeOutputData(len);
    T *y = pImpl->getOutputDataPointer();
    if (len > 1)
    {
        zscore.initialize(len, x); // Throws if all points are identical 
        zscore.apply(len, x, &y);
    }
    else
    {
        y[0] = 0;
    }
    pImpl->lfirstFilter_ = false;
}

//----------------------------------------------------------------------------//
//                                   Tapering                                 //
//----------------------------------------------------------------------------//

template<class T>
void Waveform<T>::taper(const double percentage,
                        const TaperParameters::Type window)
{
    if (!pImpl->lfirstFilter_){pImpl->overwriteInputWithOutput();}
    int len = pImpl->getLengthOfInputSignal();
    if (len < 1)
    {
        std::cerr << "No data is set on the module" << std::endl;
        return;
    }
    // Taper the data
    if (percentage < 0 || percentage > 100)
    {
        throw std::invalid_argument("Percentage = "
                                  + std::to_string(percentage)
                                  + " must be in range [0,100]");
    }
    auto windowType = RTSeis::FilterImplementations::TaperWindowType::Hamming;
    if (window == TaperParameters::Type::HAMMING)
    {
        windowType = RTSeis::FilterImplementations::TaperWindowType::Hamming;
    }
    else if (window == TaperParameters::Type::BARTLETT)
    {
        windowType = RTSeis::FilterImplementations::TaperWindowType::Bartlett;
    }
    else if (window == TaperParameters::Type::HANN)
    {
        windowType = RTSeis::FilterImplementations::TaperWindowType::Hann;
    }
    else if (window == TaperParameters::Type::BLACKMAN)
    {
        windowType = RTSeis::FilterImplementations::TaperWindowType::Blackman;
    } 
    else if (window == TaperParameters::Type::SINE)
    {
        windowType = RTSeis::FilterImplementations::TaperWindowType::Sine;
    }
    else
    {
        throw std::invalid_argument("Unhandled window type");
    }
    RTSeis::FilterImplementations::Taper<T> taper;
    taper.initialize(percentage, windowType);

    const T *x = pImpl->getInputDataPointer();
    pImpl->resizeOutputData(len);
    T *y = pImpl->getOutputDataPointer();
    taper.apply(len, x, &y);
    pImpl->lfirstFilter_ = false;
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
#ifndef NDEBUG
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
    // Classify the convolution implementation
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

FilterDesign::IIRPrototype
classifyIIRPrototype(const IIRPrototype prototype)
{
    FilterDesign::IIRPrototype ptype;
    if (prototype == IIRPrototype::BESSEL)
    {
        ptype = FilterDesign::IIRPrototype::BESSEL;
    }
    else if (prototype == IIRPrototype::BUTTERWORTH)
    {
        ptype = FilterDesign::IIRPrototype::BUTTERWORTH;
    }
    else if (prototype == IIRPrototype::CHEBYSHEV1)
    {
        ptype = FilterDesign::IIRPrototype::CHEBYSHEV1;
    }
    else if (prototype == IIRPrototype::CHEBYSHEV2)
    {
        ptype = FilterDesign::IIRPrototype::CHEBYSHEV2;
    }
    else
    {
        throw std::invalid_argument("Unsupported window = " + 
                                  std::to_string(static_cast<int> (prototype)));
    }
    return ptype;
}

FilterDesign::FIRWindow
classifyFIRWindow(const FIRWindow windowIn)
{
    FilterDesign::FIRWindow window;
    if (windowIn == FIRWindow::HAMMING)
    {
        window = FilterDesign::FIRWindow::HAMMING;
    }
    else if (windowIn == FIRWindow::HANN)
    {
        window = FilterDesign::FIRWindow::HANN;
    }
    else if (windowIn == FIRWindow::BLACKMAN_OPT)
    {
        window = FilterDesign::FIRWindow::BLACKMAN_OPT;
    }
    else if (windowIn == FIRWindow::BARTLETT)
    {
        window = FilterDesign::FIRWindow::BARTLETT;
    }
    else
    {
        throw std::invalid_argument("Unsupported window = " + 
                                   std::to_string(static_cast<int> (windowIn)));
    }
    return window;
}

///--------------------------------------------------------------------------///
///                         Template Instantiation                           ///
///--------------------------------------------------------------------------///
template class PostProcessing::SingleChannel::Waveform<double>;
