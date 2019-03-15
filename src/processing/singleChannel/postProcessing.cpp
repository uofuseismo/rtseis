#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <exception>
#include <vector>
#include <memory>
#include <algorithm>
#include <ipps.h>
#ifdef __INTEL_COMPILER
#include <pstl/execution>
#include <pstl/algorithm>
#endif
#include "rtseis/private/throw.hpp"
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/postProcessing/singleChannel/waveform.hpp"
#include "rtseis/postProcessing/singleChannel/detrend.hpp"
#include "rtseis/postProcessing/singleChannel/demean.hpp"
#include "rtseis/postProcessing/singleChannel/taper.hpp"
#include "rtseis/utilities/math/convolve.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/utilities/filterRepresentations/sos.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#include "rtseis/utilities/filterImplementations/iirFilter.hpp"
#include "rtseis/utilities/filterImplementations/iiriirFilter.hpp"
#include "rtseis/utilities/filterImplementations/sosFilter.hpp"

using namespace RTSeis::PostProcessing::SingleChannel;

static inline void reverse(std::vector<double> &x);
static inline void reverse(const std::vector<double> &x,
                           std::vector<double> &y);
static inline void copy(const std::vector<double> &x, 
                        std::vector<double> &y);

static inline void copy(const std::vector<double> &x, 
                        std::vector<double> &y)
{
#ifdef __INTEL_COMPILER
    std::copy(pstl::execution::unseq, x.begin(), x.end(), y.begin());
#else
    int len = static_cast<int> (x.size());
    ippsCopy_64f(x.data(), y.data(), len);
#endif
}

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

static inline void reverse(const std::vector<double> &x, 
                           std::vector<double> &y)
{
#ifdef __INTEL_COMPILER
    y.resize(x.size());
    std::reverse_copy(pstl::execution::unseq, x.begin(), x.end(), y.begin());
#else
    int len = static_cast<int> (x.size());
    ippsFlip_64f(x.data(), y.data(), len);
#endif
    return;
}

class Waveform::DataImpl
{
public:
    DataImpl(void)
    {
        return;
    }
    ~DataImpl(void)
    {
        return;
    }
    void resizeOutputData(const int n)
    {
        y_.resize(n);
    }
    void setData(const size_t n, const double x[])
    {
        x_.resize(n);
#ifdef __INTEL_COMPILER
        std::copy(pstl::execution::unseq, x, x+n, x_.data());
#else
        std::copy(x, x+n, x_.data());
#endif
        return;
    }
    const double *getInputDataPointer(void) const
    {
        return x_.data();
    }
    double *getOutputDataPointer(void)
    {
        return y_.data();
    }
    int getLengthOfInputSignal(void) const
    {
        return static_cast<int> (x_.size());
    }
//private:
    std::vector<double> x_; 
    std::vector<double> y_;
    /// Sampling period
    double dt = 1;
};

Waveform::Waveform(const double dt) :
    pImpl(new DataImpl()) 
{
    if (dt <= 0)
    {
        RTSEIS_THROW_IA("Sampling period = %lf must be positive", dt);
    }
    pImpl->dt = dt;
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

void Waveform::setData(const size_t n, const double x[])
{
    if (n < 1 || x == nullptr)
    {
        if (n < 1)
        {
            RTSEIS_ERRMSG("%s", "x has zero length");
            throw std::invalid_argument("x has zero length");
        }
        if (x == nullptr)
        {
            RTSEIS_ERRMSG("%s", "x is NULL");
            throw std::invalid_argument("x is NULL");
        }
        throw std::invalid_argument("Invalid arguments");
    }
    pImpl->setData(n, x);
    return;
}

void Waveform::getData(std::vector<double> &y)
{
    y = pImpl->y_;
    return;
}

void Waveform::getData(const size_t nwork, double y[]) const
{
    size_t leny = getOutputLength();
    if (nwork < leny)
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
#ifdef __INTEL_COMPILER
    std::copy(pstl::execution::unseq, pImpl->y_.begin(), pImpl->y_.end(), y);
#else
    std::copy(pImpl->y_.begin(), pImpl->y_.end(), y);
#endif 
    return;
}

size_t Waveform::getOutputLength(void) const
{
    return pImpl->y_.size();
}

//----------------------------------------------------------------------------//
//                     Convolution/Correlation/AutoCorrelation                //
//----------------------------------------------------------------------------//

void Waveform::convolve(
    const std::vector<double> &s,
    const Utilities::Math::Convolve::Mode mode,
    const Utilities::Math::Convolve::Implementation implementation)
{
    int nx = pImpl->getLengthOfInputSignal();
    int ny = static_cast<int> (s.size());
    if (nx < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    } 
    if (ny < 1)
    {
        RTSEIS_THROW_IA("%s", "No data points in s");
    }
    int ierr = Utilities::Math::Convolve::convolve(pImpl->x_, s, pImpl->y_,
                                                   mode, implementation);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute convolution");
        pImpl->y_.resize(0);
        return;
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
        throw std::invalid_argument("Algorithmic failure");
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
//                               General Filtering                            //
//----------------------------------------------------------------------------//

void Waveform::filter(const Utilities::FilterRepresentations::FIR &fir,
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
    // Zero-phase filtering needs workspace so that x isn't annihalated
    if (lremovePhase)
    {
        std::vector<double> xwork(len);
        copy(pImpl->x_, xwork);
        firFilter.apply(len, xwork.data(), pImpl->y_.data());
        reverse(pImpl->y_, xwork);
        firFilter.apply(len, pImpl->x_.data(), pImpl->y_.data());
        reverse(pImpl->y_);
    }
    else
    {
        firFilter.apply(len, pImpl->x_.data(), pImpl->y_.data());
    }
    return;
}

void Waveform::filter(const Utilities::FilterRepresentations::BA &ba,
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
        iirFilter.apply(len, pImpl->x_.data(), pImpl->y_.data());
    }
    else
    {
        RTSeis::Utilities::FilterImplementations::IIRIIRFilter iiriirFilter;
        iiriirFilter.initialize(nb, b.data(),
                                na, a.data(),
                                Precision::DOUBLE);
        pImpl->resizeOutputData(len);
        iiriirFilter.apply(len, pImpl->x_.data(), pImpl->y_.data());
    }
    return;
}

void Waveform::filter(const Utilities::FilterRepresentations::SOS &sos,
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
    // Zero-phase filtering needs workspace so that x isn't annihalated
    if (lremovePhase)
    {
        std::vector<double> xwork(len);
        copy(pImpl->x_, xwork);
        sosFilter.apply(len, xwork.data(), pImpl->y_.data());
        reverse(pImpl->y_, xwork);
        sosFilter.apply(len, pImpl->x_.data(), pImpl->y_.data());
        reverse(pImpl->y_);
    }
    else
    {
        sosFilter.apply(len, pImpl->x_.data(), pImpl->y_.data());
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
