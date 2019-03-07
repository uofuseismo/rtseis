#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#ifdef __INTEL_COMPILER
#include <pstl/execution>
#include <pstl/algorithm>
#endif
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/postProcessing/singleChannel/waveform.hpp"
#include "rtseis/postProcessing/singleChannel/detrend.hpp"
#include "rtseis/postProcessing/singleChannel/demean.hpp"
#include "rtseis/postProcessing/singleChannel/taper.hpp"
#include "rtseis/utilities/math/convolve.hpp"

using namespace RTSeis::PostProcessing::SingleChannel;

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
        throw std::invalid_argument("Sampling period = "
                                   + std::to_string(dt)
                                   + " must be postiive");
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
//#ifdef __INTEL_COMPILER
    std::copy(pstl::execution::unseq, pImpl->y_.begin(), pImpl->y_.end(), y);
//#else
    std::copy(pImpl->y_.begin(), pImpl->y_.end(), y);
//#endif 
    return;
}

size_t Waveform::getOutputLength(void) const
{
    return pImpl->y_.size();
}

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
        RTSEIS_ERRMSG("%s", "No data points in s");
        throw std::invalid_argument("s has no data points");
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
