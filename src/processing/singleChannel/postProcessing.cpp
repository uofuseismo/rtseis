#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
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
            std::copy(x, x+n, x_.data());
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
};

Waveform::Waveform(void) :
    pData_(new DataImpl()) 
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
    pData_->setData(n, x);
    return;
}

void Waveform::getData(std::vector<double> &y)
{
    y = pData_->y_;
    return;
}

void Waveform::convolve(
    const std::vector<double> &s,
    const Utilities::Math::Convolve::Mode mode,
    const Utilities::Math::Convolve::Implementation implementation)
{
    int nx = pData_->getLengthOfInputSignal();
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
    int ierr = Utilities::Math::Convolve::convolve(pData_->x_, s, pData_->y_,
                                                   mode, implementation);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute convolution");
        pData_->y_.resize(0);
        return;
    }
    return;
}

void Waveform::demean(void)
{
    int len = pData_->getLengthOfInputSignal();
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
        const double *x = pData_->getInputDataPointer();
        pData_->resizeOutputData(len);
        double  *y = pData_->getOutputDataPointer();
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
    int len = pData_->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "At least 2 data points required to detrend");
        return;
    }
    // Detrend the data
    DetrendParameters parms(RTSeis::Precision::DOUBLE);
    Detrend detrend(parms);
    const double *x = pData_->getInputDataPointer();
    pData_->resizeOutputData(len);
    double  *y = pData_->getOutputDataPointer();
    detrend.apply(len, x, y); 
    return;
}

void Waveform::taper(const double pct,
                     const TaperParameters::Type window)
{
    int len = pData_->getLengthOfInputSignal();
    if (len < 1)
    {
        RTSEIS_WARNMSG("%s", "No data is set on the module");
        return;
    }
    // Taper the data
    TaperParameters parms(pct, window, RTSeis::Precision::DOUBLE);
    Taper taper(parms);
    const double *x = pData_->getInputDataPointer();
    pData_->resizeOutputData(len);
    double  *y = pData_->getOutputDataPointer();
    taper.apply(len, x, y);
    return;
}
