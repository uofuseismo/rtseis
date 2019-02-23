#include <vector>
#include <memory>
#include <exception>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include "include/rtseis/postProcessing/singleChannel/waveform.hpp"

namespace py = pybind11;

/*
int add(int i, int j) {
    return i + j;
}
*/

class SingleChannelWaveform
{
    public:
        SingleChannelWaveform(void) :
           waveform_(new RTSeis::PostProcessing::SingleChannel::Waveform())
        {
            return;
        }
        ~SingleChannelWaveform(void){return;}
        /// Remove mean
        void demean(void)
        {
            try
            {
                waveform_->demean();
            }
            catch (const std::invalid_argument& ia) 
            {
                fprintf(stderr, "%s\n", ia.what());
                throw std::invalid_argument("Demean failed");
            }
            return;
        }
        /// Remove trend
        void detrend(void)
        {
            try
            {
                waveform_->detrend();
            }
            catch (const std::invalid_argument& ia)
            {
                fprintf(stderr, "%s\n", ia.what());
                throw std::invalid_argument("Detrend failed");
            }
            return;
        }
        /// Sets the data
        void setData(py::array_t<double> x)
        {
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
        py::array_t<double> getData(void)
        {
            std::vector<double> y;
            waveform_->getData(y);
            return py::array(y.size(), y.data());
        }
        bool isInitialized(void) const
        {
            return true;//waveform_->isInitialized();
        }
    private:
        std::unique_ptr<RTSeis::PostProcessing::SingleChannel::Waveform> waveform_;
};

PYBIND11_MODULE(libpyrtseis, m)
{
    py::class_<SingleChannelWaveform> singleChannelWaveform(m, "SingleChannelWaveform");
    //singleChannelWaveform
    singleChannelWaveform.def(py::init<>());
    singleChannelWaveform.def("setData", &SingleChannelWaveform::setData,
                              "Sets the signal to process on the class");
    singleChannelWaveform.def("getData", &SingleChannelWaveform::getData,
                              "Gets the filtered data as a NumPy array");
    singleChannelWaveform.def("demean", &SingleChannelWaveform::demean,
                              "Removes the mean from the time series");
    singleChannelWaveform.def("detrend", &SingleChannelWaveform::detrend,
                              "Removes the trend from the time series");
    singleChannelWaveform.def("isInitialized",
                              &SingleChannelWaveform::isInitialized,
                              "Checks if the class is initialized");

    //m.doc() = "pybind11 example plugin"; // optional module docstring

    //m.def("add", &add, "A function which adds two numbers");
}
