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
using namespace PBFilterRepresentations;

FIR::FIR(void) :
   fir_(new RTSeis::Utilities::FilterRepresentations::FIR()) 
{
   return;
}
/// Set the filter coefficients
void FIR::setTaps(py::array_t<double, py::array::c_style | py::array::forcecast> taps)
{
   std::vector<double> tapsVec(taps.size());
   std::memcpy(tapsVec.data(), tapsVec.data(), taps.size()*sizeof(double));
   fir_->setFilterTaps(tapsVec);
}
/// Get the filter coefficients in a form usable
py::array_t<double> FIR::getTaps(void) const
{
   std::vector<double> taps;
   taps = fir_->getFilterTaps();
   return py::array(taps.size(), taps.data());
}

/*
class BA
{
public:
    BA(void) :
       ba_(new RTSeis::Utilities::FilterRepresentations::BA())
    {
        return;
    }
 
private:
    std::unique_ptr<RTSeis::Utilities::FilterRepresentations::BA> ba_;
};
};
};

*/
