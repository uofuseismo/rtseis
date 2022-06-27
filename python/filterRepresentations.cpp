#include <vector>
#include <memory>
#include <exception>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "modules.hpp"

namespace
{
class FIR
{
public:
    FIR() :
        mFilter(std::make_unique<RTSeis::FilterRepresentations::FIR> ())
    {
    }
    FIR(const FIR &fir)
    {
        *this = fir;
    }
    FIR(FIR &&fir) noexcept
    {
        *this = std::move(fir);
    }
    ~FIR() = default;
    FIR& operator=(const FIR &fir)
    {
        if (&fir == this){return *this;}
        mFilter = std::make_unique<RTSeis::FilterRepresentations::FIR> (*fir.mFilter);
        return *this;
    }
    FIR& operator=(FIR &&fir) noexcept
    {
        if (&fir == this){return *this;}
        mFilter = std::move(fir.mFilter);
        return *this;
    }
    RTSeis::FilterRepresentations::FIR getNativeClass() const noexcept
    {
        return *mFilter;
    }
    void setFilterTaps(const std::vector<double> &b)
    {
        mFilter->setFilterTaps(b);
    }
    std::vector<double> getFilterTaps() const noexcept
    {
        return mFilter->getFilterTaps();
    }
    std::unique_ptr<RTSeis::FilterRepresentations::FIR> mFilter;
};
}

void PFilterRepresentations::initialize(pybind11::module &m)
{
    pybind11::class_<::FIR> fir(m, "FIR");
    fir.def(pybind11::init<> ());
    fir.doc() = R""""(
Defines the finite-impulse response (FIR) filter structure.

Parameters :
    taps : The filter taps.  The order of the filter is len(taps) - 1.
)"""";
    fir.def_property("taps",
                     &::FIR::getFilterTaps,
                     &::FIR::setFilterTaps);
    fir.def("__copy__", [](const ::FIR &self)
    {
        return ::FIR(self);
    });
    fir.def("__deepcopy_-", [](const ::FIR &self, pybind11::dict)
    {
        return ::FIR(self);
    });
}

