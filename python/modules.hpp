#ifndef PYRTSEIS_MODULES_HPP
#define PYRTSEIS_MODULES_HPP
#include <pybind11/pybind11.h>

void init_pp_waveform(pybind11::module &m);
namespace PTransforms
{
[[maybe_unused]]
void initialize(pybind11::module &m);
}
namespace PFilterRepresentations
{
[[maybe_unused]]
void initialize(pybind11::module &m);
}

#endif
