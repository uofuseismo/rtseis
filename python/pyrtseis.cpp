//#include "modules.hpp"
#include <rtseis/version.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

PYBIND11_MODULE(pyrtseis, modules)
{
    modules.attr("__version__") = RTSEIS_VERSION;
    modules.attr("__name__") = "pyrtseis"; 
    modules.attr("__doc__") = "The Python interface to RTSeis.";
    //------------------------------------------------------------------------//
    //                         PostProcessing Group                           //
    //------------------------------------------------------------------------// 
}
