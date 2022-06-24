#include <iostream>
#include <vector>
#include <pybind11/stl.h>
#include "rtseis/transforms/continuousWavelet.hpp"
#include "rtseis/transforms/wavelets/morlet.hpp"
#include "modules.hpp"

using namespace PTransforms;

namespace
{
class Morlet
{
public:
    std::unique_ptr<RTSeis::Transforms::Wavelets::Morlet> mWavelet; 
};

class CWT
{
public:
    /// Constructor
    CWT() :
        mCWT(std::make_unique<RTSeis::Transforms::ContinuousWavelet<double>> ())
    {
    }
    /// Destructor
    ~CWT() = default;
    std::unique_ptr<RTSeis::Transforms::ContinuousWavelet<double>> mCWT;
}; 
}

void PTransforms::init_transforms(pybind11::module &m)
{

}
