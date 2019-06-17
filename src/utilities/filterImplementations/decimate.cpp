#include <cstdio>
#include <cstdlib>
#include "rtseis/utilities/filterImplementations/decimate.hpp"
#include "rtseis/utilities/filterImplementations/multiRateFIRFilter.hpp"

using namespace RTSeis::Utilities::FilterImplementations;

class Decimate::DecimateImpl
{
public:
    
};

Decimate::Decimate() :
    pImpl(std::make_unique<DecimateImpl> ())
{
}

Decimate::~Decimate() = default;


