#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "rtseis/utilities/rotate/utilities.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace RTSeis::Utilities::Rotate;

TEST(UtilitiesRotate, ne2rt)
{
    int npts = 100;
    std::vector<double> north(npts), east(npts);
    std::vector<double> radial(npts), transverse(npts);
    
}

}
