#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <vector>
#include <chrono>
#include "rtseis/utilities/math/interpolate.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis::Utilities::Math::Interpolation;

//int test_interpolation_interpft(void)
TEST(UtilitiesInterpolation, interpft)
{
    int npts = 31;
    int npnew = npts*7;
    auto dx = 3.0*M_PI/static_cast<double> (npts - 1);
    std::vector<double> f(npts);
    for (auto i=0; i<npts; i++)
    {
        auto x = static_cast<double> (i)*dx;
        f[i] = std::pow(std::sin(x), 2)*std::cos(x);
    }
    std::vector<double> fnew;
    EXPECT_NO_THROW(fnew = interpft(f, npnew));
    int iskip = npnew/npts;
    double emax = 0;
    for (auto i=0; i<npts; i++)
    {
        emax = std::max(emax, std::abs(fnew[i*iskip] - f[i]));
        if (std::abs(fnew[i*iskip] - f[i]) > 1.e-12)
        {
            fprintf(stderr, "ft interp failed %lf %lf", f[i], fnew[i*iskip]);
        }
    }
    EXPECT_LE(emax, 1.e-12);
    double dy = dx/7.;
    emax = 0;
    for (int i=0; i<npnew; i++)
    {
        auto x = static_cast<double> (i)*dy;
        auto ftrue = std::pow(std::sin(x), 2)*std::cos(x);
        emax = std::max(emax, std::abs(ftrue - fnew[i]));
        if (std::abs(ftrue - fnew[i]) > 1.e-1)
        {
            fprintf(stderr, "Failed more aggressive check %lf %lf %lf %lf\n",
                    x, ftrue, fnew[i], ftrue - fnew[i]);
        }
    }
    EXPECT_LE(emax, 1.e-1);
}

}
