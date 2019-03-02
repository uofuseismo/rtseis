#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <vector>
#include <chrono>
#define RTSEIS_LOGGING 1
#include "utils.hpp"
#include "rtseis/utilities/math/interpolate.hpp"
#include "rtseis/log.h"

int test_interpolation_interpft(void);


int rtseis_test_utils_interpolation(void)
{
    int ierr;
    ierr = test_interpolation_interpft();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed interpft test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed interpft test");
    return EXIT_SUCCESS;
}

int test_interpolation_interpft(void)
{
    int npts = 31;
    int npnew = npts*7;
    double dx = 3.0*M_PI/static_cast<double> (npts - 1);
    std::vector<double> f(npts);
    for (int i=0; i<npts; i++)
    {
        double x = static_cast<double> (i)*dx;
        f[i] = std::pow(std::sin(x), 2)*std::cos(x);
    }   
    std::vector<double> fnew;
    //double dy = 3.0*M_PI/static_cast<double> (npnew - 1);
    try
    {
        RTSeis::Utilities::Math::Interpolate::interpft(f, npnew, fnew);
    }
    catch (std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Interpolation failed: %s", ia.what());
        return EXIT_FAILURE;
    }
    int iskip = npnew/npts;
    for (int i=0; i<npts; i++)
    {
        if (std::abs(fnew[i*iskip] - f[i]) > 1.e-12)
        {
            RTSEIS_ERRMSG("ft interp failed %lf %lf", f[i], fnew[i*iskip]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
