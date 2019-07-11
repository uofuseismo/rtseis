#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <mkl.h>
#include "rtseis/utilities/math/interpolate.hpp"
#include "rtseis/utilities/math/cubicSpline.hpp"
#include <gtest/gtest.h>

namespace
{

void spline(const std::vector<double> &x,
            const std::vector<double> &y,
            const std::vector<double> &xq);

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

TEST(UtilitiesInterpolation, cubicSpline)
{
    // Let's interpolate a sine wave
    int npts = 10;
    std::vector<double> x(npts);
    std::vector<double> y(npts);
    for (auto i=0; i<npts; ++i)
    {
        x[i] = static_cast<double> (i);
        y[i] = sin(x[i]);        
    }
    // Choose some interpolation points
    double xMin = 0.1;
    double xMax = 8.5;
    double dx = 0.1;
    int nq = static_cast<int> ((xMax - xMin)/dx + 0.5) + 1;
    std::vector<double> xq(nq);
    std::vector<double> yq(nq);
    for (auto i=0; i<nq; ++i)
    {
        xq[i] = xMin + i*dx;
        yq[i] = 0;
    }
    // Create the spline
    CubicSpline spline; 
    EXPECT_NO_THROW(spline.initialize(npts, x.data(), y.data(),
                                CubicSplineBoundaryConditionType::NOT_A_KNOT));//NATURAL));
    double *yPtr = yq.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr));
    for (auto i=0; i<nq; ++i)
    {
        printf("%lf %lf\n", sin(xq[i]), yq[i]);
    }
}

/*
void spline(const std::vector<double> &xIn,
            const std::vector<double> &yIn,
            const std::vector<double> &xq)
{
    double *x = (double *) calloc(xIn.size(), sizeof(double));
    double *y = (double *) calloc(yIn.size(), sizeof(double));
    std::copy(xIn.begin(), xIn.end(), x);
    std::copy(yIn.begin(), yIn.end(), y);
    MKL_INT nx = xIn.size();
    MKL_INT ny = 1;//yIn.size();
    DFTaskPtr task;
    auto status = dfdNewTask1D(&task, nx, x, DF_NO_HINT,
                               ny, y, DF_NO_HINT);
    if (status != DF_STATUS_OK){throw std::runtime_error("failed new task");}
    auto splineOrder = DF_PP_CUBIC;
    auto splineType  = DF_PP_NATURAL;
    auto splineBC    = DF_BC_FREE_END;
    auto splineIC    = DF_NO_IC;
    auto nwork = ny*(nx - 1)*splineOrder;
    double *bc = NULL;
    double *ic = NULL;
    double *scoeff = (double *) calloc(nwork, sizeof(double)); //, 0.0);
    status = dfdEditPPSpline1D(task, splineOrder, splineType, splineBC, bc,
                               splineIC, ic, scoeff, DF_NO_HINT);
    if (status != DF_STATUS_OK){throw std::runtime_error("failed edit spline");}
    status = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
    if (status != DF_STATUS_OK){throw std::runtime_error("construct failed");}
    MKL_INT nsite = xq.size();
    MKL_INT *cell = (MKL_INT *) calloc(nsite, sizeof(MKL_INT));
    status = dfdSearchCells1D(task, DF_METHOD_STD, nsite, xq.data(),
                              DF_NO_HINT, DF_NO_APRIORI_INFO, cell);
    if (status != DF_STATUS_OK){throw std::runtime_error("failed search");}
    MKL_INT norder = 1;
    MKL_INT dorder[1] = {0};
    double *vq = (double *) calloc(nsite, sizeof(double));
    status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
                              nsite, xq.data(),
                              DF_SORTED_DATA, norder, dorder,
                              DF_NO_APRIORI_INFO, vq,
                              DF_MATRIX_STORAGE_ROWS, cell);

    status = dfDeleteTask(&task);
    if (status != DF_STATUS_OK){throw std::runtime_error("failed delete task");}
}
*/

}
