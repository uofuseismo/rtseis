#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <vector>
#include <chrono>
#include <ipps.h>
//include<mkl.h>
#include "rtseis/utilities/interpolation/interpolate.hpp"
#include "rtseis/utilities/interpolation/cubicSpline.hpp"
#include "rtseis/utilities/interpolation/weightedAverageSlopes.hpp"
#include <gtest/gtest.h>

namespace
{

double uniformRandom(const double xMin, const double xMax)
{
    auto r = static_cast<double> (rand())/RAND_MAX; // [0,1]
    r = xMin + (xMax - xMin)*r;
    r = std::min(xMax, std::max(xMin, r));
    return r;
}
//void spline(const std::vector<double> &x,
//            const std::vector<double> &y,
//            const std::vector<double> &xq);

using namespace RTSeis::Utilities::Interpolation;

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

TEST(UtilitiesInterpolation, linearInterolation)
{

}

TEST(UtilitiesInterpolation, cubicSpline)
{
    int nq = 101;
    // Load the right answers
    const std::string interpFileName = "data/cubicSplineInterpolation.txt";
    std::vector<double> yqNotAKnotRef(nq);
    std::vector<double> yqNaturalRef(nq);
    std::vector<double> yqClampedRef(nq);
    std::vector<double> yqPeriodicRef(nq);
    std::ifstream textFile(interpFileName); 
    std::string line;
    int i = 0;
    while (std::getline(textFile, line))
    {
        double xval, knot, natural, clamped, periodic;
        std::sscanf(line.c_str(), "%lf, %lf, %lf, %lf, %lf\n",
                    &xval, &knot, &clamped, &natural, &periodic);
        yqNotAKnotRef[i] = knot;
        yqNaturalRef[i] = natural;
        yqClampedRef[i] = clamped;
        yqPeriodicRef[i] = periodic;
        i = i + 1;
    }
    textFile.close();
    EXPECT_EQ(i, nq);

    // Let's interpolate a sine wave
    int npts = 21;
    std::vector<double> x(npts);
    std::vector<double> y(npts);
    for (auto i=0; i<npts; ++i)
    {
        x[i] = static_cast<double> (i)*M_PI/5.0;
        y[i] = sin(x[i]);
    }
    y[0] = 0.01; // make life a little more exciting
    y[npts-1] = y[0]; // periodic b.c.'s are pretty stringent
    // Choose some interpolation points
    double xMin = x[0];
    double xMax = x[npts-1];
    //int nq = 101;
    double dx = (xMax - xMin)/static_cast<double> (nq - 1);
    std::vector<double> xq(nq);
    std::vector<double> yqNotAKnot(nq, 0.0);
    std::vector<double> yqNatural(nq, 0.0);
    std::vector<double> yqClamped(nq, 0.0);
    std::vector<double> yqPeriodic(nq, 0.0);
    for (auto i=0; i<nq; ++i)
    {
        xq[i] = xMin + i*dx;
    }
    xq[0] = 0;
    xq[nq-1] = xMax;
    double *yPtr = nullptr;
    // Create the spline
    CubicSpline spline; 
    EXPECT_NO_THROW(spline.initialize(npts, x.data(), y.data(),
                                CubicSplineBoundaryConditionType::NOT_A_KNOT));
    yPtr = yqNotAKnot.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr));
    // natural
    EXPECT_NO_THROW(spline.initialize(npts, x.data(),  y.data(),
                               CubicSplineBoundaryConditionType::NATURAL));
    yPtr = yqNatural.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr)); 
    // clamped
    EXPECT_NO_THROW(spline.initialize(npts, x.data(),  y.data(),
                               CubicSplineBoundaryConditionType::CLAMPED));
    yPtr = yqClamped.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr));
    // periodic
    EXPECT_NO_THROW(spline.initialize(npts, x.data(),  y.data(),
                               CubicSplineBoundaryConditionType::PERIODIC));
    yPtr = yqPeriodic.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr));

    double error;
    ippsNormDiff_Inf_64f(yqNotAKnot.data(), yqNotAKnotRef.data(), nq, &error);
    EXPECT_LE(error, 1.e-7);
    ippsNormDiff_Inf_64f(yqNatural.data(),  yqNaturalRef.data(),  nq, &error);
    EXPECT_LE(error, 1.e-7);
    ippsNormDiff_Inf_64f(yqClamped.data(),  yqClampedRef.data(),  nq, &error);
    EXPECT_LE(error, 1.e-7);
    ippsNormDiff_Inf_64f(yqPeriodic.data(), yqPeriodicRef.data(), nq, &error);
    EXPECT_LE(error, 1.e-7);

    // It checks out - copy the results and
    yqNotAKnotRef = yqNotAKnot;
    yqNaturalRef  = yqNatural;
    yqClampedRef  = yqClamped;
    yqPeriodicRef = yqPeriodic;
    // Repeat the exercise for uniform partition
    std::pair<double, double> xDomain(xMin, xMax);
    EXPECT_NO_THROW(spline.initialize(npts, xDomain, y.data(),
                                      CubicSplineBoundaryConditionType::NOT_A_KNOT));
    yPtr = yqNotAKnot.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr));
    EXPECT_NEAR(spline.integrate(std::pair<double, double> (xMin, 1)),
                0.4625836825854708, 1.e-12);
    // natural
    EXPECT_NO_THROW(spline.initialize(npts, xDomain,  y.data(),
                               CubicSplineBoundaryConditionType::NATURAL));
    yPtr = yqNatural.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr));
    EXPECT_NEAR(spline.integrate(std::pair<double, double> (xMax, 0)),
                -0.004955392017811033, 1.e-12);
    // clamped
    EXPECT_NO_THROW(spline.initialize(npts, xDomain,  y.data(),
                               CubicSplineBoundaryConditionType::CLAMPED));
    yPtr = yqClamped.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr));
    EXPECT_NEAR(spline.integrate(std::pair<double, double> (4.08, 6)),
                -1.5508914757864318, 1.e-12);
    // periodic
    EXPECT_NO_THROW(spline.initialize(npts, xDomain,  y.data(),
                               CubicSplineBoundaryConditionType::PERIODIC));
    yPtr = yqPeriodic.data();
    EXPECT_NO_THROW(spline.interpolate(nq, xq.data(), &yPtr));
    EXPECT_NEAR(spline.integrate(std::pair<double, double> (3.99, xMax)),
                -1.657622164641667, 1.e-12);

    // Compute the differences
    ippsNormDiff_Inf_64f(yqNotAKnot.data(), yqNotAKnotRef.data(), nq, &error);
    EXPECT_LE(error, 1.e-14);
    ippsNormDiff_Inf_64f(yqNatural.data(),  yqNaturalRef.data(),  nq, &error);
    EXPECT_LE(error, 1.e-14);
    ippsNormDiff_Inf_64f(yqClamped.data(),  yqClampedRef.data(),  nq, &error);
    EXPECT_LE(error, 1.e-14);
    ippsNormDiff_Inf_64f(yqPeriodic.data(), yqPeriodicRef.data(), nq, &error);
    EXPECT_LE(error, 1.e-14);

    // Now let's test the quadrature of many points
    auto nIntervals = 1000;
    srand(89023);
    std::vector<std::pair<double,double>> intervals(nIntervals);
    std::vector<double> integralsRef(nIntervals);
    for (auto i=0; i<nIntervals; ++i)
    {
        auto r1 = uniformRandom(xMin, xMax);
        auto r2 = uniformRandom(xMin, xMax);
        if (i%20 == 0){r2 = r1;}
        std::pair<double, double> r(r1, r2);
        integralsRef[i] = spline.integrate(r);
        intervals[i] = r;
    }
    std::vector<double> integrals(nIntervals, 0);
    double *intPtr = integrals.data();
    spline.integrate(nIntervals, intervals.data(), &intPtr);
    ippsNormDiff_Inf_64f(integrals.data(), integralsRef.data(),
                         nIntervals, &error);
    EXPECT_LE(error, 1.e-14);
}

TEST(UtilitiesInterpolation, weighedAverageSlopes)
{
    // Load the data
    double dt = 1.0/200.0;
    double dtnew = 1.0/250.0;
    std::string dataFile = "data/gse2.txt";    
    std::string refFile = "data/wigint.txt";
    // Get input signal
    std::ifstream textFile(dataFile);
    std::string line;
    int i = 0;
    std::vector<double> y;
    y.reserve(12000);
    while (std::getline(textFile, line))
    {
        double yval;
        std::sscanf(line.c_str(), "%lf\n", &yval);
        y.push_back(yval);
        i = i + 1;
    }
    textFile.close();
    EXPECT_EQ(static_cast<int> (y.capacity()), i);
    EXPECT_EQ(y.size(), y.capacity());
    // Get reference signal
    textFile.open(refFile);
    std::vector<double> xq, yIntRef;
    int nqref = 14999;
    xq.reserve(nqref);
    yIntRef.reserve(nqref);
    i = 0;
    while (std::getline(textFile, line))
    {
        double xin, yin;
        std::sscanf(line.c_str(), "%lf, %lf\n", &xin, &yin);
        xq.push_back(xin);
        yIntRef.push_back(yin);
        i = i + 1;
    }
    textFile.close();
    EXPECT_EQ(yIntRef.capacity(), i);
    EXPECT_EQ(yIntRef.size(), yIntRef.capacity());

    std::pair<double, double> xInterval(0, dt*(y.size() - 1));
    // Initialize
    WeightedAverageSlopes<double> slopes;
    EXPECT_NO_THROW(slopes.initialize(y.size(), xInterval, y.data()));
    // Compute number of new points 
    auto tmax = static_cast<double> (y.size() - 1)*dt;
    auto nq = static_cast<int> (tmax/dtnew + 0.5);
    EXPECT_EQ(nq, nqref);
    std::vector<double> yInt(nq);
    double *yIntPtr = yInt.data();
    slopes.interpolate(nq, xq.data(), &yIntPtr); 
    // Tabulate the error
    double error = 0;
    ippsNormDiff_Inf_64f(yIntRef.data(), yInt.data(), nq, &error); 
    EXPECT_LE(error, 1.e-8);
    // Interpolate on the evenly spaced interval
    std::pair<double, double> xIntervalNew(0, dtnew*(nq - 1));
    slopes.interpolate(nq, xIntervalNew, &yIntPtr);
    ippsNormDiff_Inf_64f(yIntRef.data(), yInt.data(), nq, &error);
    EXPECT_LE(error, 1.e-8);
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
