#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <chrono>
#include <ipps.h>
#include "rtseis/utilities/normalization/minMax.hpp"
#include "rtseis/utilities/normalization/signBit.hpp"
#include "rtseis/utilities/normalization/zscore.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis::Utilities::Normalization;

//int test_normalization_signBit(void)
TEST(UtilitiesNormalization, signBit)
{
    SignBit signBit;
    signBit.initialize();
    // Create a random signal
    int npts = 6500;
    std::vector<double> x(npts);
    for (auto i=0; i<npts; i++)
    {   
        x[i] = (static_cast<double> (rand())/RAND_MAX - 0.5)*100;
    }
    std::vector<double> y(npts);
    std::vector<double> yref(npts);
    auto timeStart = std::chrono::high_resolution_clock::now();
    double *yPtr = y.data();
    EXPECT_NO_THROW(signBit.apply(npts, x.data(), &yPtr));
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double emax;
    for (auto i=0; i<npts; i++)
    {
        double sign =-1;
        if (x[i] >= +0){sign =+1;}
        yref[i] = sign;
    }
    ippsNormDiff_Inf_64f(y.data(), yref.data(), npts, &emax);
    EXPECT_LE(emax, 1.e-10);
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    SignBit signBitRT = signBit; 
    // Do a real-time test
    std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                 1000, 1024, 1200, 2048, 4000, 4096, 5000});
    for (auto job=0; job<2; job++)
    {
        for (auto ip=0; ip<static_cast<int> (packetSize.size()); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0;
            int nptsPass = 0;
            while (nxloc < npts)
            {
                nptsPass = packetSize[ip];
                if (job == 1)
                {
                     nptsPass = std::max(1, nptsPass + rand()%50 - 25);
                }
                nptsPass = std::min(nptsPass, npts - nxloc);
                double *yPtr = static_cast<double *> (y.data()) + nxloc;
                EXPECT_NO_THROW(signBitRT.apply(nptsPass, &x[nxloc], &yPtr));//&y[nxloc]));
                nxloc = nxloc + nptsPass;
            }
            signBitRT.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            ippsNormDiff_Inf_64f(yref.data(), y.data(), npts, &emax);
            EXPECT_LE(emax, 1.e-10);
            if (job == 0)
            {
                fprintf(stdout,
                        "Passed signBit filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed onebit filter random in %.8e (s)\n", tdif.count());
            }
        }
    }
}

TEST(UtilitiesNormalization, zscore)
{
    ZScore zscore;
    // Matlab will do z=(x-mean(x))/std(x) where std is standard deviation
    int npts = 500;
    std::vector<double> x(npts);
    double mean = 0;
    for (auto i=0; i<npts; ++i)
    {
        x[i] = static_cast<double>(rand())/RAND_MAX;
        mean = mean + x[i];
    }
    mean = mean/static_cast<double> (npts);
    double var = 0;
    for (auto i=0; i<npts; ++i)
    {
        var = var + (x[i] - mean)*(x[i] - mean);
    }
    var = var/static_cast<double> (npts - 1);
    auto std = std::sqrt(var);
    // Initialize
    EXPECT_NO_THROW(zscore.initialize(mean, std));
    EXPECT_TRUE(zscore.isInitialized());
    std::vector<double> y(npts);
    double *yPtr = y.data();
    EXPECT_NO_THROW(zscore.apply(npts, x.data(), &yPtr));
    // Verify
    double error = 0;
    for (auto i=0; i<npts; ++i)
    {
        double z = (x[i] - mean)/std;
        error = std::max(error, std::abs(z - y[i]));
    }
    EXPECT_LE(error, 1.e-14);
    // Let's do this the other way
    auto zref = y;
    ZScore zscoreAlternate;
    EXPECT_NO_THROW(zscoreAlternate.initialize(npts, x.data()));
    zscore = zscoreAlternate;
    std::fill(y.begin(), y.end(), 0.0);
    yPtr = y.data();
    EXPECT_NO_THROW(zscore.apply(npts, x.data(), &yPtr));
    ippsNormDiff_Inf_64f(y.data(), zref.data(), npts, &error);
    EXPECT_LE(error, 1.e-14);
}

TEST(UtilitiesNormalization, minMax)
{
    MinMax minMax;
    int npts = 500;
    std::vector<double> x(npts), y(npts);
    // Generate uniform random numbers in range [0,1]
    for (auto i=0; i<npts; ++i)
    {
        x[i] = static_cast<double>(rand())/RAND_MAX;
    }
    // Fix ends points of data interval to be in range [0,1]
    x[0] = 0;
    x[1] = 1;
    // Transform from [0,1] to [-2,4] 
    std::pair<double, double> dataRange(0.0, 1.0);
    std::pair<double, double> targetRange(-2.0, 4.0);
    EXPECT_NO_THROW(minMax.initialize(dataRange, targetRange));
    double *yPtr = y.data();
    EXPECT_NO_THROW(minMax.apply(npts, x.data(), &yPtr));
    double error = 0;
    for (auto i=0; i<npts; ++i)
    {
        auto yr = -2.0 + 6.0*x[i];
        error = std::max(error, std::abs(y[i] - yr));
    }
    EXPECT_LE(error, 1.e-14);
    // Change x in [0,1] to [-5,5]
    std::vector<double> x5(npts);
    for (auto i=0; i<npts; ++i)
    {
        x5[i] = -5.0 + 10.0*x[i];
    }
    // Rescale back to [0,1]
    MinMax minMax5;
    EXPECT_NO_THROW(minMax5.initialize(npts, x5.data(),
                                       std::make_pair<double, double> (0,1)));
    minMax = minMax5;
    yPtr = y.data();
    EXPECT_NO_THROW(minMax5.apply(npts, x5.data(), &yPtr));
    error = 0;
    for (auto i=0; i<npts; ++i)
    {
         auto yr = x[i];
         error = std::max(error, std::abs(y[i] - yr));
    }
    EXPECT_LE(error, 1.e-14);
}

}
