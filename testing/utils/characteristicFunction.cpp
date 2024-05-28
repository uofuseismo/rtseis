#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <numeric>
#include <cmath>
#include <vector>
#include <chrono>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/utilities/characteristicFunction/classicSTALTA.hpp"
#include "rtseis/utilities/characteristicFunction/carlSTALTA.hpp"
#include <gtest/gtest.h>

namespace
{

std::vector<double> readTextFile(const std::string &fileName);
std::vector<double> computeCarlSTALTA(const int n,
                                      const double x[],
                                      const int nsta,
                                      const int nlta,
                                      const double ratio,
                                      const double quiet);

using namespace RTSeis::Utilities::CharacteristicFunction;

TEST(UtilitiesCharacteristicFunction, classicSTALTA)
{
    int nlta = 2000; //static_cast<int> (10./dt + 0.5);
    int nsta = 1000; //static_cast<int> (5./dt + 0.5);
    auto x = readTextFile("data/gse2.txt");
    auto yRef = readTextFile("data/classicSTALTA_ref.txt");
    ASSERT_TRUE(x.size() > 0);
    ASSERT_TRUE(x.size() == yRef.size());

    PostProcessing::ClassicSTALTA<double> stalta;
    EXPECT_NO_THROW(stalta.initialize(nsta, nlta));
    EXPECT_TRUE(stalta.isInitialized());
    std::vector<double> y(x.size());
    auto yPtr = y.data();
    EXPECT_NO_THROW(stalta.apply(x.size(), x.data(), &yPtr));
    double error;
    ippsNormDiff_Inf_64f(yRef.data(), y.data(), yRef.size(), &error);
    EXPECT_LT(error, 1.e-8);
    // Ensure resetting is okay
    EXPECT_NO_THROW(stalta.apply(x.size(), x.data(), &yPtr);
    ippsNormDiff_Inf_64f(yRef.data(), y.data(), yRef.size(), &error));
    EXPECT_LT(error, 1.e-8);
    // Test copy c'tor
    PostProcessing::ClassicSTALTA<double> staltaCopy(stalta);
    EXPECT_TRUE(staltaCopy.isInitialized());
    EXPECT_NO_THROW(staltaCopy.apply(x.size(), x.data(), &yPtr));
    EXPECT_LT(error, 1.e-8);
    
    // Test the real-time module
}

TEST(UtilitiesCharacteristicFunction, carlSTALTA)
{
    const double dt = 1./200;;
    double ratio = 0.2;
    double quiet = 2;
    auto nlta = static_cast<int> (10./dt + 0.5);
    auto nsta = static_cast<int> (1./dt + 0.5);
    auto x = readTextFile("data/gse2.txt");
    //auto yRef = readTextFile("data/classicSTALTA_ref.txt");
    auto yRef = computeCarlSTALTA(x.size(), x.data(), nsta, nlta, ratio, quiet);
    ASSERT_TRUE(x.size() > 0);
    //ASSERT_TRUE(x.size() == yRef.size());

    CarlSTALTA<double, RTSeis::ProcessingMode::POST_PROCESSING> stalta;
    EXPECT_NO_THROW(stalta.initialize(nsta, nlta, ratio, quiet));
    EXPECT_TRUE(stalta.isInitialized());
    std::vector<double> y(x.size());
    auto yPtr = y.data();
    EXPECT_NO_THROW(stalta.apply(x.size(), x.data(), &yPtr));

    double error;
    ippsNormDiff_Inf_64f(yRef.data(), y.data(), yRef.size(), &error);
    EXPECT_LT(error, 1.e-8);

    // Test the real-time module

/*
    std::ofstream outfl("carlSTALTA.txt");
    for (int i=0; i<static_cast<int> (x.size()); ++i)
    {
        outfl << dt*i << " " << y[i] << " " << yRef[i] << " " << x[i] << std::endl;
    }
    outfl.close(); 
*/
}

std::vector<double> computeCarlSTALTA(const int n,
                                      const double x[],
                                      const int nsta,
                                      const int nlta,
                                      const double ratio, 
                                      const double quiet)
{
    const double zero = 0;
    auto npad = std::max(nsta, nlta);
    std::vector<double> xwork(npad + n, 0);
    std::vector<double> sta(npad + n, 0);
    std::vector<double> lta(npad + n, 0); 
    std::copy(x, x + n, xwork.data() + npad);
    // Compute STA and LTA signals
    for (int i=0; i<n; ++i)
    {
        auto ix = i + npad;
        sta[ix] = std::accumulate(xwork.data() + ix - nsta + 1,
                                  xwork.data() + ix + 1, zero);
        sta[ix] = sta[ix]/static_cast<double> (nsta);
        lta[ix] = std::accumulate(xwork.data() + ix - nlta + 1,
                                  xwork.data() + ix + 1, zero);
        lta[ix] = lta[ix]/static_cast<double> (nlta);
    }
    // Compute rectified average signals
    std::vector<double> result(n, 0);
    for (int i=0; i<n; ++i)
    {
         auto ix = i + npad;
         double star = 0;
         #pragma omp simd reduction(+:star)
         for (int j=0; j<nsta; ++j)
         {
             star = star + std::abs(xwork[ix-j] - lta[ix-j]);
         }
         star = star/static_cast<double> (nsta);
         double ltar = 0;
         #pragma omp simd reduction(+:ltar)
         for (int j=0; j<nlta; ++j)
         {
             ltar = ltar + std::abs(xwork[ix-j] - lta[ix-j]);
         }
         ltar = ltar/static_cast<double> (nlta);
         // Compute the result
         result[i] = star - ratio*ltar - std::abs(sta[ix] - lta[ix]) - quiet;
     }
     return result;
}

std::vector<double> readTextFile(const std::string &fileName)
{
    std::vector<double> x;
    char line[64];
    FILE *fl = fopen(fileName.c_str(), "r");
    int npts = 0; 
    while (fscanf(fl, "%s", line) != EOF) 
    {    
        npts = npts + 1; 
    }    
    rewind(fl);
    if (npts < 1) 
    {    
        fprintf(stderr, "%s", "No data points in file\n");
        return x;
    }
    x.resize(npts);
    for (int i=0; i<npts; i++) 
    {
        memset(line, 0, 64*sizeof(char));
        fgets(line, 64, fl);
        sscanf(line, "%lf\n", &x[i]);
    }
    fclose(fl);
    return x;
}

}
