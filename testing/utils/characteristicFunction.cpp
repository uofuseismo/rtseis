#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <chrono>
#include <ipps.h>
#include "rtseis/utilities/characteristicFunction/classicSTALTA.hpp"
#include <gtest/gtest.h>

namespace
{

std::vector<double> readTextFile(const std::string &fileName);

using namespace RTSeis::Utilities::CharacteristicFunction;

TEST(UtilitiesCharacteristicFunction, classicSTALTA)
{
    const double dt = 1/200.;
    auto nlta = 2000; //static_cast<int> (10./dt + 0.5);
    auto nsta = 1000; //static_cast<int> (5./dt + 0.5);
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
