#include <fstream>
#include <cmath>
#include <vector>
#include "rtseis/amplitude/timeDomainWoodAndersonParameters.hpp"
#include "rtseis/amplitude/timeDomainWoodAnderson.hpp"
#include <gtest/gtest.h>

using namespace RTSeis::Amplitude;

namespace
{
std::vector<double> readSeismogram(const std::string &fileName)
{
    std::ifstream infl(fileName, std::ios::in);
    std::vector<double> x;
    if (infl.is_open())
    {   
        std::string line;
        x.reserve(14000);
        while (std::getline(infl, line))
        {
            double t, xi; 
            sscanf(line.c_str(), "%lf, %lf\n", &t, &xi);
            x.push_back(xi);
        }
        infl.close();
    }   
    return x;
}

}

TEST(Amplitude, TimeDomainWoodAndersonParameters)
{
    TimeDomainWoodAndersonParameters parms;
    EXPECT_TRUE(parms.isSamplingRateSupported(20));
    EXPECT_TRUE(parms.isSamplingRateSupported(40));
    EXPECT_TRUE(parms.isSamplingRateSupported(50));
    EXPECT_TRUE(parms.isSamplingRateSupported(80));
    EXPECT_TRUE(parms.isSamplingRateSupported(100));
    EXPECT_TRUE(parms.isSamplingRateSupported(200));
    EXPECT_TRUE(parms.isSamplingRateSupported(250));
    EXPECT_TRUE(parms.isSamplingRateSupported(500));
    EXPECT_FALSE(parms.isSamplingRateSupported(49));
    // Check the constants
    parms.setInputUnits(InputUnits::Velocity);
    EXPECT_EQ(parms.getInputUnits(), InputUnits::Velocity);

    const double tol = 1.e-8;

    EXPECT_NO_THROW(parms.setSamplingRate(20));
    EXPECT_NEAR(parms.getSamplingRate(), 20, 1.e-10);

    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.63382E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.14094E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.33662E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(40));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.74306E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.13422E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.31859E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(50));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.75820E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.13251E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.31203E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(80));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.77721E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12980E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.30093E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(100));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.78262E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12886E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.29694E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(200));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.79251E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12692E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.28865E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(250));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.79384E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12656E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.28695E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(500));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.79702E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12578E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.28349E+04, tol);

    parms.setInputUnits(InputUnits::Acceleration);
    EXPECT_EQ(parms.getInputUnits(), InputUnits::Acceleration);
    EXPECT_NO_THROW(parms.setSamplingRate(20));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.68405E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.13829E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.34011E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(40));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.75422E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.13357E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.31913E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(50));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.76517E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.13210E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.31235E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(80));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.77980E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12963E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.30101E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(100));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.78427E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12876E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.29700E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(200));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.79251E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12692E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.28865E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(250));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.79409E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12654E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.28695E+04, tol);
    EXPECT_NO_THROW(parms.setSamplingRate(500));
    EXPECT_NEAR(parms.getOptimizedDampingConstant(),         0.79711E+00, tol);
    EXPECT_NEAR(parms.getOptimizedNaturalAngularFrequency(), 0.12578E+01, tol);
    EXPECT_NEAR(parms.getOptimizedGain(),                    0.28349E+04, tol);

    const double gain = 321168.428435*100;
    EXPECT_NO_THROW(parms.setSimpleResponse(gain));
    EXPECT_NEAR(parms.getSimpleResponse(), gain, 1.e-8);

    const WoodAndersonGain waGain = WoodAndersonGain::WA_2080;
    parms.setWoodAndersonGain(waGain);
    EXPECT_EQ(parms.getWoodAndersonGain(), waGain);

    const DetrendType detrendType = DetrendType::RemoveTrend;
    parms.setDetrendType(detrendType);
    EXPECT_EQ(parms.getDetrendType(), detrendType);

    const WindowType taperType = WindowType::None;
    const double pct = 5;
    EXPECT_NO_THROW(parms.setTaper(pct, taperType));
    EXPECT_EQ(parms.getWindowType(), taperType);
    EXPECT_NEAR(parms.getTaperPercentage(), pct, 1.e-10);

    // Test copy
    TimeDomainWoodAndersonParameters parmsCopy(parms);
    EXPECT_EQ(parmsCopy.getInputUnits(), InputUnits::Acceleration);
    EXPECT_NEAR(parmsCopy.getSamplingRate(), 500, 1.e-10);
    EXPECT_NEAR(parmsCopy.getOptimizedDampingConstant(),
                0.79711E+00, tol);
    EXPECT_NEAR(parmsCopy.getOptimizedNaturalAngularFrequency(),
                0.12578E+01, tol);
    EXPECT_NEAR(parmsCopy.getOptimizedGain(),
                0.28349E+04, tol);
    EXPECT_NEAR(parmsCopy.getSimpleResponse(), gain, 1.e-8);
    EXPECT_EQ(parmsCopy.getWoodAndersonGain(), waGain);
    EXPECT_EQ(parmsCopy.getDetrendType(), detrendType);
    EXPECT_NEAR(parmsCopy.getTaperPercentage(), pct, 1.e-10);
}

TEST(Amplitude, TimeDomainWoodAndersonVelocity)
{
    double dt = 0.01;
    double df = std::round(1/dt);
    double gain = 1274800764.8712056/100;
    double pct = 0.15;
    auto x = readSeismogram("data/UU.SPU.HHN.01.txt");
    ASSERT_EQ(static_cast<int> (x.size()), 7201);

    TimeDomainWoodAndersonParameters parameters;
    parameters.setInputUnits(InputUnits::Velocity);
    parameters.setSamplingRate(df);
    parameters.setSimpleResponse(gain);
    parameters.setDetrendType(DetrendType::RemoveMean);
    parameters.setTaper(pct, WindowType::Sine);
    parameters.setWoodAndersonGain(WoodAndersonGain::WA_2800);

    TimeDomainWoodAnderson<RTSeis::ProcessingMode::POST, double> filter;
    EXPECT_NO_THROW(filter.initialize(parameters));
    EXPECT_TRUE(filter.isInitialized());
    EXPECT_TRUE(filter.isVelocityFilter());
    std::vector<double> y(x.size(), 0);
    double *yPtr = y.data();
    filter.apply(x.size(), x.data(), &yPtr);

    std::ofstream ofl("waVel.txt");
    for (int i = 0; i < y.size(); ++i)
    {
        ofl << dt*i << " " << y[i] << std::endl;
    }
    ofl.close();
}

