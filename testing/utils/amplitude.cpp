#include <fstream>
#include <cmath>
#include <vector>
#include "rtseis/amplitude/timeDomainWoodAndersonParameters.hpp"
#include "rtseis/amplitude/timeDomainWoodAnderson.hpp"
#include "rtseis/amplitude/tauPParameters.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
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

    const double q = 0.9;
    EXPECT_NO_THROW(parms.setHighPassRCFilter(q, HighPassFilter::Yes));
    EXPECT_NEAR(parms.getHighPassFilterQ(), q, 1.e-10);
    EXPECT_EQ(parms.getHighPassFilter(), HighPassFilter::Yes);

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
    EXPECT_NEAR(parmsCopy.getHighPassFilterQ(), q, 1.e-10);
    EXPECT_EQ(parmsCopy.getHighPassFilter(), HighPassFilter::Yes);
}

TEST(Amplitude, TimeDomainWoodAndersonVelocity)
{
    double q = 0.998;
    double dt = 0.01;
    double df = std::round(1/dt);
    double gain = 1274800764.8712056/100;
    double accGain = 321168.428435/100;
    double pct = 5;
    auto x = readSeismogram("data/UU.SPU.HHN.01.txt");
    ASSERT_EQ(static_cast<int> (x.size()), 7201);
    auto xAcc = readSeismogram("data/UU.TMU.ENN.01.txt");
    ASSERT_EQ(static_cast<int> (xAcc.size()), 13489L);

    TimeDomainWoodAndersonParameters parameters;
    parameters.setInputUnits(InputUnits::Velocity);
    parameters.setSamplingRate(df);
    parameters.setSimpleResponse(gain);
    parameters.setDetrendType(DetrendType::RemoveMean);
    parameters.setTaper(pct, WindowType::None);
    parameters.setWoodAndersonGain(WoodAndersonGain::WA_2800);

    TimeDomainWoodAnderson<RTSeis::ProcessingMode::POST, double> filter;
    EXPECT_NO_THROW(filter.initialize(parameters));
    EXPECT_TRUE(filter.isInitialized());
    EXPECT_TRUE(filter.isVelocityFilter());
    std::vector<double> y(x.size(), 0);
    double *yPtr = y.data();
    filter.apply(x.size(), x.data(), &yPtr);

    std::ofstream ofl("waVel.txt");
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {
        ofl << dt*i << " " << y[i] << std::endl;
    }
    ofl.close();

    TimeDomainWoodAndersonParameters accParameters;
    accParameters.setInputUnits(InputUnits::Acceleration);
    accParameters.setSamplingRate(df);
    accParameters.setSimpleResponse(accGain);
    accParameters.setDetrendType(DetrendType::RemoveMean);
    accParameters.setTaper(pct, WindowType::Sine);
    accParameters.setWoodAndersonGain(WoodAndersonGain::WA_2800);
    accParameters.setHighPassRCFilter(q, HighPassFilter::Yes);
    TimeDomainWoodAnderson<RTSeis::ProcessingMode::POST, double> accFilter;
    accFilter.initialize(accParameters);
    std::vector<double> yAcc(xAcc.size());
    yPtr = yAcc.data();
    accFilter.apply(xAcc.size(), xAcc.data(), &yPtr);

    ofl.open("waAcc.txt");
    for (int i =0; i < static_cast<int> (y.size()); ++i)
    {
        ofl << dt*i << " " << yAcc[i] << std::endl;
    }
    ofl.close();
}

TEST(Amplitude, TauPParameters)
{
    TauPParameters parameters;
    const double samplingRate{50};
    const double alphaDefault = 1. - 1./samplingRate;
    const double alpha = 0.9;
    const double gain = 450;
    const double q = 0.8;
    const InputUnits units = InputUnits::Acceleration;
    const DetrendType detrendType = DetrendType::RemoveTrend;
    const WindowType taperType = WindowType::None;
    const double pct = 8;
    EXPECT_NO_THROW(parameters.setSamplingRate(samplingRate));
    parameters.setInputUnits(units);
    EXPECT_NO_THROW(parameters.setSimpleResponse(gain));
    EXPECT_NO_THROW(parameters.setSmoothingParameter(alphaDefault));
    EXPECT_NEAR(parameters.getSmoothingParameter(), alphaDefault, 1.e-10); 
    EXPECT_NO_THROW(parameters.setSmoothingParameter(alpha));
    EXPECT_NO_THROW(parameters.setFilterConstantQ(q));
    parameters.setDetrendType(detrendType);
    EXPECT_NO_THROW(parameters.setTaper(pct, taperType));
    

    TauPParameters copy(parameters);
    EXPECT_EQ(copy.getInputUnits(), units);
    EXPECT_EQ(copy.getDetrendType(), detrendType);
    EXPECT_NEAR(copy.getSamplingRate(), samplingRate, 1.e-10);
    EXPECT_NEAR(copy.getSmoothingParameter(), alpha, 1.e-10);
    EXPECT_NEAR(copy.getSimpleResponse(), gain, 1.e-10);
    EXPECT_NEAR(copy.getFilterConstantQ(), q, 1.e-10);
    EXPECT_EQ(copy.getWindowType(), taperType);
    EXPECT_NEAR(copy.getTaperPercentage(), pct, 1.e-10);

    auto ba = copy.getVelocityFilter();
    std::vector<double> bRef{0.02785977,  0.05571953,  0.02785977};
    std::vector<double> aRef{1,          -1.47548044,  0.58691951};
    auto b = ba.getNumeratorCoefficients();
    auto a = ba.getDenominatorCoefficients();
    EXPECT_EQ(b.size(), bRef.size());
    double resMax = 0;
    for (int i = 0; i < static_cast<int> (b.size()); ++i)
    { 
        resMax = std::max(resMax, std::abs(b[i] - bRef[i]));
        resMax = std::max(resMax, std::abs(a[i] - aRef[i]));
    }
    EXPECT_NEAR(resMax, 0, 1.e-6);

/*
    ba = copy.getAccelerationFilter();
    std::vector<double> bRef2{0.9877613892768257,
                             -3.9510455571073027,
                              5.926568335660954,
                             -3.9510455571073027,
                              0.9877613892768257};
    std::vector<double> aRef2{1.0,
                             -3.97537191256092,
                              5.926418555965423,
                             -3.926719197756784,
                              0.975672562146085};
*/
/* 
    std::vector<double> bRef2{0.98776139, -1.97552278,  0.98776139,
                              1.        , -2.        ,  1.};
    std::vector<double> aRef2{1.00000000e+00, -1.9826478,  0.99281262,
                              1.00000000e+00, -1.99272411, 0.99281262};
*/
/*
    b = ba.getNumeratorCoefficients();
    a = ba.getDenominatorCoefficients();
    EXPECT_EQ(b.size(), bRef2.size());
    resMax = 0;
    for (int i = 0; i < static_cast<int> (b.size()); ++i)
    {
        resMax = std::max(resMax, std::abs(b[i] - bRef2[i]));
        resMax = std::max(resMax, std::abs(a[i] - aRef2[i]));
    }
*/
}

TEST(Amplitude, TauP)
{

}

}
