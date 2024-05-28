#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <string>
#include <vector>
#include "rtseis/trigger/waterLevel.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis::Trigger;

TEST(UtilitiesTrigger, waterLevel)
{
    WaterLevel<RTSeis::ProcessingMode::POST_PROCESSING, double> trigger;
    double triggerOn = 0.8;
    double triggerOff = 0.2;
    // Create an oscillatory function
    double dt = 0.01;
    double freq = 1;  // 1 per second
    double tlen = 10; // Should be 10
    auto len = static_cast<int> (tlen/dt) + 1;
    std::vector<double> x(len);
    std::vector<bool> refIsOn(len, false);
    for (int i=0; i<len; ++i)
    {
        auto t = dt*i;
        x[i] = std::sin(2*M_PI*freq*t);
        if (x[i] > triggerOn && x[i] > triggerOff){refIsOn[i] = true;}
    }
    EXPECT_NO_THROW(trigger.initialize(triggerOn, triggerOff));
    EXPECT_TRUE(trigger.isInitialized());
    EXPECT_NO_THROW(trigger.apply(x.size(), x.data()));
    EXPECT_EQ(trigger.getNumberOfWindows(), 10);
//printf("%d\n", trigger.getNumberOfWindows());
    auto triggers = trigger.getWindows();
    // Compute the reference soluion 
    std::vector<int> isOn(len, false); 
    for (int i=0; i<static_cast<int> (triggers.size()); ++i)
    {
        auto j1 = std::max(0, triggers[i].first - 1);
        auto j2 = triggers[i].second - 1;
        EXPECT_TRUE(x[j1] < triggerOn  && x[j1+1] > triggerOn);
        EXPECT_TRUE(x[j2] > triggerOff && x[j2+1] < triggerOff);
    }
/*
 FILE *fout = fopen("trigger.txt", "w");
 for (int i=0; i<static_cast<int> (isOn.size()); ++i)
 {
   fprintf(fout, "%lf, %d, %d, %lf\n", i*dt, (int) refIsOn[i], (int) isOn[i], x[i]);
 }
fclose(fout);
*/
}

}
