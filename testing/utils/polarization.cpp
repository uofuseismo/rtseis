#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <cmath>
#include <ipps.h>
#include "rtseis/utilities/polarization/eigenPolarizer.hpp"
#include "rtseis/utilities/polarization/svdPolarizer.hpp"
#include "rtseis/utilities/rotate/utilities.hpp"
#include <gtest/gtest.h> 

namespace
{

using namespace RTSeis::Utilities::Polarization;
namespace Rotate = RTSeis::Utilities::Rotate;

TEST(UtilitiesPolarization, svdPolarizer)
{
    SVDPolarizer<double> svd;
}

TEST(UtilitiesPolarization, eigenPolarizer)
{
    EigenPolarizer<double> eigen;
    // Make an (L,Q,T) P-wave trace and rotate it to (Z,N,E)
    double freq = 1.0;
    double omega = 2.0*M_PI*freq;
    double T = 1./freq;
    double dt = 1./40.;
    int nt = static_cast<int> (T/dt + 0.5) + 1; 
    std::vector<double> longitudinal(nt), radial(nt), transverse(nt);
    for (int i=0; i<nt; ++i)
    {
        longitudinal[i] = std::sin(omega*(i*dt));
        radial[i] = 0;
        transverse[i] = 0;
    }
    // Initialize class
    EXPECT_NO_THROW(eigen.initialize(nt));
    // Vertical incidence will result in undefined azimuth.  This is
    // because the azimuth for a plane wave the station from directly below 
    // will have an azimuth that is valid for all angles.
    std::vector<double> aois({0.1, 20, 45, 60, 90});
    std::vector<double> bazs({0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 359});
    for (auto &aoi : aois)
    {
        for (auto &baz : bazs)
        {
            double az = baz + 180;
            if (az > 360){az = az - 360;}
            std::vector<double> vertical(nt), east(nt), north(nt);
            // Rotate LQT to ZNE
            double *vPtr = vertical.data();
            double *nPtr = north.data();
            double *ePtr = east.data();
            Rotate::longitudinalRadialTransverseToVerticalNorthEast(
                        nt, baz*M_PI/180, aoi*M_PI/180,
                        longitudinal.data(), radial.data(), transverse.data(),
                        &vPtr, &nPtr, &ePtr); 
            // Now, recover the azimuth and incidence angle
            EXPECT_TRUE(eigen.isInitialized());
            EXPECT_EQ(eigen.getNumberOfSamples(), nt);
            EXPECT_NO_THROW(eigen.setSignals(nt, vertical.data(),
                            north.data(), east.data()));
            EXPECT_NEAR(eigen.getRectilinearity(), 1.0, 1.e-6);
            bool lwantRadians = false;
            EXPECT_NEAR(eigen.getAzimuth(lwantRadians), az, 1.e-6);
            EXPECT_NEAR(eigen.getBackAzimuth(lwantRadians), baz, 1.e-6);
            EXPECT_NEAR(eigen.getIncidenceAngle(lwantRadians), aoi, 1.e-6);
            lwantRadians = true;
            EXPECT_NEAR(eigen.getAzimuth(lwantRadians), az*M_PI/180, 1.e-6);
            EXPECT_NEAR(eigen.getBackAzimuth(lwantRadians), baz*M_PI/180, 1.e-6);
            EXPECT_NEAR(eigen.getIncidenceAngle(lwantRadians), aoi*M_PI/180, 1.e-6);
            //printf("%lf %lf %lf %lf\n", eigen.getAzimuth(false), eigen.getBackAzimuth(false), eigen.getIncidenceAngle(false), eigen.getRectilinearity());
        }
    }
}

}
