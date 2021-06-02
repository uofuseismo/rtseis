#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>
#include <ipps.h>
#include "rtseis/utilities/polarization/eigenPolarizer.hpp"
#include "rtseis/utilities/polarization/svdPolarizer.hpp"
#include "rtseis/rotate/utilities.hpp"
#include <gtest/gtest.h> 

namespace
{

using namespace RTSeis::Utilities::Polarization;
namespace Rotate = RTSeis::Rotate;

void load3C(const std::string &fileName,
            std::vector<double> &z, std::vector<double> &n, 
            std::vector<double> &e);
void loadSVDresults(const std::string &fileName,
                    std::vector<double> &klz,
                    std::vector<double> &kln,
                    std::vector<double> &kle,
                    std::vector<double> &rect,
                    std::vector<double> &incl,
                    std::vector<double> &pz,
                    std::vector<double> &pn,
                    std::vector<double> &pe,
                    std::vector<double> &sz,
                    std::vector<double> &sn,
                    std::vector<double> &se);


TEST(UtilitiesPolarization, svdPolarizer)
{
    // Load the data
    const std::string dataFile = "data/pb.b206.eh.windowed.txt"; 
    const std::string refFile = "data/svdReference.txt";
    std::vector<double> zTrace, nTrace, eTrace;
    load3C(dataFile, zTrace, nTrace, eTrace);
    auto npts = static_cast<int> (zTrace.size());
    ASSERT_EQ(npts, 5001);
    std::vector<double> klzRef, klnRef, kleRef, rectRef, incRef;
    std::vector<double> pzRef, pnRef, peRef, szRef, snRef, seRef;
    loadSVDresults(refFile, klzRef, klnRef, kleRef, rectRef, incRef,
                   pzRef, pnRef, peRef, szRef, snRef, seRef);
    ASSERT_EQ(klzRef.size(), zTrace.size());
    // Initialize the polarizer
    SVDPolarizer<double> svd;
    RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING; 
    double noise = 0.8;    // Minimize size for innovation
    double decayP = 0.99;  // 99/100 1s window at 100 samples per second
    //double decayS = 0.998; // 499/500 5s window at 100 samples per second
    EXPECT_NO_THROW(svd.initialize(decayP, noise, mode));
    // Let's fire this thing off
    EXPECT_TRUE(svd.isInitialized());
    std::vector<double> klz(npts, 0), kln(npts, 0), kle(npts, 0);
    std::vector<double> rect(npts, 0), inc(npts, 0);
    auto klzPtr = klz.data();
    auto klnPtr = kln.data();
    auto klePtr = kle.data();
    auto rectPtr = rect.data();
    auto incPtr = inc.data();
    svd.polarize(zTrace.size(),
                 zTrace.data(), nTrace.data(), eTrace.data(),
                 &klzPtr, &klnPtr, &klePtr,
                 &incPtr, &rectPtr);
    double emax;
    ippsNormDiff_Inf_64f(klzRef.data(), klz.data(), klzRef.size(), &emax);
    EXPECT_LT(emax, 1.e-6);
    ippsNormDiff_Inf_64f(klnRef.data(), kln.data(), klnRef.size(), &emax);
    EXPECT_LT(emax, 1.e-6);
    ippsNormDiff_Inf_64f(kleRef.data(), kle.data(), kleRef.size(), &emax);
    EXPECT_LT(emax, 1.e-6);
    ippsNormDiff_Inf_64f(incRef.data(), inc.data(), incRef.size(), &emax);
    EXPECT_LT(emax, 1.e-10);
    ippsNormDiff_Inf_64f(rectRef.data(), rect.data(), rectRef.size(), &emax);
    EXPECT_LT(emax, 1.e-10);
    // Try the other interfaces - KL transform
    svd.polarize(zTrace.size(),
                 zTrace.data(), nTrace.data(), eTrace.data(),
                 &klzPtr, &klnPtr, &klePtr);
    ippsNormDiff_Inf_64f(klzRef.data(), klz.data(), klzRef.size(), &emax);
    EXPECT_LT(emax, 1.e-6);
    ippsNormDiff_Inf_64f(klnRef.data(), kln.data(), klnRef.size(), &emax);
    EXPECT_LT(emax, 1.e-6);
    ippsNormDiff_Inf_64f(kleRef.data(), kle.data(), kleRef.size(), &emax);
    EXPECT_LT(emax, 1.e-6);
    // cosIncidence angle, rectilinearity only
    svd.polarize(zTrace.size(),
                 zTrace.data(), nTrace.data(), eTrace.data(),
                 &incPtr, &rectPtr);
    ippsNormDiff_Inf_64f(incRef.data(), inc.data(), incRef.size(), &emax);
    EXPECT_LT(emax, 1.e-10);
    ippsNormDiff_Inf_64f(rectRef.data(), rect.data(), rectRef.size(), &emax);
    EXPECT_LT(emax, 1.e-10);
    // Modulate the signals
    std::vector<double> pzMod(npts, 0), pnMod(npts, 0), peMod(npts, 0);
    std::vector<double> szMod(npts, 0), snMod(npts, 0), seMod(npts, 0);
    auto zModPtr = pzMod.data();
    auto nModPtr = pnMod.data();
    auto eModPtr = peMod.data();
    modulateP(npts,
              klz.data(), kln.data(), kle.data(),
              inc.data(), rect.data(),
              &zModPtr, &nModPtr, &eModPtr);  
    ippsNormDiff_Inf_64f(pzRef.data(), pzMod.data(), pzMod.size(), &emax);
    EXPECT_LT(emax, 1.e-4);
    ippsNormDiff_Inf_64f(pnRef.data(), pnMod.data(), pnMod.size(), &emax); 
    EXPECT_LT(emax, 1.e-4);
    ippsNormDiff_Inf_64f(peRef.data(), peMod.data(), peMod.size(), &emax);
    EXPECT_LT(emax, 1.e-4);

    zModPtr = szMod.data();
    nModPtr = snMod.data();
    eModPtr = seMod.data();
    modulateS(npts,
              klz.data(), kln.data(), kle.data(),
              inc.data(), rect.data(),
              &zModPtr, &nModPtr, &eModPtr);  
    ippsNormDiff_Inf_64f(szRef.data(), szMod.data(), szMod.size(), &emax);
    EXPECT_LT(emax, 1.e-4);
    ippsNormDiff_Inf_64f(snRef.data(), snMod.data(), snMod.size(), &emax); 
    EXPECT_LT(emax, 1.e-4);
    ippsNormDiff_Inf_64f(seRef.data(), seMod.data(), seMod.size(), &emax);
    EXPECT_LT(emax, 1.e-4);

    // Try real-time example. Initialize the polarizer
    SVDPolarizer<double> svdrt;
    mode = RTSeis::ProcessingMode::REAL_TIME;
    EXPECT_NO_THROW(svdrt.initialize(decayP, noise, mode));
    EXPECT_TRUE(svdrt.isInitialized());

    std::vector<double> klzRT(npts, 0), klnRT(npts, 0), kleRT(npts, 0);
    std::vector<double> rectRT(npts, 0), incRT(npts, 0);
    int packetSize = 100;
    for (int i=0; i<npts; i=i+packetSize)
    {
        auto i1 = i;
        auto i2 = std::min(npts, i + packetSize);
        auto nproc = i2 - i1;

        klzPtr = klzRT.data() + i1;
        klnPtr = klnRT.data() + i1;
        klePtr = kleRT.data() + i1;
        rectPtr = rectRT.data() + i1;
        incPtr = incRT.data() + i1;

        svdrt.polarize(nproc,
                       zTrace.data()+i1, nTrace.data()+i1, eTrace.data()+i1,
                       &klzPtr, &klnPtr, &klePtr, &incPtr, &rectPtr);
    }
    ippsNormDiff_Inf_64f(klz.data(), klzRT.data(), klzRT.size(), &emax);
    EXPECT_LT(emax, 1.e-14);
    ippsNormDiff_Inf_64f(kln.data(), klnRT.data(), klnRT.size(), &emax);
    EXPECT_LT(emax, 1.e-14);
    ippsNormDiff_Inf_64f(kle.data(), kleRT.data(), kleRT.size(), &emax);
    EXPECT_LT(emax, 1.e-14);
    ippsNormDiff_Inf_64f(inc.data(), incRT.data(), incRT.size(), &emax);
    EXPECT_LT(emax, 1.e-14);
    ippsNormDiff_Inf_64f(rect.data(), rectRT.data(), rectRT.size(), &emax);
    EXPECT_LT(emax, 1.e-14);

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

void load3C(const std::string &fileName,
            std::vector<double> &z, std::vector<double> &n,
            std::vector<double> &e)
{
    // Null out result
    z.resize(0);
    n.resize(0);
    e.resize(0);
    // Open and read file
    std::ifstream inFile;
    inFile.open(fileName);
    if (inFile.is_open())
    {
        // Space query
        std::string line;
        int npts = 0;
        while (std::getline(inFile, line))
        {
            npts = npts + 1;
        }
        inFile.close();
        inFile.open(fileName);
        // Allocate space
        z.resize(npts);
        n.resize(npts);
        e.resize(npts);
        // Load data
        int i = 0;
        while (std::getline(inFile, line))
        {
            double t;
            sscanf(line.c_str(), "%lf, %lf, %lf, %lf\n",
                    &t, &z[i], &n[i], &e[i]);
            i = i + 1;
        }
    }
    inFile.close();
}

void loadSVDresults(const std::string &fileName,
                    std::vector<double> &klz,
                    std::vector<double> &kln,
                    std::vector<double> &kle,
                    std::vector<double> &rect,
                    std::vector<double> &inc,
                    std::vector<double> &pz,
                    std::vector<double> &pn,
                    std::vector<double> &pe,
                    std::vector<double> &sz,
                    std::vector<double> &sn,
                    std::vector<double> &se)
{
    // Null out result
    klz.resize(0);
    kln.resize(0);
    kle.resize(0);
    rect.resize(0);
    inc.resize(0);
    pz.resize(0);
    pn.resize(0);
    pe.resize(0);
    sz.resize(0);
    sn.resize(0);
    se.resize(0);
    // Open and read file
    std::ifstream inFile;
    inFile.open(fileName);
    if (inFile.is_open())
    {   
        // Space query
        std::string line;
        int npts = 0;
        while (std::getline(inFile, line))
        {
            npts = npts + 1;
        }
        inFile.close();
        inFile.open(fileName);
        // Allocate space
        klz.resize(npts);
        kln.resize(npts);
        kle.resize(npts);
        rect.resize(npts);
        inc.resize(npts);
        pz.resize(npts);
        pn.resize(npts);
        pe.resize(npts);
        sz.resize(npts);
        sn.resize(npts);
        se.resize(npts);
        // Load the data
        // Load data
        int i = 0;
        while (std::getline(inFile, line))
        {
            sscanf(line.c_str(), "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
                    &klz[i], &kln[i], &kle[i], &rect[i], &inc[i],
                    &pz[i], &pn[i], &pe[i], &sz[i], &sn[i], &se[i]);
            i = i + 1;
        }
    }
    inFile.close();
}

}
