#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/rotate/utilities.hpp"
#include <gtest/gtest.h>

namespace
{
namespace Rotate = RTSeis::Rotate;

int loadData(const std::string &fileName, // = "data/rotate_zne_rt_lqt.txt",
             std::vector<double> *verticalRef,
             std::vector<double> *northRef,
             std::vector<double> *eastRef,
             std::vector<double> *radialRef,
             std::vector<double> *transverseRef,
             std::vector<double> *lRef,
             std::vector<double> *qRef,
             std::vector<double> *tRef);

TEST(UtilitiesRotate, ne2rt)
{
    // Load the reference data
    std::vector<double> verticalRef, northRef, eastRef,
                        radialRef, transverseRef,
                        lRef, qRef, tRef;
    loadData("data/rotate_zne_rt_lqt.txt",
             &verticalRef, &northRef, &eastRef,
             &radialRef, &transverseRef,
             &lRef, &qRef, &tRef); 
    EXPECT_EQ(static_cast<int> (verticalRef.size()), 100);
    int npts = 100;
    std::vector<double> vertical(npts), north(npts), east(npts);
    std::vector<double> longitudinal(npts), radial(npts), transverse(npts),
                        q(npts), t(npts), zeroVertical(npts, 0);
    // Test north/east and radial/transverse
    double baz = 131*M_PI/180;
    double *rPtr = radial.data();
    double *tPtr = transverse.data();
    EXPECT_NO_THROW(
    Rotate::northEastToRadialTransverse(npts, baz, 
                                        northRef.data(), eastRef.data(),
                                        &rPtr, &tPtr)
    );
    double error;
    ippsNormDiff_Inf_64f(radial.data(), radialRef.data(), npts, &error); 
    EXPECT_LE(error, 1.e-7);
    ippsNormDiff_Inf_64f(transverse.data(), transverseRef.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    // Go back the other way
    double *nPtr = north.data();
    double *ePtr = east.data(); 
    EXPECT_NO_THROW(
    Rotate::radialTransverseToNorthEast(npts, baz,
                                        radial.data(), transverse.data(),
                                        &nPtr, &ePtr)
    );
    ippsNormDiff_Inf_64f(north.data(), northRef.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    ippsNormDiff_Inf_64f(east.data(), eastRef.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    // Check the LQT rotation
    double aoi = 72*M_PI/180; 
    baz = 333*M_PI/180; 
    double *lPtr = longitudinal.data(); 
    rPtr = q.data();
    tPtr = t.data();
    EXPECT_NO_THROW(
    Rotate::verticalNorthEastToLongitudinalRadialTransverse(
        npts, baz, aoi,
        verticalRef.data(), northRef.data(), eastRef.data(),
        &lPtr, &rPtr, &tPtr)
    );
    ippsNormDiff_Inf_64f(lRef.data(), longitudinal.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    // Note that the reference solution was computed with Obspy which uses 
    // (Up, South, East).  I am using (Up, North, East) so negate q
    std::vector<double> qSouth(npts); 
    for (int i=0; i<npts; ++i){qSouth[i] =-q[i];}
    ippsNormDiff_Inf_64f(qRef.data(), qSouth.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    ippsNormDiff_Inf_64f(tRef.data(), t.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    // Now go the other way
    double *vPtr = vertical.data();
    nPtr = north.data();
    ePtr = east.data();
    EXPECT_NO_THROW(
    Rotate::longitudinalRadialTransverseToVerticalNorthEast(
        npts, baz, aoi,
        longitudinal.data(), q.data(), t.data(),
        &vPtr, &nPtr, &ePtr)
    );
    ippsNormDiff_Inf_64f(verticalRef.data(), vertical.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    ippsNormDiff_Inf_64f(northRef.data(), north.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    ippsNormDiff_Inf_64f(eastRef.data(), east.data(), npts, &error);
    EXPECT_LE(error, 1.e-7);
    // Verify 2D and 3D are consistent
    std::vector<double> bazs({0, 45, 90, 135, 180, 225, 270, 335});
    for (auto bazDeg : bazs)
    {
        baz = bazDeg*M_PI/180;
        aoi = 0;
        lPtr = longitudinal.data();
        rPtr = q.data();
        tPtr = t.data();
        EXPECT_NO_THROW(
        Rotate::verticalNorthEastToLongitudinalRadialTransverse(
            npts, baz, aoi,
            zeroVertical.data(), northRef.data(), eastRef.data(), 
            &lPtr, &rPtr, &tPtr)
        );
        rPtr = radial.data();
        tPtr = transverse.data();
        EXPECT_NO_THROW(
        Rotate::northEastToRadialTransverse(npts, baz,
                                            northRef.data(), eastRef.data(),
                                            &rPtr, &tPtr)
        );
        ippsNormDiff_Inf_64f(q.data(), radial.data(), npts, &error);
        EXPECT_LE(error, 1.e-7);
        ippsNormDiff_Inf_64f(t.data(), transverse.data(), npts, &error);
        EXPECT_LE(error, 1.e-7);
    }
}

int loadData(const std::string &fileName, 
             std::vector<double> *verticalRef,
             std::vector<double> *northRef,
             std::vector<double> *eastRef,
             std::vector<double> *radialRef,
             std::vector<double> *transverseRef,
             std::vector<double> *lRef,
             std::vector<double> *qRef,
             std::vector<double> *tRef)
{
    verticalRef->resize(0);
    northRef->resize(0);
    eastRef->resize(0);
    radialRef->resize(0);
    transverseRef->resize(0);
    lRef->resize(0);
    qRef->resize(0);
    tRef->resize(0);
    std::ifstream textFile(fileName);
    if (textFile.is_open())
    {
        std::string line;
        while (std::getline(textFile, line))
        {
            double v, n, e, r, t, l1, q1, t1;
            sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf\n",
                   &v, &n, &e, &r, &t, &l1, &q1, &t1);
            verticalRef->push_back(v);
            northRef->push_back(n);
            eastRef->push_back(e);
            radialRef->push_back(r);
            transverseRef->push_back(t);
            lRef->push_back(l1);
            qRef->push_back(q1);
            tRef->push_back(t1);
        }
    }
    else
    {
        fprintf(stderr, "Can't find text file %s\n", fileName.c_str());
    }
    return -1;
}

}
