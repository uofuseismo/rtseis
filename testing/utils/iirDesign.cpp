#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterDesign/analogPrototype.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis;
using namespace RTSeis::FilterRepresentations;
using namespace RTSeis::FilterDesign;

///  IIR analog prototype design. 
//int rtseis_test_utils_design_iir_ap(void)
TEST(UtilitiesDesignIIR, analogPrototype)
{
    double k;
    ZPK zpk, zpkRef;
    std::vector<std::complex<double>> pref;
    std::vector<std::complex<double>> zref;
    // Test butterworth order 1
    EXPECT_NO_THROW(zpk = IIR::AnalogPrototype::butter(1));
    k = 1;
    pref.push_back(std::complex<double> (-1, 0));
    zpkRef = ZPK(zref, pref, k);
    EXPECT_EQ(zpkRef, zpk);
    zpkRef.clear();
    // Test butterworth order 4
    EXPECT_NO_THROW(zpk = IIR::AnalogPrototype::butter(5));
    pref.clear();
    zref.clear();
    pref.resize(5);
    k = 1;
    pref[0] = std::complex<double> (-0.3090169943749474,  0.95105651629515353);
    pref[1] = std::complex<double> (-0.80901699437494745, 0.58778525229247303);
    pref[2] = std::complex<double> (-1, 0);
    pref[3] = std::complex<double> (-0.80901699437494745,-0.58778525229247303);
    pref[4] = std::complex<double> (-0.30901699437494751,-0.95105651629515353);
    zpkRef = ZPK(zref, pref, k);
    EXPECT_EQ(zpkRef, zpk);
    zpkRef.clear();
    // Test order 1 cheby1
    EXPECT_NO_THROW(zpk = IIR::AnalogPrototype::cheb1ap(1, 2.2));
    zref.clear();
    pref.clear();
    pref.resize(1);
    k = 1.2313003041963828;
    pref[0] = std::complex<double> (-1.2313003041963828, 0);
    zpkRef = ZPK(zref, pref, k);
    EXPECT_EQ(zpkRef, zpk);
    // Test order 6 cheby1
    EXPECT_NO_THROW(zpk = IIR::AnalogPrototype::cheb1ap(6, 0.994));
    zref.clear();
    pref.clear();
    pref.resize(6);
    k = 0.061620501119488615;
    pref[0] = std::complex<double> (-0.062314231644038744,+0.9935274525619241); 
    pref[1] = std::complex<double> (-0.17024564688613014,+0.72731257398980598);
    pref[2] = std::complex<double> (-0.23255987853016891,+0.26621487857211812);
    pref[3] = std::complex<double> (-0.23255987853016891,-0.26621487857211795);
    pref[4] = std::complex<double> (-0.17024564688613017,-0.72731257398980587);
    pref[5] = std::complex<double> (-0.062314231644038813,-0.99352745256192398);
    zpkRef = ZPK(zref, pref, k); 
    EXPECT_EQ(zpkRef, zpk);
    // Test order 2 cheby2 
    EXPECT_NO_THROW(zpk = IIR::AnalogPrototype::cheb2ap(1, 1.1));
    zref.clear();
    pref.clear();
    pref.resize(1);
    k = 1.862583192806328;
    pref[0] = std::complex<double> (-1.862583192806328,0);
    zpkRef = ZPK(zref, pref, k);
    EXPECT_EQ(zpkRef, zpk);
    // Test order 6 cheby2 
    EXPECT_NO_THROW(zpk = IIR::AnalogPrototype::cheb2ap(6, 1.2));
    k = 0.8709635899560811;
    zref.resize(6);
    pref.resize(6);
    zref[0] = std::complex<double> (0, -1.035276180410083);
    zref[1] = std::complex<double> (0, -1.4142135623730951);
    zref[2] = std::complex<double> (0, -3.8637033051562737);
    zref[3] = std::complex<double> (0,  3.8637033051562737);
    zref[4] = std::complex<double> (0,  1.4142135623730951);
    zref[5] = std::complex<double> (0,  1.035276180410083); 
    pref[0] = std::complex<double> (-0.024686186266327684,-1.0305393933278832);
    pref[1] = std::complex<double> (-0.12492582633346083,-1.3973824335027194);
    pref[2] = std::complex<double> (-1.1553327165440157,-3.462761441343769);
    pref[3] = std::complex<double> (-1.1553327165440157,+3.462761441343769);
    pref[4] = std::complex<double> (-0.12492582633346083,+1.3973824335027194);
    pref[5] = std::complex<double> (-0.024686186266327684,+1.0305393933278832);
    zpkRef = ZPK(zref, pref, k);
    EXPECT_EQ(zpkRef, zpk);
}
//============================================================================//
//int rtseis_test_utils_design_zpk2tf(void)
TEST(UtilitiesDesignIIR, zpk2tf)
{
    std::vector<std::complex<double>> zref;
    std::vector<std::complex<double>> pref;
    double kref = 0.8709635899560811;
    zref.resize(6);
    pref.resize(6);
    zref[0] = std::complex<double> (0, -1.035276180410083);
    zref[1] = std::complex<double> (0, -1.4142135623730951);
    zref[2] = std::complex<double> (0, -3.8637033051562737);
    zref[3] = std::complex<double> (0,  3.8637033051562737);
    zref[4] = std::complex<double> (0,  1.4142135623730951);
    zref[5] = std::complex<double> (0,  1.035276180410083);
    pref[0] = std::complex<double> (-0.024686186266327684,-1.0305393933278832);
    pref[1] = std::complex<double> (-0.12492582633346083,-1.3973824335027194);
    pref[2] = std::complex<double> (-1.1553327165440157,-3.462761441343769);
    pref[3] = std::complex<double> (-1.1553327165440157,+3.462761441343769);
    pref[4] = std::complex<double> (-0.12492582633346083,+1.3973824335027194);
    pref[5] = std::complex<double> (-0.024686186266327684,+1.0305393933278832);
    ZPK zpkref, zpk;
    EXPECT_NO_THROW(zpkref = ZPK(zref, pref, kref));

    BA ba;
    EXPECT_NO_THROW(ba = IIR::zpk2tf(zpkref));
    EXPECT_NO_THROW(zpk = IIR::tf2zpk(ba));
    // Poles aren't gauranteed to be in any order so need to look
    std::vector<std::complex<double>> z = zpk.getZeros();
    std::vector<std::complex<double>> p = zpk.getPoles();
    double k = zpk.getGain();
    EXPECT_NEAR(k, kref, 1.e-13); //std::abs(k - kref) > 1.e-13)
    EXPECT_EQ(z.size(), zref.size());
    for (auto i=0; i<static_cast<int> (z.size()); i++)
    {
        bool lfound = false;
        for (auto j=0; j<static_cast<int> (zref.size()); j++)
        {
            if (std::abs(zref[j] - z[i]) < 1.e-13){lfound = true;}
        }
        EXPECT_TRUE(lfound);
    }
    assert(p.size() == pref.size());
    for (auto i=0; i<static_cast<int> (p.size()); i++)
    {
        bool lfound = false;
        for (auto j=0; j<static_cast<int> (pref.size()); j++)
        {
            if (std::abs(pref[j] - p[i]) < 1.e-13){lfound = true;}
        }
        EXPECT_TRUE(lfound);
    }
}
TEST(UtilitiesDesignIIR, tf2sos)
{
    std::vector<double> b{0.9877613892768257,
                         -3.9510455571073027,
                          5.926568335660954,
                         -3.9510455571073027,
                          0.9877613892768257};
    std::vector<double> a{1.0,
                         -3.97537191256092,
                          5.926418555965423,
                         -3.926719197756784,
                          0.975672562146085};
    BA ba;
    ba.setNumeratorCoefficients(b);
    ba.setDenominatorCoefficients(a);
    auto sos = IIR::tf2sos(ba,  SOSPairing::NEAREST);

    std::vector<double> bRef{//0.98776139, -1.97552278,  0.98776139,
                             0.9877613892768257, 
                            -1.9755228135715934,
                             0.987761379675909,
                             1.        , -1.9999999645481794,  1.000000009719877};
    std::vector<double> aRef{1.0, -1.982647799569223, 0.9827358586551028,
                             1.0, -1.9927241129916946,0.9928126195387987};
    SOS sosRef;
    sosRef.setSecondOrderSections(2, bRef, aRef);

    sosRef.setEqualityTolerance(1.e-7);
    sos.setEqualityTolerance(1.e-7);
    EXPECT_EQ(sosRef, sos);
}

//============================================================================//
//int rtseis_test_utils_design_iir(void)
TEST(UtilitiesDesignIIR, butterworth)
{
    int n;
    BA ba;
    IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL;
    double rp = 5, rs = 60;
    double Wn[2];
    //------------------------------------------------------------------------//
    //                         Butterworth filter designs                     //
    //------------------------------------------------------------------------//
    std::vector<double> bref_blp({0.000000025260,   0.000000227339,
                                  0.000000909357,   0.000002121833,
                                  0.000003182750,   0.000003182750,
                                  0.000002121833,   0.000000909357,
                                  0.000000227339,   0.000000025260});
    std::vector<double> aref_blp({1.000000000000,  -7.191438395125,
                                  23.136179228821, -43.681913446689,
                                  53.315329241295, -43.608590547581,
                                  23.895557703800,  -8.456015196801,
                                  1.753099824111,   -0.162195478751});
    BA butterlp_ref(bref_blp, aref_blp);
    n = 9, Wn[0] = 0.1; Wn[1] = 0; rp = 5; rs = 60;
    EXPECT_NO_THROW(
       ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::LOWPASS,
                                    IIRPrototype::BUTTERWORTH,
                                    ldigital);
    );
    EXPECT_EQ(ba, butterlp_ref);

    std::vector<double> bref_bhp({0.402734998170,  -3.624614983529,
                                  14.498459934115, -33.829739846269,
                                  50.744609769404, -50.744609769404, 
                                  33.829739846269, -14.498459934115,
                                  3.624614983529,  -0.402734998170}); 
    std::vector<double> aref_bhp({1.000000000000,  -7.191438395125,
                                  23.136179228821, -43.681913446689,
                                  53.315329241295, -43.608590547581, 
                                  23.895557703800,  -8.456015196801,
                                  1.753099824111,  -0.162195478751});
    BA butterhp_ref(bref_bhp, aref_bhp);
    n = 9, Wn[0] = 0.1; Wn[1] = 0; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs,
                                    Bandtype::HIGHPASS,
                                    IIRPrototype::BUTTERWORTH,
                                    ldigital);
    );
    EXPECT_EQ(ba, butterhp_ref);

    std::vector<double> bref_bbp({0.001065394524,   0.000000000000,
                                 -0.009588550712,   0.000000000000,
                                  0.038354202850,  -0.000000000000,
                                 -0.089493139982,   0.000000000000,
                                  0.134239709973,  -0.000000000000,
                                 -0.134239709973,  -0.000000000000,
                                  0.089493139982,  -0.000000000000,
                                 -0.038354202850,   0.000000000000,
                                  0.009588550712,   0.000000000000,
                                 -0.001065394524});
    std::vector<double> aref_bbp({1.000000000000,  -4.122017284406,
                                  9.504415355480,  -15.991234305338,
                                  22.198678051760, -26.228218162255,
                                  26.787904695666, -23.899415459352,
                                  18.869452192123, -13.188721760666,
                                  8.151696207859,  -4.434500672760,
                                  2.117697034965,  -0.874531051117,
                                  0.307669850093,  -0.089310128165,
                                  0.020735426453,  -0.003439408018,
                                  0.000355580604});
    BA butterbp_ref(bref_bbp, aref_bbp);
    n = 9, Wn[0] = 0.2; Wn[1] = 0.6; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs,
                                    Bandtype::BANDPASS,
                                    IIRPrototype::BUTTERWORTH,
                                    ldigital);
    );
    ba.setEqualityTolerance(1.e-10);
    EXPECT_EQ(ba , butterbp_ref);

    std::vector<double> bref_bbs({0.018886917953,  -0.129854892873,
                                  0.566783505349,  -1.746140555257,
                                  4.268033050318,  -8.498908785998,
                                  14.287136169364, -20.461288100967,
                                  25.338462554529, -27.159003561885,
                                  25.338462554529, -20.461288100967,
                                  14.287136169364,  -8.498908785998,
                                   4.268033050318,  -1.746140555257,
                                   0.566783505349,  -0.129854892873,
                                  0.018886917953});
    std::vector<double> aref_bbs({1.000000000000,  -4.122017284406,
                                   9.504415355480, -15.991234305338,
                                  22.198678051760, -26.228218162255,
                                  26.787904695666, -23.899415459351, 
                                  18.869452192123, -13.188721760666,
                                  8.151696207860,  -4.434500672759,
                                  2.117697034965,  -0.874531051117,
                                  0.307669850093,  -0.089310128165,
                                  0.020735426453,  -0.003439408018,
                                  0.000355580604});
    BA butterbs_ref(bref_bbs, aref_bbs);
    n = 9, Wn[0] = 0.2; Wn[1] = 0.6; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::BANDSTOP,
                                    IIRPrototype::BUTTERWORTH,
                                    ldigital);
    );
    ba.setEqualityTolerance(1.e-10); 
    EXPECT_EQ(ba, butterbs_ref);
}

TEST(UtilitiesDesignIIR, bessel)
{
    int n;
    BA ba; 
    IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL;
    double rp = 5, rs = 60; 
    double Wn[2];
    //------------------------------------------------------------------------//
    //                           Bessel filter designs                        //
    //------------------------------------------------------------------------//
    std::vector<double> bref_belp({
         0.000000022600,   0.000000203401,   0.000000813606,
         0.000001898413,   0.000002847620,   0.000002847620,
         0.000001898413,   0.000000813606,   0.000000203401,
         0.000000022600});
    std::vector<double> aref_belp({
         1.000000000000,  -6.995866412951,  21.889109813443,
       -40.190364280954,  47.708322046868, -37.959619208469,
        20.239343756074,  -6.971442215827,   1.407371385616,
        -0.126843312520});
    BA bessellp_ref(bref_belp, aref_belp);
    n = 9, Wn[0] = 0.1; Wn[1] = 0; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::LOWPASS,
                                    IIRPrototype::BESSEL,
                                    ldigital)
    );
    EXPECT_EQ(ba, bessellp_ref);

    std::vector<double> bref_behp({
         0.1237642230273, -1.1138780072456,  4.4555120289825, 
       -10.3961947342925, 15.5942921014388, -15.5942921014388, 
        10.3961947342925, -4.4555120289825,  1.1138780072456, 
       -0.1237642230273});
    std::vector<double> aref_behp({
        1.0000000000000, -5.1068260795845, 11.8775397999909, 
      -16.4395829343521, 14.8781830025086, -9.1106255014978, 
        3.7683662953178, -1.0138777253587, 0.1608308466516, 
       -0.0114500047115});
    BA besselhp_ref(bref_behp, aref_behp);
    n = 9, Wn[0] = 0.2; Wn[1] = 0; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs,
                                    Bandtype::HIGHPASS,
                                    IIRPrototype::BESSEL,
                                    ldigital)
    );
    EXPECT_EQ(ba, besselhp_ref);

    std::vector<double> bref_bebp({
          0.000784173933,   0.000000000000,  -0.007057565393,
          0.000000000000,   0.028230261570,  -0.000000000000,
         -0.065870610331,   0.000000000000,   0.098805915497,
         -0.000000000000,  -0.098805915497,  -0.000000000000,
          0.065870610331,  -0.000000000000,  -0.028230261570,
          0.000000000000,   0.007057565393,   0.000000000000,
         -0.000784173933});
    std::vector<double> aref_bebp({
          1.000000000000,  -4.048548350453,   8.990217360521,
        -14.294084391341,  18.444991832552, -20.230752051045,
         19.287967813235, -16.132032791176,  11.949696014765,
         -7.854853158848,   4.584510285738,  -2.362764193344,
          1.070076604012,  -0.420034501333,   0.140912159193,
         -0.039089313988,   0.008661875932,  -0.001375232998,
          0.000136336255});
    BA besselbp_ref(bref_bebp, aref_bebp);
    n = 9, Wn[0] = 0.2; Wn[1] = 0.6; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::BANDPASS,
                                    IIRPrototype::BESSEL,
                                    ldigital)
    );
    ba.setEqualityTolerance(1.e-10);
    EXPECT_EQ(ba, besselbp_ref);

    std::vector<double> bref_bebs({
        0.013644933989,  -0.093814218171,   0.409475147608,
       -1.261506650882,   3.083458581241,  -6.140072702909,
       10.321802128384, -14.782344380738,  18.305879752594,
      -19.621137325685,  18.305879752593, -14.782344380738,
       10.321802128384,  -6.140072702909,   3.083458581241,
       -1.261506650882,   0.409475147608,  -0.093814218171,
       0.013644933989});
    std::vector<double> aref_bebs({
       1.000000000000,  -4.089137154831,   9.232489406157,
      -14.875414921667,  19.249735800831, -20.865925048848,
       19.355134362319, -15.546782196294,  10.943020320099,
       -6.769919175677,   3.679974288098,  -1.750286695420,
        0.725970419042,  -0.259177483102,   0.078325672086,
       -0.019434146756,   0.003825453352,  -0.000536408493,
        0.000045365626});
    BA besselbs_ref(bref_bebs, aref_bebs);
    n = 9, Wn[0] = 0.2; Wn[1] = 0.6; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
       ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::BANDSTOP,
                                    IIRPrototype::BESSEL,
                                    ldigital)
    );
    ba.setEqualityTolerance(1.e-10); 
    EXPECT_EQ(ba, besselbs_ref);
}

TEST(UtilitiesDesignIIR, chebyshev1)
{
    int n;
    BA ba; 
    IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL;
    double rp = 5, rs = 60; 
    double Wn[2];
    //------------------------------------------------------------------------//
    //                           ChebyI filter designs                        //
    //------------------------------------------------------------------------//
    std::vector<double> bref_c1lp({
         0.000000000148,   0.000000001330,   0.000000005320,
         0.000000012414,   0.000000018621,   0.000000018621,
         0.000000012414,   0.000000005320,   0.000000001330,
         0.000000000148});
    std::vector<double> aref_c1lp({
         1.000000000000,  -8.652517885380,  33.487739440892,
       -76.084279249681, 111.825664229386, -110.255835235509,
        72.922994322042, -31.198411239773,   7.834505870746,
        -0.879860177057});
    BA cheby1lp_ref(bref_c1lp, aref_c1lp);
    n = 9, Wn[0] = 0.1; Wn[1] = 0; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::LOWPASS,
                                    IIRPrototype::CHEBYSHEV1,
                                    ldigital)
    );
    EXPECT_EQ(ba, cheby1lp_ref);

    std::vector<double> bref_c1hp({
         0.156526610694, -1.56526610694, 7.04369748122,
         -18.7831932833, 32.8705882457, -39.4447058949,
          32.8705882457, -18.7831932833, 7.04369748122,
         -1.56526610694, 0.156526610694});
    std::vector<double> aref_c1hp({
          1.0, -7.0698034649, 23.0809859182, -45.9927395117,
          62.5198783862, -61.5183650023, 45.3099617797,
         -25.1907008265, 10.241928725, -2.7419928524,
          0.362045627246});
    BA cheby1hp_ref(bref_c1hp, aref_c1hp);
    n = 10, Wn[0] = 0.1; Wn[1] = 0; rp = 5; rs = 60;
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs,
                                    Bandtype::HIGHPASS,
                                    IIRPrototype::CHEBYSHEV1,
                                    ldigital)
    );
    ba.setEqualityTolerance(1.e-10);
    EXPECT_EQ(ba, cheby1hp_ref);

    std::vector<double> bref_c1bp({
                        0.000042514612,   0.000000000000,  -0.000382631506,
                       -0.000000000000,   0.001530526024,   0.000000000000,
                       -0.003571227390,  -0.000000000000,   0.005356841085,
                        0.000000000000,  -0.005356841085,  -0.000000000000,
                        0.003571227390,  -0.000000000000,  -0.001530526024,
                       -0.000000000000,   0.000382631506,   0.000000000000,
                       -0.000042514612});
    std::vector<double> aref_c1bp({
                        1.000000000000,  -5.414904195421,  18.519819169954,
                      -46.188454029386,  93.937838026904, -160.515786539183,
                      238.218301075754, -309.983372667261, 358.532632489002,
                      -368.996247380907, 339.611173225014, -278.041979149337,
                      202.224208725540, -128.835879574125,  71.213949242754,
                      -33.010956002641,  12.465255168265,  -3.424310067682,
                        0.599734743788});
    BA cheby1bp_ref(bref_c1bp, aref_c1bp);
    n = 9, Wn[0] = 0.2; Wn[1] = 0.6; rp = 5; rs = 60;
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs,
                       Bandtype::BANDPASS, IIRPrototype::CHEBYSHEV1,
                       ldigital)
    );
    ba.setEqualityTolerance(1.e-10);
    EXPECT_EQ(ba, cheby1bp_ref);

    std::vector<double> bref_c1bs({
                        0.002146776265,  -0.014759920207,   0.064423289170,
                       -0.198474579554,   0.485124787140,  -0.966026098468,
                        1.623943350788,  -2.325726609324,   2.880089284293,
                       -3.087020570499,   2.880089284293,  -2.325726609324,
                        1.623943350788,  -0.966026098468,   0.485124787140,
                       -0.198474579554,   0.064423289170,  -0.014759920207,
                        0.002146776265});
    std::vector<double> aref_c1bs({
                        1.000000000000,  -2.980992943560,   3.101499838780,
                       -2.451457479529,   5.029662024748,  -5.795331906920,
                        1.377753306580,  -0.991118154963,   4.253677611093,
                       -1.078796688361,  -2.559091332204,  -0.423920836411,
                        1.177610357720,   2.186576639770,  -2.039232376476,
                        0.188604144355,  -0.765448786084,   1.249442240016,
                       -0.464975668856});
    BA cheby1bs_ref(bref_c1bs, aref_c1bs);
    n = 9, Wn[0] = 0.2; Wn[1] = 0.6; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::BANDSTOP,
                                    IIRPrototype::CHEBYSHEV1,
                                    ldigital)
    );
    ba.setEqualityTolerance(1.e-10);
    EXPECT_EQ(ba, cheby1bs_ref);
}

TEST(UtilitiesDesignIIR, chebyshev2)
{
    int n;
    BA ba; 
    IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL;
    double rp = 5, rs = 60; 
    double Wn[2];
    //------------------------------------------------------------------------//
    //                           ChebyII filter designs                       //
    //------------------------------------------------------------------------//
    std::vector<double> bref_c2lp({
          0.00410229901970, -0.00377904884003, 0.01083454827137, 
         -0.00281635851952,  0.01021847244343, 0.00076410015492, 
          0.01021847244343, -0.00281635851952, 0.01083454827137, 
         -0.00377904884003,  0.00410229901970});
    std::vector<double> aref_c2lp({
          1.00000000000000, -4.61066319266803, 10.47263927023236, 
        -14.88830441791697, 14.51826097876597, -10.05409602415944, 
          4.98526570054316, -1.73895205154971, 0.40768256556736, 
         -0.05771936784473,  0.00377046393483});
    BA cheby2lp_ref(bref_c2lp, aref_c2lp);
    n = 10, Wn[0] = 0.3; Wn[1] = 0; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::LOWPASS,
                                    IIRPrototype::CHEBYSHEV2,
                                    ldigital)
    );
    EXPECT_EQ(ba, cheby2lp_ref);

    std::vector<double> bref_c2hp({
         0.43293356378981,  -4.22273302625180,   18.63839282119448, 
       -49.02204942735300,  85.08305568194051, -101.81919921855661, 
        85.08305568194051, -49.02204942735300,   18.63839282119448, 
        -4.22273302625180,  0.43293356378981});
    std::vector<double> aref_c2hp({
         1.00000000000000,  -8.11342805009665,   29.93816281599581, 
       -66.11854688628588,  96.73306698067843,  -97.91572875594167, 
        69.41876516736002, -34.02544243130928,   11.03134173680539, 
        -2.13561396448704,  0.18743147065574});
    BA cheby2hp_ref(bref_c2hp, aref_c2hp);
    n = 10, Wn[0] = 0.1; Wn[1] = 0; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::HIGHPASS,
                                    IIRPrototype::CHEBYSHEV2,
                                    ldigital)
    );
    ba.setEqualityTolerance(1.e-10);
    EXPECT_EQ(ba, cheby2hp_ref);

    std::vector<double> bref_c2bp({
        0.01780621293436, -0.06696882576862,  0.06108344898157, 
        0.01004581026441,  0.07349345972810, -0.17060184129854, 
       -0.05121355317387,  0.10282917034156,  0.20601137119488, 
       -0.06846583336487, -0.22803865239622, -0.06846583336489, 
        0.20601137119492,  0.10282917034153, -0.05121355317385, 
       -0.17060184129855,  0.07349345972810,  0.01004581026441, 
        0.06108344898157, -0.06696882576862,  0.01780621293436});
    std::vector<double> aref_c2bp({
        1.00000000000000,  -6.84255089779208,  22.40616706864133, 
      -48.24813970119118,  79.20429150293148, -106.74490358686622, 
      121.38665124501451, -117.95186383396000, 99.45754979982192, 
      -73.49753937532221,  47.45468657563394, -26.70866962293626, 
       13.18874996326441, -5.66506565650582,   2.05161013482263, 
       -0.62948790279196,  0.16999211512400,  -0.03407552184754, 
       0.00318343314108,  -0.00074355290758,  0.00034509546998});
    BA cheby2bp_ref(bref_c2bp, aref_c2bp);
    n = 10, Wn[0] = 0.1; Wn[1] = 0.6; rp = 5; rs = 60; 
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::BANDPASS,
                                    IIRPrototype::CHEBYSHEV2,
                                    ldigital)
    );
    EXPECT_EQ(ba, cheby2bp_ref);

    std::vector<double> bref_c2bs({
             0.03403410673803, -0.16827387273319,  0.54220130005646, 
            -1.26791509475874,  2.42988115437137, -3.89063719715933, 
             5.40220598811858, -6.52370289233002,  6.96037794524616, 
            -6.52370289233001,  5.40220598811858, -3.89063719715933, 
             2.42988115437137, -1.26791509475874,  0.54220130005646, 
            -0.16827387273319,  0.03403410673803});
    std::vector<double> aref_c2bs({
            1.00000000000000,  -3.13888032989080,  4.71448903765691, 
           -5.40209994509873,   6.50241192219140, -7.14166953696394, 
            6.24508961724931,  -4.78833659631079,  3.60582442601859, 
           -2.42823290403137,   1.36134353524197, -0.67602590382599, 
            0.31219005162189,  -0.11821905447571,  0.03444463360583, 
           -0.00759384336523,   0.00122982022910});
    BA cheby2bs_ref(bref_c2bs, aref_c2bs);
    n = 8, Wn[0] = 0.2; Wn[1] = 0.6; rp = 5; rs = 60;
    EXPECT_NO_THROW(
        ba = IIR::designBAIIRFilter(n, Wn, rp, rs, 
                                    Bandtype::BANDSTOP,
                                    IIRPrototype::CHEBYSHEV2,
                                    ldigital);
    );
    ba.setEqualityTolerance(1.e-10);
    EXPECT_EQ(ba, cheby2bs_ref);
}

/// IIR SOS analog prototype design.
//int rtseis_test_utils_design_zpk2sos(void)
TEST(UtilitiesDesignIIR, zpk2sos)
{
    SOS sos;
    int n = 4;
    double Wn[1] = {0.1};
    const std::vector<double> bsRef1({
         4.16599204e-04,   8.33198409e-04,   4.16599204e-04,
         1.00000000e+00,   2.00000000e+00,   1.00000000e+00});
    const std::vector<double> asRef1({
         1.,        -1.47967422,  0.55582154,
         1.,        -1.70096433,  0.78849974});
    SOS sosRef1(2, bsRef1, asRef1);
    const std::vector<double> bsRefEll({
         0.0014154,   0.00248707,  0.0014154,
         1.,          0.72965193,  1., 
         1.,          0.17594966,  1.});
    const std::vector<double> asRefEll({
         1.,        -1.32543251,  0.46989499,
         1.,        -1.26117915,  0.6262586,
         1.,        -1.25707217,  0.86199667});
    SOS sosRefEll(3, bsRefEll, asRefEll);
    IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL;
    EXPECT_NO_THROW(
         sos = IIR::designSOSIIRFilter(n, Wn, 5, 60,  
                                       Bandtype::LOWPASS,
                                       IIRPrototype::BUTTERWORTH,
                                       ldigital, SOSPairing::NEAREST)
    );  
    sos.setEqualityTolerance(1.e-5);
    EXPECT_EQ(sos, sosRef1);

    //  
    std::vector<std::complex<double>> zell;
    zell.push_back(std::complex<double> (-0.878578886634,  0.47759725707));
    zell.push_back(std::complex<double> (-0.364825965978,  0.93107572976));
    zell.push_back(std::complex<double> (-0.0879748281791, 0.996122698068));
    zell.push_back(std::complex<double> (-0.878578886634, -0.47759725707));
    zell.push_back(std::complex<double> (-0.364825965978, -0.93107572976));
    zell.push_back(std::complex<double> (-0.0879748281791,-0.996122698068));
    std::vector<std::complex<double>> pell;
    pell.push_back(std::complex<double> (0.662716257451,-0.175220303539)); 
    pell.push_back(std::complex<double> (0.630589576722,-0.478137415927));
    pell.push_back(std::complex<double> (0.62853608609, -0.683329390963));
    pell.push_back(std::complex<double> (0.662716257451,+0.175220303539));
    pell.push_back(std::complex<double> (0.630589576722,+0.478137415927));
    pell.push_back(std::complex<double> (0.62853608609, +0.683329390963));
    double kell = 0.00141539634442;
    ZPK zpk(zell, pell, kell);
    EXPECT_NO_THROW(sos = IIR::zpk2sos(zpk, SOSPairing::NEAREST));
    sos.setEqualityTolerance(1.e-5);
    EXPECT_EQ(sos, sosRefEll);

    // Test this edge case so i can remove a debugging statement
    // iirfilter(2, 10/(100), rp=60, btype='lowpass', ftype='cheby1', output='sos')
    std::vector<double> bsRefCheb1 = {1.2386078193258956e-05,
                                      2.477215638651791e-05,
                                      1.2386078193258956e-05};  
    std::vector<double> asRefCheb1 = {1.0, -1.9502344968431102, 0.999778809616146};
    SOS sosRefCheb1(1, bsRefCheb1, asRefCheb1);
    double Wn2[2] = {10/100., 0}; 
    EXPECT_NO_THROW(
         sos = IIR::designSOSIIRFilter(2, Wn2, 60, 0,
                                       Bandtype::LOWPASS,
                                       IIRPrototype::CHEBYSHEV1,
                                       ldigital, SOSPairing::NEAREST)
    );  
    sos.setEqualityTolerance(1.e-12);
    EXPECT_EQ(sos, sosRefCheb1);
}

}
