#include <cstdio>
#include <cstdlib>
#include <string>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/filterDesign/response.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis;
using namespace RTSeis::FilterRepresentations;
using namespace RTSeis::FilterDesign;

//int rtseis_test_utils_design_freqs(void)
TEST(UtilitiesResponse, Freqs)
{
    const int nw = 50;
    std::vector<double> bs({1611.7315706,  0.,  0.,  0.,  0.});
    std::vector<double> as({1.00000000e+00,   8.57530241e+00,   1.57767906e+02,
                            7.98628595e+02,   4.76375068e+03,   7.98628595e+03,
                            1.57767906e+04,   8.57530241e+03,   1.00000000e+04});
    BA ba(bs, as);
    std::vector<double> bz({0.056340000000000, -0.000935244000000,
                           -0.000935244000000,  0.056340000000000});
    std::vector<double> az({1.000000000000000, -2.129100000000000,
                            1.783386300000000, -0.543463100000000});
    BA baz(bz, az);
    // make solutions for freqs
    std::vector<std::complex<double>> href1; href1.resize(nw);
    href1[0] = std::complex<double> (+0.03488577413550,-0.02699025100928);
    href1[1] = std::complex<double> (+0.05348973603688,-0.05152025023515);
    href1[2] = std::complex<double> (+0.07995180799247,-0.10404220047713);
    href1[3] = std::complex<double> (+0.10278736966951,-0.22594753568789);
    href1[4] = std::complex<double> (+0.00848918600791,-0.50531182099015);
    href1[5] = std::complex<double> (-0.62270477129993,-0.66247486038261);
    href1[6] = std::complex<double> (-0.98269958366215,-0.00383789980467);
    href1[7] = std::complex<double> (-0.80405788502503,+0.42332819774255);
    href1[8] = std::complex<double> (-0.60615492623016,+0.65515371685383);
    href1[9] = std::complex<double> (-0.39658757122898,+0.83087903105862);
    href1[10] = std::complex<double> (-0.13970324532235,+0.95304016024704);
    href1[11] = std::complex<double> (+0.15497598941797,+0.98179877285411);
    href1[12] = std::complex<double> (+0.43290796318079,+0.89967494452506);
    href1[13] = std::complex<double> (+0.64293416921331,+0.73874321858729);
    href1[14] = std::complex<double> (+0.77441733458820,+0.54927348734769);
    href1[15] = std::complex<double> (+0.84574252312641,+0.36381879901862);
    href1[16] = std::complex<double> (+0.87952203818244,+0.19215869533035);
    href1[17] = std::complex<double> (+0.89095834680252,+0.03087283491702);
    href1[18] = std::complex<double> (+0.88609894228890,-0.12857036725930);
    href1[19] = std::complex<double> (+0.86230825135038,-0.29537391207665);
    href1[20] = std::complex<double> (+0.80798359282187,-0.47566481744457);
    href1[21] = std::complex<double> (+0.70276882681403,-0.66602904861868);
    href1[22] = std::complex<double> (+0.52426157690959,-0.84353059262555);
    href1[23] = std::complex<double> (+0.26865301465227,-0.96221985462147);
    href1[24] = std::complex<double> (-0.02714612573254,-0.97737881098598);
    href1[25] = std::complex<double> (-0.30297302526253,-0.88677677485093);
    href1[26] = std::complex<double> (-0.52857685878989,-0.72830931064414);
    href1[27] = std::complex<double> (-0.72551867160719,-0.52450849935703);
    href1[28] = std::complex<double> (-0.92781800750244,-0.20388405324537);
    href1[29] = std::complex<double> (-0.88469336203972,+0.44215940372882);
    href1[30] = std::complex<double> (-0.15481996440243,+0.64414889509246);
    href1[31] = std::complex<double> (+0.09562472100872,+0.31084903769236);
    href1[32] = std::complex<double> (+0.09120258029299,+0.13942195037108);
    href1[33] = std::complex<double> (+0.06291517417591,+0.06720800519134);
    href1[34] = std::complex<double> (+0.04122820755716,+0.03451868927268);
    href1[35] = std::complex<double> (+0.02685204970343,+0.01856698156418);
    href1[36] = std::complex<double> (+0.01756607670780,+0.01032696165320);
    href1[37] = std::complex<double> (+0.01156809207517,+0.00588806178611);
    href1[38] = std::complex<double> (+0.00766837172005,+0.00342072639359);
    href1[39] = std::complex<double> (+0.00511306514726,+0.00201622740725);
    href1[40] = std::complex<double> (+0.00342630435821,+0.00120189184245);
    href1[41] = std::complex<double> (+0.00230562492063,+0.00072288480789);
    href1[42] = std::complex<double> (+0.00155691392926,+0.00043789173559);
    href1[43] = std::complex<double> (+0.00105437478338,+0.00026677845884);
    href1[44] = std::complex<double> (+0.00071575337674,+0.00016328428015);
    href1[45] = std::complex<double> (+0.00048684418034,+0.00010031589633);
    href1[46] = std::complex<double> (+0.00033168506020,+0.00006181943458);
    href1[47] = std::complex<double> (+0.00022628105324,+0.00003819156150);
    href1[48] = std::complex<double> (+0.00015454510105,+0.00002364289385);
    href1[49] = std::complex<double> (+0.00010564846092,+0.00001466105711);

    const double x1 =-1;
    const double x2 = 1;
    const double dx = (x2 - x1)/static_cast<double> (nw - 1);
    std::vector<double> w(nw); //w.resize(nw);
    for (int i=0; i<nw; i++)
    {
        w[i] = 2.0*M_PI*std::pow(10, x1 + static_cast<double> (i)*dx);
    }
    std::vector<std::complex<double>> h;
    EXPECT_NO_THROW(h = Response::freqs(ba, w));
    EXPECT_EQ(h.size(), static_cast<size_t> (nw));
    double emax = 0;
    for (auto i=0; i<static_cast<int> (h.size()); i++)
    {
        auto res = std::abs(h[i] - href1[i]);
        emax = std::max(res, emax);
        if (res > 1.e-10) //std::abs(h[i] - href1[i]) > 1.e-10)
        {
            fprintf(stderr, "Failed freqs test (%lf,%lf) (%lf,%lf)",
                    std::real(h[i]), std::imag(h[i]),
                    std::real(href1[i]), std::imag(href1[i]));
        }
    }
    EXPECT_LE(emax, 1.e-10);
}

TEST(UtilitiesResponse, Freqz)
{
    const int nf = 41; 
    std::vector<double> bz({0.056340000000000, -0.000935244000000,
                           -0.000935244000000,  0.056340000000000});
    std::vector<double> az({1.000000000000000, -2.129100000000000,
                            1.783386300000000, -0.543463100000000});
    BA baz(bz, az);
    // make solution for freqz
    std::vector<std::complex<double>> href2(nf);
    href2[0] = std::complex<double> (+0.99987648795559,+0.00000000000000);
    href2[1] = std::complex<double> (+0.95601820836871,-0.24635538126662);
    href2[2] = std::complex<double> (+0.84182473416863,-0.45234494448287);
    href2[3] = std::complex<double> (+0.69218147536286,-0.60548822657119);
    href2[4] = std::complex<double> (+0.53073903024459,-0.72046222089002);
    href2[5] = std::complex<double> (+0.35389724345662,-0.82273951176198);
    href2[6] = std::complex<double> (+0.11450731879359,-0.92765406032472);
    href2[7] = std::complex<double> (-0.30356986746249,-0.94904121363739);
    href2[8] = std::complex<double> (-0.74535194596939,-0.48906787904061);
    href2[9] = std::complex<double> (-0.52345405359234,+0.03883239252739);
    href2[10] = std::complex<double> (-0.23022242683644,+0.11340130414619);
    href2[11] = std::complex<double> (-0.09247116985374,+0.07419998418252);
    href2[12] = std::complex<double> (-0.03163529655289,+0.03366165116973);
    href2[13] = std::complex<double> (-0.00365030620855,+0.00475347063792);
    href2[14] = std::complex<double> (+0.00945480661491,-0.01445905891768);
    href2[15] = std::complex<double> (+0.01535492883430,-0.02690992479787);
    href2[16] = std::complex<double> (+0.01759436581746,-0.03478100394682);
    href2[17] = std::complex<double> (+0.01792468519841,-0.03954002587212);
    href2[18] = std::complex<double> (+0.01725815161918,-0.04216037245807);
    href2[19] = std::complex<double> (+0.01608849455770,-0.04329285168849);
    href2[20] = std::complex<double> (+0.01468989855489,-0.04337885254581);
    href2[21] = std::complex<double> (+0.01321670410931,-0.04272230968619);
    href2[22] = std::complex<double> (+0.01175579491817,-0.04153550307532);
    href2[23] = std::complex<double> (+0.01035521965306,-0.03996857117507);
    href2[24] = std::complex<double> (+0.00904033577320,-0.03812884046514);
    href2[25] = std::complex<double> (+0.00782315309451,-0.03609370155880);
    href2[26] = std::complex<double> (+0.00670785329979,-0.03391932843019);
    href2[27] = std::complex<double> (+0.00569410163719,-0.03164667442229);
    href2[28] = std::complex<double> (+0.00477905564053,-0.02930565474428);
    href2[29] = std::complex<double> (+0.00395859088629,-0.02691810244844);
    href2[30] = std::complex<double> (+0.00322804947486,-0.02449988289502);
    href2[31] = std::complex<double> (+0.00258269445972,-0.02206242321480);
    href2[32] = std::complex<double> (+0.00201798188745,-0.01961383022046);
    href2[33] = std::complex<double> (+0.00152971946493,-0.01715971571506);
    href2[34] = std::complex<double> (+0.00111415501448,-0.01470381187907);
    href2[35] = std::complex<double> (+0.00076802196457,-0.01224843497540);
    href2[36] = std::complex<double> (+0.00048855920133,-0.00979483895752);
    href2[37] = std::complex<double> (+0.00027351634834,-0.00734348911731);
    href2[38] = std::complex<double> (+0.00012115155313,-0.00489427799722);
    href2[39] = std::complex<double> (+0.00003022628902,-0.00244670032556);
    href2[40] = std::complex<double> (-0.00000000000000,-0.00000000000000);
    std::vector<double> w(nf);
    double df = (M_PI - 0)/static_cast<double> (nf - 1);
    for (auto i=0; i<nf; i++)
    {
        w[i] = 0 + static_cast<double> (i)*df;
    }  
    std::vector<std::complex<double>> h;
    EXPECT_NO_THROW(h = Response::freqz(baz, w));
    EXPECT_EQ(h.size(), static_cast<size_t> (nf));
    double emax = 0;
    for (auto i=0; i<static_cast<int> (h.size()); i++)
    {
        auto res = std::abs(h[i] - href2[i]);
        emax = std::max(res, emax);
        if (std::abs(h[i] - href2[i]) > 1.e-10)
        {
            fprintf(stderr, "Failed freqz test (%lf,%lf) (%lf,%lf)",
                    std::real(h[i]), std::imag(h[i]),
                    std::real(href2[i]), std::imag(href2[i]));
        }
    }
    EXPECT_LE(emax, 1.e-10);
}

//int rtseis_test_utils_design_groupDelay(void)
TEST(UtilitiesResponse, GroupDelay)
{
    const int nf = 20;
    const std::vector<double> b({4.88711377891e-05, 0.000195484551156,
                                 0.000293226826734, 0.000195484551156,
                                 4.88711377891e-05});
    const std::vector<double> a({1.0, -3.77274406567, 5.43219217707,
                                 -3.53522638245, 0.877168775657});
    const std::vector<double> gdR({5.85741709682, 5.86011546331, 5.86821742633,
                                   5.88174360001, 5.90072841926, 5.92522025663,
                                   5.95528158551, 5.99098918908, 6.03243441439,
                                   6.07972347054, 6.13297776959, 6.19233430846,
                                   6.25794608936, 6.32998257563, 6.40863017903,
                                   6.49409277321, 6.58659222689, 6.68636894831,
                                   6.79368243081, 6.90881178634});
    BA ba(b, a);
    std::vector<double> w(nf);
    double di = ((2*M_PI)/180)/static_cast<double> (nf - 1);
    for (int i=0; i<nf; i++)
    {
        w[i] =  0 + di*static_cast<double> (i);
    }
    std::vector<double> gd;
    EXPECT_NO_THROW(gd = Response::groupDelay(ba, w));
    double emax;
    EXPECT_EQ(gd.size(), static_cast<size_t> (nf));
    ippsNormDiff_Inf_64f(gd.data(), gdR.data(), nf, &emax);
    EXPECT_LE(emax, 1.e-7);

/*
    // group_delay((b, a), 20, whole=True)
    ierr = Response::groupDelay(ba, gd, 20, true);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute group delay");
        return EXIT_FAILURE;
    }
    for (int i=0; i<20; i++)
    {
        printf("%lf\n", gd[i]);
    }
getchar();
*/
}

}
