#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ipps.h>
#include "rtseis/utilities/filterDesign/fir.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis;
using namespace RTSeis::Utilities::FilterDesign;

/*! Test FIR-window based filter design (fir1 in matlab). */
//int rtseis_test_utils_design_fir_fir1(void)
TEST(UtilitiesDesignFIR, lowpass)
{
    // b = fir1(13, 0.6); % default hamming
    std::vector<double> ref1({-0.001213632932920, -0.006228195759795,
                               0.015988039487026,  0.013651616475176,
                              -0.089746334361735,  0.058133494212324,
                               0.509415012879924,  0.509415012879924,
                               0.058133494212324, -0.089746334361735,
                               0.013651616475176,  0.015988039487026,
                              -0.006228195759795, -0.001213632932920});
    FilterRepresentations::FIR fir;
    // Test 1
    int order = 13;
    double r = 0.6;
    EXPECT_NO_THROW(fir = FIR::FIR1Lowpass(order, r, FIRWindow::HAMMING));
    auto ptaps = fir.getFilterTaps();
    EXPECT_EQ(ptaps.size(), ref1.size());
    double error;
    ippsNormDiff_Inf_64f(ptaps.data(), ref1.data(), ptaps.size(), &error);
    ASSERT_LE(error, 1.e-12);
}

TEST(UtilitiesDesignFIR, highpass)
{
    // fir1(16, 0.45, 'high', hann(16+1));
    std::vector<double> ref2({ 0,                  0.000785256397236,
                              -0.006281697554181, -0.013886226413618,
                               0.023373298499376,  0.065319627573863,
                              -0.041954095979927, -0.302244991506668,
                               0.549672322171091, -0.302244991506668,
                              -0.041954095979927,  0.065319627573863,
                               0.023373298499376, -0.013886226413618,
                              -0.006281697554181,  0.000785256397236,
                               0});
    FilterRepresentations::FIR fir;
    int order = 16;
    double r = 0.45;
    EXPECT_NO_THROW(fir = FIR::FIR1Highpass(order, r, FIRWindow::HANN));
    auto ptaps = fir.getFilterTaps();
    EXPECT_EQ(ptaps.size(), ref2.size());
    double error;
    ippsNormDiff_Inf_64f(ptaps.data(), ref2.data(), ptaps.size(), &error);
    ASSERT_LE(error, 1.e-12);
}

TEST(UtilitiesDesignFIR, bandpass)
{
    // fir1(11, [0.2, 0.8], 'bandpass', bartlett(11+1));
    std::vector<double> ref3({ 0,                 -0.018100893326756,
                              -0.008171956220958, -0.077570731498438,
                              -0.240779937850341,  0.415028961385482,
                               0.415028961385482, -0.240779937850341,
                              -0.077570731498438, -0.008171956220958,
                              -0.018100893326756,  0});
    std::pair<double,double> r2;
    int order = 11;
    r2 = std::make_pair(0.2, 0.8);
    FilterRepresentations::FIR fir;
    EXPECT_NO_THROW(fir = FIR::FIR1Bandpass(order, r2, FIRWindow::BARTLETT));
    auto ptaps = fir.getFilterTaps();
    EXPECT_EQ(ptaps.size(), ref3.size());
    double error;
    ippsNormDiff_Inf_64f(ptaps.data(), ref3.data(), ptaps.size(), &error);
    ASSERT_LE(error, 1.e-12);
}

TEST(UtilitiesDesignFIR, bandstop)
{
    // IPP uses optimal blackman
    // alpha = -0.5/(1 + cos(2*pi/20));
    // fir1(20, [0.15,0.85], 'stop', (alpha+1)/2 - 0.5*cos(2*pi*linspace(0,20,21)/20) - alpha/2*cos(4*pi*linspace(0,20,21)/20))
    std::vector<double> ref4({ 0.000000000000000,  0,
                              -0.000380306566042,  0.000000000000000,
                               0.004359748553493,  0.000000000000000,
                               0.074832232270797,  0.000000000000000,
                               0.245754965804324,                  0,
                               0.350866719874856,                  0,
                               0.245754965804324,  0.000000000000000,
                               0.074832232270797,  0.000000000000000,
                               0.004359748553493,  0.000000000000000,
                              -0.000380306566042, -0.000000000000000,
                               0});
    std::pair<double,double> r2;
    int order = 20;
    r2 = std::make_pair(0.15, 0.85);
    FilterRepresentations::FIR fir;
    EXPECT_NO_THROW(fir = FIR::FIR1Bandstop(order, r2, FIRWindow::BLACKMAN_OPT));
    auto ptaps = fir.getFilterTaps();
    EXPECT_EQ(ptaps.size(), ref4.size());
    double error;
    ippsNormDiff_Inf_64f(ptaps.data(), ref4.data(), ptaps.size(), &error);
    ASSERT_LE(error, 1.e-12);
}

TEST(UtilitiesDesignFIR, hilbert)
{
    // Edge case
    int order = 0;
    std::pair<FilterRepresentations::FIR, FilterRepresentations::FIR> hilbertFIR;
    EXPECT_NO_THROW(hilbertFIR = FIR::HilbertTransformer(order, 8));
    auto realTaps = hilbertFIR.first.getFilterTaps();
    auto imagTaps = hilbertFIR.second.getFilterTaps();
    EXPECT_EQ(static_cast<int> (realTaps.size()), 1);
    EXPECT_EQ(static_cast<int> (imagTaps.size()), 1);
    EXPECT_NEAR(realTaps[0], 1, 1.e-14);
    EXPECT_NEAR(imagTaps[0], 0, 1.e-14);

    // NOTE that there are two even cases to test.  The first or second elemetn
    // of the imaginary part can be non-zero for the TypeIII filter.
#ifdef __STDCPP_MATH_SPEC_FUNCS__ // __cplusplus > 201402L
    double tol = 1.e-14;
#else
    double tol = 1.e-6;
#endif
    const std::vector<double> hfiltR16 = {-0.000000000000000,   0.000000000000000,
       -0.000000000000000,   0.000000000000000,  -0.000000000000000,
        0.000000000000000,  -0.000000000000000,   0.000000000000000,
        1.000000000000000,   0.000000000000000,  -0.000000000000000,
        0.000000000000000,  -0.000000000000000,   0.000000000000000,
       -0.000000000000000,   0.000000000000000,  -0.000000000000000};
    const std::vector<double> hfiltI16 = {-0.000000000000000,  -0.002154949073886,
       -0.000000000000000,  -0.025048021740503,  -0.000000000000000,
       -0.123110554263836,  -0.000000000000000,  -0.600346453738133,
                        0,   0.600346453738133,   0.000000000000000,
        0.123110554263836,   0.000000000000000,   0.025048021740503,
        0.000000000000000,   0.002154949073886,   0.000000000000000};
    order = 16; 
    EXPECT_NO_THROW(hilbertFIR = FIR::HilbertTransformer(order, 8));
    realTaps = hilbertFIR.first.getFilterTaps();
    imagTaps = hilbertFIR.second.getFilterTaps();
    EXPECT_EQ(static_cast<int> (realTaps.size()), order + 1);
    EXPECT_EQ(static_cast<int> (imagTaps.size()), order + 1);
    double error;
    ippsNormDiff_Inf_64f(realTaps.data(), hfiltR16.data(),
                         realTaps.size(), &error);
    ASSERT_LE(error, tol);
    ippsNormDiff_Inf_64f(imagTaps.data(), hfiltI16.data(), 
                         imagTaps.size(), &error);
    ASSERT_LE(error, tol);

    const std::vector<double> hfiltR19 = {-0.000078362668030, 0.000687065652501,
        -0.002495774178985,   0.006621258126556,  -0.014693594258722,
         0.029092392815675,  -0.053803945082584,   0.097839237445465,
        -0.193197503788858,   0.630029225936982,   0.630029225936982,
        -0.193197503788858,   0.097839237445465,  -0.053803945082584,
         0.029092392815675,  -0.014693594258722,   0.006621258126556,
        -0.002495774178985,   0.000687065652501,  -0.000078362668030};
    const std::vector<double> hfiltI19 = {-0.000078362668030,  -0.000687065652501,
        -0.002495774178985,  -0.006621258126556,  -0.014693594258722,
        -0.029092392815675,  -0.053803945082584,  -0.097839237445464,
        -0.193197503788858,  -0.630029225936981,   0.630029225936981,
         0.193197503788858,   0.097839237445464,   0.053803945082584,
         0.029092392815675,   0.014693594258722,   0.006621258126556,
         0.002495774178985,   0.000687065652501,   0.000078362668030};
    order = 19; 
    EXPECT_NO_THROW(hilbertFIR = FIR::HilbertTransformer(order, 8));
    realTaps = hilbertFIR.first.getFilterTaps();
    imagTaps = hilbertFIR.second.getFilterTaps();
    EXPECT_EQ(static_cast<int> (realTaps.size()), order + 1);
    EXPECT_EQ(static_cast<int> (imagTaps.size()), order + 1);
    ippsNormDiff_Inf_64f(realTaps.data(), hfiltR19.data(),
                         realTaps.size(), &error);
    ASSERT_LE(error, tol);
    ippsNormDiff_Inf_64f(imagTaps.data(), hfiltI19.data(),
                         imagTaps.size(), &error);
    ASSERT_LE(error, tol);
}

}
