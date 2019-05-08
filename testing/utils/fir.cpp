#include <cstdio>
#include <cstdlib>
#include <cassert>
#define RTSEIS_LOGGING 1
#include <cmath>
#include "rtseis/utilities/design/fir.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/log.h"
#include "utils.hpp"

using namespace RTSeis::Utilities;
using namespace RTSeis::Utilities::FilterDesign;

/*! Test FIR-window based filter design (fir1 in matlab). */
int rtseis_test_utils_design_fir_fir1(void)
{
    int order = 13;
    double r = 0.6;
    // b = fir1(13, 0.6); % default hamming
    std::vector<double> ref1({-0.001213632932920, -0.006228195759795,
                               0.015988039487026,  0.013651616475176,
                              -0.089746334361735,  0.058133494212324,
                               0.509415012879924,  0.509415012879924,
                               0.058133494212324, -0.089746334361735,
                               0.013651616475176,  0.015988039487026,
                              -0.006228195759795, -0.001213632932920});
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
    // fir1(11, [0.2, 0.8], 'bandpass', bartlett(11+1));
    std::vector<double> ref3({ 0,                 -0.018100893326756,
                              -0.008171956220958, -0.077570731498438,
                              -0.240779937850341,  0.415028961385482,
                               0.415028961385482, -0.240779937850341,
                              -0.077570731498438, -0.008171956220958,
                              -0.018100893326756,  0});
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
    FilterRepresentations::FIR fir;
    std::vector<double> ptaps;
    // Test 1
    order = 13; r = 0.6; 
    try
    {
        fir = FIR::FIR1Lowpass(order, r, FIRWindow::HAMMING);
    }
    catch (std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to generate lowpass filter %s", ia.what());
        return EXIT_FAILURE;
    }
    ptaps = fir.getFilterTaps();
    if (ptaps.size() != ref1.size())
    {
        RTSEIS_ERRMSG("%s", "Inconsistent sizes");
        return EXIT_FAILURE;
    }
    for (size_t i=0; i<ref1.size(); i++)
    {
        if (std::abs(ptaps[i] - ref1[i]) > 1.e-12)
        {
            RTSEIS_ERRMSG("Failed design %lf %lf", ref1[i], ptaps[i]);
            return EXIT_FAILURE;
        }
    }
    // Test 2
    order = 16; r = 0.45;
    try
    {
        fir = FIR::FIR1Highpass(order, r, FIRWindow::HANN);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to generate highpass filter %s", ia.what());
        return EXIT_FAILURE;
    }
    ptaps = fir.getFilterTaps();
    if (ptaps.size() != ref2.size())
    {
        RTSEIS_ERRMSG("%s", "Inconsistent sizes");
        return EXIT_FAILURE;
    }
    for (size_t i=0; i<ref2.size(); i++)
    {
        if (std::abs(ptaps[i] - ref2[i]) > 1.e-12)
        {
            RTSEIS_ERRMSG("Failed design %lf %lf", ref2[i], ptaps[i]);
            return EXIT_FAILURE;
        }
    }
    // Test 3
    std::pair<double,double> r2;
    order = 11;
    r2 = std::make_pair(0.2, 0.8);
    try
    {
        fir = FIR::FIR1Bandpass(order, r2, FIRWindow::BARTLETT); 
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to generate bandpass filter %s", ia.what());
        return EXIT_FAILURE;
    }
    ptaps = fir.getFilterTaps();
    if (ptaps.size() != ref3.size())
    {
        RTSEIS_ERRMSG("%s", "Inconsistent sizes");
        return EXIT_FAILURE;
    }
    for (size_t i=0; i<ref3.size(); i++)
    {
        if (std::abs(ptaps[i] - ref3[i]) > 1.e-12)
        {
            RTSEIS_ERRMSG("Failed design %lf %lf", ref3[i], ptaps[i]);
            return EXIT_FAILURE;
        }
    }
    // Test 4
    order = 20;
    r2 = std::make_pair(0.15, 0.85);
    try
    {
        fir = FIR::FIR1Bandstop(order, r2, FIRWindow::BLACKMAN_OPT);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed to generate bandstop filter %s", ia.what());
        return EXIT_FAILURE;
    }
    ptaps = fir.getFilterTaps();
    if (ptaps.size() != ref4.size())
    {
        RTSEIS_ERRMSG("%s", "Inconsistent sizes");
        return EXIT_FAILURE;
    }
    for (size_t i=0; i<ref4.size(); i++)
    {
        if (std::abs(ptaps[i] - ref4[i]) > 1.e-12)
        {
            RTSEIS_ERRMSG("Failed design %lf %lf", ref4[i], ptaps[i]);
            return EXIT_FAILURE;
        }
    }
    // Test 5 - Hilbert
    order = 0;
    std::vector<double> realTaps, imagTaps;
    std::pair<FilterRepresentations::FIR, FilterRepresentations::FIR> hilbertFIR;
    try
    {
        hilbertFIR = FIR::HilbertTransformer(order, 8);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed hilbert length 1 design %s", ia.what());
        return EXIT_FAILURE;
    }
    realTaps = hilbertFIR.first.getFilterTaps();
    imagTaps = hilbertFIR.second.getFilterTaps();
    if (realTaps.size() != 1 || imagTaps.size() != 1)
    {
        RTSEIS_ERRMSG("%s", "Hilbert is wrong length");
        return EXIT_FAILURE;
    }
    if (std::abs(realTaps[0] - 1) > 1.e-14 || std::abs(imagTaps[0]) > 1.e-14)
    {
        RTSEIS_ERRMSG("%s", "Hilbert order 0 is wrong");
        return EXIT_FAILURE;
    }
   
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
    try 
    {
        hilbertFIR = FIR::HilbertTransformer(order, 8); 
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed hilbert length 9 design %s", ia.what()); 
        return EXIT_FAILURE;
    }
    realTaps = hilbertFIR.first.getFilterTaps();
    imagTaps = hilbertFIR.second.getFilterTaps();
    if (realTaps.size() != order + 1 || imagTaps.size() != order + 1)
    {
        RTSEIS_ERRMSG("%s", "Hilbert is wrong length");
        return EXIT_FAILURE;
    }
    for (int i=0; i<order+1; i++)
    {   
        if (std::abs(hfiltR16[i] - realTaps[i]) > tol ||
            std::abs(hfiltI16[i] - imagTaps[i]) > tol)
        {   
            RTSEIS_ERRMSG("hilbert failed %d (%e,%e), (%e,%e)\n", i,
                           hfiltR16[i], realTaps[i], hfiltI16[i], imagTaps[i]);

        }   
    }
 
    const std::vector<double> hfiltR18 = {0.000000000000000, -0.000000000000000,
         0.000000000000000,  -0.000000000000000,   0.000000000000000,
        -0.000000000000000,   0.000000000000000,  -0.000000000000000,
         0.000000000000000,   1.000000000000000,   0.000000000000000,
        -0.000000000000000,   0.000000000000000,  -0.000000000000000,
         0.000000000000000,  -0.000000000000000,   0.000000000000000,
        -0.000000000000000,   0.000000000000000};
    const std::vector<double> hfiltI18 = {-0.000165438416514,  -0.000000000000000,
        -0.005942591074159,  -0.000000000000000,  -0.036400686183285,
        -0.000000000000000,  -0.138460275070318,  -0.000000000000000,
        -0.607805580229393,                   0,   0.607805580229393,
         0.000000000000000,   0.138460275070318,   0.000000000000000,
         0.036400686183285,   0.000000000000000,   0.005942591074159,
         0.000000000000000,   0.000165438416514};
    order = 18;
    try 
    {   
        hilbertFIR = FIR::HilbertTransformer(order, 8); 
    }   
    catch (const std::invalid_argument &ia)
    {   
        RTSEIS_ERRMSG("Failed hilbert length 9 design %s", ia.what()); 
        return EXIT_FAILURE;
    }
    realTaps = hilbertFIR.first.getFilterTaps();
    imagTaps = hilbertFIR.second.getFilterTaps();
    if (realTaps.size() != order + 1 || imagTaps.size() != order + 1)
    {   
        RTSEIS_ERRMSG("%s", "Hilbert is wrong length");
        return EXIT_FAILURE;
    }
    for (int i=0; i<order+1; i++)
    {   
        if (std::abs(hfiltR18[i] - realTaps[i]) > tol ||
            std::abs(hfiltI18[i] - imagTaps[i]) > tol)
        {   
            RTSEIS_ERRMSG("hilbert failed %d (%e,%e), (%e,%e)\n", i,
                           hfiltR18[i], realTaps[i], hfiltI18[i], imagTaps[i]);

        }   
    }

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
    try
    {
        hilbertFIR = FIR::HilbertTransformer(order, 8);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Failed hilbert length 10 design %s", ia.what()); 
        return EXIT_FAILURE;
    }
    realTaps = hilbertFIR.first.getFilterTaps();
    imagTaps = hilbertFIR.second.getFilterTaps();
    if (realTaps.size() != order + 1 || imagTaps.size() != order + 1)
    {
        RTSEIS_ERRMSG("%s", "Hilbert is wrong length");
        return EXIT_FAILURE;
    }
    for (int i=0; i<order+1; i++)
    {
        if (std::abs(hfiltR19[i] - realTaps[i]) > tol ||
            std::abs(hfiltI19[i] - imagTaps[i]) > tol)
        {
            RTSEIS_ERRMSG("hilbert failed %d (%e,%e), (%e,%e)\n", i,
                           hfiltR19[i], realTaps[i], hfiltI19[i], imagTaps[i]);

        }
    }
    return EXIT_SUCCESS;
}
