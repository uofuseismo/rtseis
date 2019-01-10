#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include <cmath>
#include "rtseis/utils/design.hpp"
#include "rtseis/log.h"
#include "utils.hpp"

using namespace RTSeis::Utils::FilterDesign;

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
    std::vector<double> ptaps;
    // Test 1
    order = 13; r = 0.6; 
    int ierr = FIR::FIR1Lowpass(order, r, ptaps, FIR::Window::HAMMING);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to generate lowpass filter");
        return EXIT_FAILURE;
    }
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
    ierr = FIR::FIR1Highpass(order, r, ptaps, FIR::Window::HANN);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to generate highpass filter");
        return EXIT_FAILURE;
    }
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
    order = 11; double r2[2] = {0.2, 0.8};
    ierr = FIR::FIR1Bandpass(order, r2, ptaps, FIR::Window::BARTLETT); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to generate bandpass filter");
        return EXIT_FAILURE;
    }
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
    r2[0] = 0.15; r2[1] = 0.85;
    ierr = FIR::FIR1Bandstop(order, r2, ptaps, FIR::Window::BLACKMAN_OPT);
    if (ierr != 0)  
    {
        RTSEIS_ERRMSG("%s", "Failed to generate bandstop filter");
        return EXIT_FAILURE;
    }
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
    return EXIT_SUCCESS;
}
