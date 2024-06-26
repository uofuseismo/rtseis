#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>
#include <stdexcept>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/postProcessing/singleChannel/waveform.hpp"
#include "rtseis/filterDesign/fir.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterImplementations/decimate.hpp"
#include "rtseis/filterImplementations/sosFilter.hpp"
#include "rtseis/filterImplementations/firFilter.hpp"
#include "rtseis/filterImplementations/iiriirFilter.hpp"
#include "rtseis/filterImplementations/iirFilter.hpp"

const std::string dataDir = "data/";
const std::string taperSolns100FileName = dataDir + "taper100.all.txt";
const std::string taperSolns101FileName = dataDir + "taper101.all.txt";
const std::string gse2FileName = dataDir + "gse2.txt";
const std::string firSolnsRef = dataDir + "firReference.txt";

using namespace RTSeis;
using namespace RTSeis::PostProcessing;
using namespace RTSeis::PostProcessing::SingleChannel;

int testDemean(void);
int testDetrend(void);
int testDownsample(const std::vector<double> &x);
int testDecimate(const std::vector<double> &x);
int testInterpolate(const std::vector<double> &x);
int testNormalization();
int testFilter(const std::vector<double> &x);
int testBandSpecificSOSFilters(const std::vector<double> &x);
int testBandSpecificIIRFilters(const std::vector<double> &x);
int testBandSpecificFIRFilters(const std::vector<double> &x);
int testTaper(void);
void readData(const std::string &fname, std::vector<double> &x);

int main(void)
{
    rtseis_utils_verbosity_setLoggingLevel(RTSEIS_SHOW_ALL);
    std::vector<double> gse2;
    readData(gse2FileName, gse2);
    if (gse2.size() < 1)
    {
        std::cerr << "No data" << std::endl;
        return EXIT_FAILURE;
    }
    int ierr;
    ierr = testDemean();
    if (ierr != EXIT_SUCCESS)
    { 
        std::cerr << "Failed demean test" << std::endl;
        return EXIT_FAILURE;
    } 
    std::cout << "Passed demean test" << std::endl;

    ierr = testDetrend();
    if (ierr != EXIT_SUCCESS)
    {
        std::cout << "Failed detrend test" << std::endl;
        return EXIT_FAILURE;
    } 
    std::cout << "Passed detrend test" << std::endl;

    ierr = testDownsample(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed downsample test" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Passed downsample test" << std::endl;

    ierr = testDecimate(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed decimate test" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Passed decimation test" << std::endl;

    ierr = testInterpolate(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed interp dft test" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Passed interp dft test" << std::endl;

    ierr = testFilter(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed filter test" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Passed generic filter test" << std::endl;

    ierr = testBandSpecificSOSFilters(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed sos band-specific tests" << std::endl;
        return EXIT_FAILURE;
    }
    std::cerr << "Passed SOS band-specific filter tests" << std::endl;

    ierr = testBandSpecificIIRFilters(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed iir band-specific tests" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Passed IIR band-specific filter tests" << std::endl;

    ierr = testBandSpecificFIRFilters(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed fir band-specific tests" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Passed FIR band-specific filter tests" << std::endl;

    ierr = testNormalization();
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed normalization tests" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "Passed normalization tests" << std::endl;

    ierr = testTaper();
    if (ierr != EXIT_SUCCESS)
    {
        std::cerr << "Failed taper test" << std::endl;
        return EXIT_FAILURE;
    }
    std::cerr << "Passed taper test" << std::endl;
    return EXIT_SUCCESS; 
}

//============================================================================//

int testFilter(const std::vector<double> &x)
{
    std::vector<double> yref;
    readData(firSolnsRef, yref);
    int len = static_cast<int> (yref.size());
    assert(len > 0);
    // Generate an FIR filter with SciPy: firwin(50, 0.4, 'hamming')
    int order = 50;
    const double fc = 0.4;
    std::vector<double> y;
    RTSeis::FilterRepresentations::FIR fir;
    try
    {
        fir = FilterDesign::FIR::FIR1Lowpass(order, fc,
                                             FilterDesign::FIRWindow::HAMMING);
        PostProcessing::SingleChannel::Waveform waveform; 
        waveform.setData(x);
        waveform.firFilter(fir);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    // Compare
    double l1Norm = 0;
    ippsNormDiff_L1_64f(yref.data(), y.data(), len, &l1Norm);
    if (l1Norm > 1.e-8)
    {
        std::cerr << "Failed to print fir filter" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

//============================================================================//

int testBandSpecificSOSFilters(const std::vector<double> &x)
{
    constexpr double tol = 1.e-10;
    double l1Norm;
    int npts = static_cast<int> (x.size());
    RTSeis::FilterRepresentations::SOS sos;
    RTSeis::FilterRepresentations::FIR fir;
    RTSeis::FilterRepresentations::BA ba;
    const double fsamp = 200;      // Sampling rate
    const double dt = 1/fsamp;     // Sampling period
    const double fnyq = 1./(2*dt); // Nyquist frequency
    std::vector<double> y;
    y.reserve(npts);
    // By this point the core filters work.  By this point we are verifying
    // that the higher order functions correctly call the low level library. 
    double fcV[2] = {0, 0};
    fcV[0] = 10/fnyq;
    sos = FilterDesign::IIR::designSOSIIRFilter(2, fcV, 5, 0, 
                      FilterDesign::Bandtype::LOWPASS,
                      FilterDesign::IIRPrototype::CHEBYSHEV1,
                      FilterDesign::IIRFilterDomain::DIGITAL);
    int ns = sos.getNumberOfSections();
    std::vector<double> bs = sos.getNumeratorCoefficients();
    std::vector<double> as = sos.getDenominatorCoefficients();
    FilterImplementations::SOSFilter
        <RTSeis::ProcessingMode::POST, double> sosFilt; 
    std::vector<double> ysosRef(npts);
    sosFilt.initialize(ns, bs.data(), as.data());
    double *yptr = ysosRef.data();
    sosFilt.apply(npts, x.data(), &yptr);
    sosFilt.clear();
    // SOS filter
//! [ppSCSOSLowpassExample]
    PostProcessing::SingleChannel::Waveform waveform;
    waveform.setSamplingPeriod(dt); // Sampling period is 1/200
    try
    {
        // Design a 2nd order Chebyshev filter with 10 Hz cutoff
        int order = 2;           // 2nd order (2 poles)
        double fc = 10;          // 10 Hz cutoff
        double ripple = 5;       // 5 dB ripple in passband
        bool lzeroPhase = false; // Apply causally (will retain nonlinear phase response)
        waveform.setData(x);
        waveform.sosLowpassFilter(order, fc,
                                  SingleChannel::IIRPrototype::CHEBYSHEV1,
                                  ripple, lzeroPhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
//! [ppSCSOSLowpassExample]
    // Compare
    ippsNormDiff_L1_64f(ysosRef.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed sos filter test with error="
                  <<  l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                            Highpass filter                             //
    //------------------------------------------------------------------------//
    fcV[0] = 10/fnyq;
    sos = FilterDesign::IIR::designSOSIIRFilter(2, fcV, 0, 0,
                      FilterDesign::Bandtype::HIGHPASS,
                      FilterDesign::IIRPrototype::BESSEL,
                      FilterDesign::IIRFilterDomain::DIGITAL);
    ns = sos.getNumberOfSections();
    bs = sos.getNumeratorCoefficients();
    as = sos.getDenominatorCoefficients();
    std::vector<double> ytemp(npts);
    sosFilt.initialize(ns, bs.data(), as.data());
    yptr = ytemp.data();
    sosFilt.apply(npts, x.data(), &yptr);
    std::reverse(ytemp.begin(), ytemp.end());
    yptr = ysosRef.data();
    sosFilt.apply(npts, ytemp.data(), &yptr); //ysosRef.data());
    std::reverse(ysosRef.begin(), ysosRef.end());
    sosFilt.clear();
    // Hammer on the filter designer
    for (int j=0; j<5; j++)
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int order = 2;           // 2nd order (2 poles)
        double fc = 10;          // 10 Hz cutoff
        double ripple = 0;       // N/A
        bool lzeroPhase = true;  // Zero-phase
        waveform.setData(x);
        waveform.sosHighpassFilter(order, fc,
                                   SingleChannel::IIRPrototype::BESSEL,
                                   ripple, lzeroPhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ysosRef.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed sos filter test with error="
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                              Bandpass filter                           //
    //------------------------------------------------------------------------//
    fcV[0] = 1/fnyq;
    fcV[1] = 10/fnyq;
    sos = FilterDesign::IIR::designSOSIIRFilter(2, fcV, 0, 0,
                      FilterDesign::Bandtype::BANDPASS,
                      FilterDesign::IIRPrototype::BUTTERWORTH,
                      FilterDesign::IIRFilterDomain::DIGITAL);
    ns = sos.getNumberOfSections();
    bs = sos.getNumeratorCoefficients();
    as = sos.getDenominatorCoefficients();
    sosFilt.initialize(ns, bs.data(), as.data());
    yptr = ytemp.data();
    sosFilt.apply(npts, x.data(), &yptr); //ytemp.data());
    std::reverse(ytemp.begin(), ytemp.end());
    yptr = ysosRef.data();
    sosFilt.apply(npts, ytemp.data(), &yptr); //ysosRef.data());
    std::reverse(ysosRef.begin(), ysosRef.end());
    sosFilt.clear();
    for (int j=0; j<1; j++)
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int order = 2;                     // 2nd order (becomes 4 poles)
        std::pair<double,double> fc(1,10); // Passband is 1 to 10 Hz
        double ripple = 0;                 // N/A
        bool lzeroPhase = true;            // Zero-phase
        waveform.setData(x);
        waveform.sosBandpassFilter(order, fc,
                                   SingleChannel::IIRPrototype::BUTTERWORTH,
                                   ripple, lzeroPhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ysosRef.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed sos bp filter test with error="
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                              Bandstop filter                           //
    //------------------------------------------------------------------------//
    fcV[0] = 1/fnyq;
    fcV[1] = 10/fnyq;
    sos = FilterDesign::IIR::designSOSIIRFilter(3, fcV, 0, 0,
                          FilterDesign::Bandtype::BANDSTOP,
                          FilterDesign::IIRPrototype::BESSEL,
                          FilterDesign::IIRFilterDomain::DIGITAL);
    ns = sos.getNumberOfSections();
    bs = sos.getNumeratorCoefficients();
    as = sos.getDenominatorCoefficients();
    sosFilt.initialize(ns, bs.data(), as.data());
    yptr = ysosRef.data();
    sosFilt.apply(npts, x.data(), &yptr); //ysosRef.data());
    sosFilt.clear();
    for (int j=0; j<1; j++)
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int order = 3;                     // 2nd order (becomes 4 poles)
        std::pair<double,double> fc(1,10); // Stopband is 1 to 10 Hz
        double ripple = 0;                 // N/A
        bool lzeroPhase = false;           // Not Zero-phase
        waveform.setData(x);
        waveform.sosBandstopFilter(order, fc,
                                   SingleChannel::IIRPrototype::BESSEL,
                                   ripple, lzeroPhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ysosRef.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed sos bs filter test with error=" 
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS; 
}

//============================================================================//

int testBandSpecificIIRFilters(const std::vector<double> &x)
{
    constexpr double tol = 1.e-10;
    double l1Norm;
    int npts = static_cast<int> (x.size());
    RTSeis::FilterRepresentations::BA ba;
    const double fsamp = 200;      // Sampling rate
    const double dt = 1/fsamp;     // Sampling period
    const double fnyq = 1./(2*dt); // Nyquist frequency
    std::vector<double> y;
    y.reserve(npts);
    // By this point the core filters work.  By this point we are verifying
    // that the higher order functions correctly call the low level library. 
    double fcV[2] = {0, 0};
    fcV[0] = 5/fnyq;
    ba = FilterDesign::IIR::designBAIIRFilter(2, fcV, 5, 0,
                      FilterDesign::Bandtype::LOWPASS,
                      FilterDesign::IIRPrototype::CHEBYSHEV1,
                      FilterDesign::IIRFilterDomain::DIGITAL);
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients();
    int nb = static_cast<int> (b.size());
    int na = static_cast<int> (a.size());
    FilterImplementations::IIRFilter<ProcessingMode::POST, double> iirFilt; 
    FilterImplementations::IIRIIRFilter<double> iiriirFilt;
    std::vector<double> yiirRef(npts);
    iirFilt.initialize(nb, b.data(), na, a.data());
    double *yptr = yiirRef.data();
    iirFilt.apply(npts, x.data(), &yptr); //yiirRef.data());
    iirFilt.clear();
    // IIR filter
//! [ppSCIIRLowpassExample]
    PostProcessing::SingleChannel::Waveform waveform;
    waveform.setSamplingPeriod(dt); // Sampling period is 1/200
    try
    {
        // Design a 2nd order Chebyshev filter with 10 Hz cutoff
        int order = 2;           // 2nd order (2 poles)
        double fc = 5;           // 5 Hz cutoff
        double ripple = 5;       // 5 dB ripple in passband
        bool lzeroPhase = false; // Apply causally (will retain nonlinear phase response)
        waveform.setData(x);
        waveform.iirLowpassFilter(order, fc,
                                  SingleChannel::IIRPrototype::CHEBYSHEV1,
                                  ripple, lzeroPhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
//! [ppSCIIRLowpassExample]
    // Compare
    ippsNormDiff_L1_64f(yiirRef.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed iir filter test with error="
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                            Highpass filter                             //
    //------------------------------------------------------------------------//
    fcV[0] = 15/fnyq;
    ba = FilterDesign::IIR::designBAIIRFilter(2, fcV, 0, 0,
                      FilterDesign::Bandtype::HIGHPASS,
                      FilterDesign::IIRPrototype::BESSEL,
                      FilterDesign::IIRFilterDomain::DIGITAL);
    b = ba.getNumeratorCoefficients();
    a = ba.getDenominatorCoefficients();
    nb = static_cast<int> (b.size());
    na = static_cast<int> (a.size());
    std::vector<double> ytemp(npts);
    iiriirFilt.initialize(nb, b.data(), na, a.data());
    yptr = ytemp.data();
    iiriirFilt.apply(npts, x.data(), &yptr); //ytemp.data());
    iiriirFilt.clear();
    // Hammer on the filter designer
    for (int j=0; j<5; j++)
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int order = 2;           // 2nd order (2 poles)
        double fc = 15;          // 15 Hz cutoff
        double ripple = 0;       // N/A
        bool lzeroPhase = true;  // Zero-phase
        waveform.setData(x);
        waveform.iirHighpassFilter(order, fc,
                                   SingleChannel::IIRPrototype::BESSEL,
                                   ripple, lzeroPhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ytemp.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed iiriir filter test with error="
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                              Bandpass filter                           //
    //------------------------------------------------------------------------//
    fcV[0] = 0.5/fnyq;
    fcV[1] = 8.5/fnyq;
    ba = FilterDesign::IIR::designBAIIRFilter(2, fcV, 0, 0,
                      FilterDesign::Bandtype::BANDPASS,
                      FilterDesign::IIRPrototype::BUTTERWORTH,
                      FilterDesign::IIRFilterDomain::DIGITAL);
    b = ba.getNumeratorCoefficients();
    a = ba.getDenominatorCoefficients();
    nb = static_cast<int> (b.size());
    na = static_cast<int> (a.size());
    iiriirFilt.initialize(nb, b.data(), na, a.data());
    yptr = ytemp.data();
    iiriirFilt.apply(npts, x.data(), &yptr); //ytemp.data());
    iiriirFilt.clear();
    for (int j=0; j<1; j++)
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int order = 2;                        // 2nd order (becomes 4 poles)
        std::pair<double,double> fc(0.5,8.5); // Passband is 0.5 to 8.5 Hz
        double ripple = 0;                    // N/A
        bool lzeroPhase = true;               // Zero-phase
        waveform.setData(x);
        waveform.iirBandpassFilter(order, fc,
                                   SingleChannel::IIRPrototype::BUTTERWORTH,
                                   ripple, lzeroPhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ytemp.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed iir bp filter test with error="
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                              Bandstop filter                           //
    //------------------------------------------------------------------------//
    fcV[0] = 3/fnyq;
    fcV[1] = 15/fnyq;
    ba = FilterDesign::IIR::designBAIIRFilter(3, fcV, 0, 0,
                          FilterDesign::Bandtype::BANDSTOP,
                          FilterDesign::IIRPrototype::BESSEL,
                          FilterDesign::IIRFilterDomain::DIGITAL);
    b = ba.getNumeratorCoefficients();
    a = ba.getDenominatorCoefficients();
    nb = static_cast<int> (b.size());
    na = static_cast<int> (a.size());
    iirFilt.initialize(nb, b.data(), na, a.data());
    yptr = ytemp.data();
    iirFilt.apply(npts, x.data(), &yptr); //ytemp.data());
    iirFilt.clear();
    for (int j=0; j<1; j++)
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int order = 3;                     // 2nd order (becomes 4 poles)
        std::pair<double,double> fc(3,15); // Stopband is 3 to 15 Hz
        double ripple = 0;                 // N/A
        bool lzeroPhase = false;           // Not Zero-phase
        waveform.setData(x);
        waveform.iirBandstopFilter(order, fc,
                                   SingleChannel::IIRPrototype::BESSEL,
                                   ripple, lzeroPhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ytemp.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed sos bs filter test with error="
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS; 
}

//===========================================================================//

int testBandSpecificFIRFilters(const std::vector<double> &x)
{
    constexpr double tol = 1.e-10;
    double l1Norm;
    int npts = static_cast<int> (x.size());
    RTSeis::FilterRepresentations::FIR fir;
    const double fsamp = 200;      // Sampling rate
    const double dt = 1/fsamp;     // Sampling period
    const double fnyq = 1./(2*dt); // Nyquist frequency
    std::vector<double> y;
    y.reserve(npts);
    // By this point the core filters work.  By this point we are verifying
    // that the higher order functions correctly call the low level library. 
    double fcd = 10/fnyq;
    int order  = 100;
    fir = FilterDesign::FIR::FIR1Lowpass(order, fcd,
                                         FilterDesign::FIRWindow::HAMMING);
    int ntaps = fir.getNumberOfFilterTaps();
    std::vector<double> b = fir.getFilterTaps();
    FilterImplementations::FIRFilter<RTSeis::ProcessingMode::POST, double> firFilt;
    std::vector<double> yfirRef(npts);
    firFilt.initialize(ntaps, b.data());
    double *yptr = yfirRef.data();
    firFilt.apply(npts, x.data(), &yptr); //yfirRef.data());
    firFilt.clear();
    // FIR filter
//! [ppSCFIRLowpassExample]
    PostProcessing::SingleChannel::Waveform waveform;
    waveform.setSamplingPeriod(dt); // Sampling period is 1/200
    try
    {
        // Design a 2nd order Chebyshev filter with 10 Hz cutoff
        int ntaps = 101;           // Number of filter taps
        double fc = 10;            // 10 Hz cutoff
        bool lremovePhase = false; // Keep linear phase shift
        waveform.setData(x);
        waveform.firLowpassFilter(ntaps, fc,
                                  SingleChannel::FIRWindow::HAMMING,
                                  lremovePhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
//! [ppSCFIRLowpassExample]
    // Compare
    ippsNormDiff_L1_64f(yfirRef.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed filter test with error=" << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                            Highpass filter                             //
    //------------------------------------------------------------------------//
    fcd = 5/fnyq;
    order  = 90;
    int nextra = (order+1)/2;
    fir = FilterDesign::FIR::FIR1Highpass(order, fcd,
                                          FilterDesign::FIRWindow::BLACKMAN_OPT);
    ntaps = fir.getNumberOfFilterTaps();
    b = fir.getFilterTaps();
    std::vector<double> xpad = x;
    xpad.resize(x.size()+nextra, 0);
    std::vector<double> ytemp(xpad.size());
    firFilt.initialize(ntaps, b.data());
    yptr = ytemp.data();
    firFilt.apply(npts+nextra, xpad.data(), &yptr);
    firFilt.clear();
    // Hammer on the filter designer
    for (int j=0; j<5; j++)
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int ntaps = order + 1;     // Length of FIR filter
        double fc = 5;             // 5 Hz cutoff
        bool lremovePhase = true;  // Remove phase shift
        waveform.setData(x);
        waveform.firHighpassFilter(ntaps, fc,
                                   SingleChannel::FIRWindow::BLACKMAN_OPT,
                                   lremovePhase);
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    assert(static_cast<int> (y.size()) == npts); 
    double *ycomp = ytemp.data();
    ippsNormDiff_L1_64f(&ycomp[nextra], y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed fir filter test with error=" << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                              Bandpass filter                           //
    //------------------------------------------------------------------------//
    std::pair<double,double> fcbp(2/fnyq, 10/fnyq);
    order  = 150;
    nextra = (order+1)/2;
    fir = FilterDesign::FIR::FIR1Bandpass(order, fcbp,
                                          FilterDesign::FIRWindow::HANN);
    ntaps = fir.getNumberOfFilterTaps();
    b = fir.getFilterTaps();
    firFilt.initialize(ntaps, b.data());
    xpad = x; 
    xpad.resize(x.size()+nextra, 0);
    ytemp.resize(xpad.size());
    yptr = ytemp.data();
    firFilt.apply(npts+nextra, xpad.data(), &yptr);
    firFilt.clear();
    // Hammer on the filter designer
    for (int j=0; j<5; j++) 
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int ntaps = order + 1;              // Length of FIR filter
        std::pair<double,double> fc(2, 10); // 2-10 Hz bandpass 
        bool lremovePhase = true;           // Remove phase shift
        waveform.setDataPointer(x.size(), x.data());
        waveform.firBandpassFilter(ntaps, fc,
                                   SingleChannel::FIRWindow::HANN,
                                   lremovePhase);
        waveform.releaseDataPointer();
        y = waveform.getData();
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    assert(static_cast<int> (y.size()) == npts);
    ycomp = ytemp.data();
    ippsNormDiff_L1_64f(&ycomp[nextra], y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {    
        std::cerr << "Failed fir filter test with error="
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    //                              Bandstop filter                           //
    //------------------------------------------------------------------------//
    std::pair<double,double> fcbs = std::make_pair(1/fnyq, 11/fnyq);
    order  = 200; 
    fir = FilterDesign::FIR::FIR1Bandstop(order, fcbs,
                                          FilterDesign::FIRWindow::BARTLETT);
    ntaps = fir.getNumberOfFilterTaps();
    b = fir.getFilterTaps();
    firFilt.initialize(ntaps, b.data());
    yptr = ytemp.data();
    firFilt.apply(npts, x.data(), &yptr);
    firFilt.clear();
    // Hammer on the filter designer
    for (int j=0; j<5; j++) 
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int ntaps = order + 1;              // Length of FIR filter
        std::pair<double,double> fc(1, 11); // 2-10 Hz bandpass 
        bool lremovePhase = false;          // Remove phase shift
        waveform.setData(x);
        waveform.firBandstopFilter(ntaps, fc,
                                   SingleChannel::FIRWindow::BARTLETT,
                                   lremovePhase);
        y = waveform.getData();
    }    
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    assert(static_cast<int> (y.size()) == npts);
    ippsNormDiff_L1_64f(ytemp.data(), y.data(), npts, &l1Norm);
    if (l1Norm > tol)
    {
        std::cerr << "Failed fir filter test with error="
                  << l1Norm << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS; 
}

//============================================================================//

int testDownsample(const std::vector<double> &x)
{
    const int nq = 7;
    int npts = static_cast<int> (x.size());
    const double dt = 1.0/200.0;
    PostProcessing::SingleChannel::Waveform waveform;
    waveform.setSamplingPeriod(dt);
    for (int iq=1; iq<nq+1; iq++)
    {
        std::vector<double> y;
        try
        {
            waveform.setData(x);
            waveform.downsample(iq);
            y = waveform.getData();
        }
        catch (std::exception &e)
        {
            std::cerr << "Error in downsampling "
                      << iq << " " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
        // Verify
        int j = 0; 
        for (int i=0; i<npts; i=i+iq)
        {
            if (std::abs(x[i] - y.at(j)) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed downsample %d %d %lf %lf", 
                              iq, i, x[i], y[j]); 
                return EXIT_FAILURE;
            }
            j = j + 1;
        }
        if (j != static_cast<int> (y.size()))
        {
            std::cerr << "Size mismatch" << std::endl;
            return EXIT_FAILURE;
        }
        double dtNew = waveform.getSamplingPeriod();
        double dtTarg = static_cast<double> (iq)*dt;
        if (std::abs(dtNew - dtTarg) > 1.e-12)
        {
            std::cerr << "Failed to update sampling period "
                      << dtTarg << " " << dtNew << std::endl;
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int testDecimate(const std::vector<double> &x)
{
    const int nq = 7;
    bool removePhaseShift = true;
    int npts = static_cast<int> (x.size());
    int firLen = 33;
    const double dt = 1.0/200.0;
    PostProcessing::SingleChannel::Waveform waveform;
    waveform.setSamplingPeriod(dt);
    for (int iq=2; iq<nq+1; iq++)
    {
        std::vector<double> y;
        try
        {
            waveform.setData(x);
            waveform.decimate(iq, firLen);
            y = waveform.getData();
        }
        catch (std::exception &e)
        {
            std::cerr << "Error in downsampling "
                      << iq << " " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
        // Verify
        FilterImplementations::Decimate
          <RTSeis::ProcessingMode::POST, double> decim;
        decim.initialize(iq, firLen, removePhaseShift);
        std::vector<double> yref(npts);
        double *yRefPtr = yref.data();
        int lenRef;
        decim.apply(x.size(), x.data(), npts, &lenRef, &yRefPtr); 
        if (lenRef != static_cast<int> (y.size()))
        {
             std::cerr << "Inconsistent sizes" << std::endl;
             return EXIT_FAILURE;
        }
        double error;
        ippsNormDiff_L1_64f(yref.data(), y.data(), lenRef, &error);
        if (error > 1.e-12)
        {
            std::cerr << "decimation failed with error " << error << std::endl;
            return EXIT_FAILURE;
        }
        double dtNew = waveform.getSamplingPeriod();
        double dtTarg = static_cast<double> (iq)*dt;
        if (std::abs(dtNew - dtTarg) > 1.e-12)
        {
            std::cerr << "Failed to update sampling period " << dtTarg 
                      << " " << dtNew << std::endl;
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int testInterpolate(const std::vector<double> &x)
{
    double df = 200;
    std::vector<double> dfNews({200, 300, 400, 450, 500});
    // This is computed: 
    //   destStride = lcm/df 
    // where lcm is the least common multiple defined by:
    // (df*dfNew)/(GCD(df, dfNew)) where GCD is the greatest common denominator.
    std::vector<int> srcStrides({1, 2, 1, 4, 2}); 
    std::vector<int> destStrides({1, 3, 2, 9, 5}); 
    for (int k=0; k<static_cast<int> (dfNews.size()); ++k)
    {
        auto dfNew = dfNews[k];
        double dt = 1.0/df;
        double dtNew = 1.0/dfNew;
        PostProcessing::SingleChannel::Waveform waveform;
        waveform.setSamplingPeriod(dt);
        waveform.setData(x); 
        waveform.detrend(); // Remove trend
        auto xDetrend = waveform.getData(); // Get detrended data
        waveform.interpolate(dtNew, InterpolationMethod::DFT);
        auto yint = waveform.getData(); // Get resampled data

        double error = 0;
        int ic = 0;
        for (auto i=0; i<static_cast<int> (xDetrend.size()); i=i+srcStrides[k])
        {
            auto j = ic*destStrides[k];
            error = std::max(error, std::abs(yint[j] - x[i]) );
            ic = ic + 1;
        }
        // Basically on the same order b/c max is 1000 
        if (error > 1.e-1)
        {
            std::cerr << "Failed interfpt on iteration "
                      << ic << " " << error << std::endl;
            return EXIT_FAILURE; 
        }
     }
     return EXIT_SUCCESS;
}

//============================================================================//

int testDemean()
{
    // Do two passes
    for (auto j=0; j<2; j++)
    {
    // Create a reference signal from which to remove the mean
    std::vector<double> x;
    int npts = 5000 + j;
    x.resize(npts);
    double xmean = 5.1 + static_cast<double> (j);
    std::fill(x.begin(), x.end(), xmean);
    std::vector<double> y;
//! [ppSCDemeanExample]
    PostProcessing::SingleChannel::Waveform waveform;
    try
    {
        waveform.setData(x);    // Set waveform data to demean
        waveform.demean();      // Demean the data
        y = waveform.getData(); // Get the demeaned data in a vector, y
    }
    catch (std::runtime_error &rt)
    {
        fprintf(stderr, "Demean failed %s\n", rt.what()); 
        return EXIT_FAILURE;
    }
//! [ppSCDemeanExample]
    // Verify
    if (y.size() != x.size())
    {
        std::cerr << "Inconsistent sizes" << std::endl;
        return EXIT_FAILURE;
    }
    double maxAbs = 0;
    ippsMaxAbs_64f(y.data(), npts, &maxAbs);
    if (maxAbs > 1.e-13)
    {
        std::cerr << "Demean failed " << maxAbs << std::endl;
        return EXIT_FAILURE;
    }
    } // Loop on passes

    return EXIT_SUCCESS;
}
//============================================================================//
int testDetrend()
{
    for (auto j=0; j<2; j++)
    {
    std::vector<double> x;
    std::vector<double> y;
    int npts = 50001 + j;
    x.resize(npts);
    #pragma omp simd
    for (int i=0; i<npts; i++){x[i] = 1.1 + 0.3*static_cast<double> (i);}
//! [ppSCDetrendExample]
    PostProcessing::SingleChannel::Waveform waveform;
    try 
    {
        waveform.setData(x);     // Set the data to detrend
        waveform.detrend();      // Detrend
        y = waveform.getData();  // Get the detrended data in a vector, y
    }   
    catch (std::runtime_error &rt)
    {
        fprintf(stderr, "Detrend failed %s\n", rt.what()); 
        return EXIT_FAILURE;
    }   
//! [ppSCDetrendExample]
    double maxAbs = 0;
    ippsMaxAbs_64f(y.data(), npts, &maxAbs);
    // ~50,000 points loses about 5 digits which is sensible 
    if (maxAbs > 1.e-8)
    {
        std::cerr << "Demean failed " << maxAbs << std::endl;
        return EXIT_FAILURE;
    }
    }; // Loop on npts
    return EXIT_SUCCESS;
}
//============================================================================//
int testNormalization()
{
    /// Normalization
    int npts = 1001;
    std::vector<double> x(npts), y;
    // Put in range [-1,2] then compute mean and stanard deviation
    double mean = 0;
    for (auto i=0; i<npts; ++i)
    {
        x[i] = -1.0 + 3*static_cast<double> (rand())/RAND_MAX;
        if (i == 0){x[i] =-1;}
        if (i == 1){x[i] = 2;}
        mean = mean + x[i];
    }
    mean = mean/static_cast<double> (npts);
    double var = 0;
    for (auto i=0; i<npts; ++i)
    {
        var = var + std::pow(x[i] - mean, 2);
    }
    auto std = std::sqrt(var/static_cast<double> (npts - 1));
    // sign-bit normalization
    PostProcessing::SingleChannel::Waveform waveform;
    try
    {
        waveform.setData(x);
        waveform.normalizeSignBit();
        y = waveform.getData();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    double error = 0;
    for (auto i=0; i<npts; ++i)
    {
        double yest =-1;
        if (x[i] >= 0){yest = 1;}
        error = std::max(error, std::abs(yest - y[i]));
    }
    if (error > 1.e-14)
    {
        std::cerr << "signBit normalization failed " << error << std::endl;
        return EXIT_FAILURE;
    }
    // z-score normalization
    try
    {
        waveform.setData(x);
        waveform.normalizeZScore();
        y = waveform.getData();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    error = 0;
    for (auto i=0; i<npts; ++i)
    {
        error = std::max(error, std::abs(y[i] - (x[i] - mean)/std));
    }
    if (error >= 1.e-14)
    {
        std::cerr << "zscore normalization failed: " << error << std::endl;
        return EXIT_FAILURE;
    }
    // transform to [-1,1]
    std::pair<double, double> targetRange(-1, 1);
    try
    {
        waveform.setData(x);
        waveform.normalizeMinMax(targetRange);
        y = waveform.getData();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    error = 0;
    for (auto i=0; i<npts; ++i)
    {
        auto yp = -1.0 + (1.0 - -1.0)*(x[i] - -1.0)/(2.0 - -1.0);
        error = std::max(error, std::abs(y[i] - yp));
    }
    if (error >= 1.e-14)
    {
        std::cerr << "minMax normalization failed: " << error << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
//============================================================================//
int testTaper()
{
    // Load the reference solutions
    std::vector<double> yHamming100Ref(100);
    std::vector<double> yHanning100Ref(100);
    std::vector<double> ySine100Ref(100);
    std::vector<double> yHamming101Ref(101);
    std::vector<double> yHanning101Ref(101);
    std::vector<double> ySine101Ref(101);
    std::string line100, line101;
    std::ifstream taper100File(taperSolns100FileName);
    std::ifstream taper101File(taperSolns101FileName);
    for (int i=0; i<100; i++)
    {
        if (!std::getline(taper100File, line100))
        {
            std::cerr << "Premature end of file " << taperSolns100FileName
                      << std::endl;
            return EXIT_FAILURE;
        }
        std::sscanf(line100.c_str(), "%lf, %lf, %lf\n",
                    &yHamming100Ref[i], &yHanning100Ref[i], &ySine100Ref[i]);
        if (!std::getline(taper101File, line101))
        {
            std::cerr << "Premature end of file "
                      << taperSolns101FileName << std::endl;
            return EXIT_FAILURE;
        }
        std::sscanf(line101.c_str(), "%lf, %lf, %lf\n",
                    &yHamming101Ref[i], &yHanning101Ref[i], &ySine101Ref[i]);
    }
    if (!std::getline(taper101File, line101))
    {
        std::cerr << "Premature end of file "
                  << taperSolns101FileName << std::endl;
        return EXIT_FAILURE;
    }
    std::sscanf(line101.c_str(), "%lf, %lf, %lf\n",
                &yHamming101Ref[100], &yHanning101Ref[100], &ySine101Ref[100]);

    taper100File.close();
    taper101File.close();

    namespace SC = PostProcessing::SingleChannel;
    PostProcessing::SingleChannel::Waveform waveform;
    std::vector<double> x100(100, 1);
    std::vector<double> yHamming100;
    std::vector<double> yHann100;
    std::vector<double> ySine100;
    try
    {
        double pct = 100*(0.2*2); // SAC to RTSeis; 40 pct of signal
        waveform.setData(x100);
        waveform.taper(pct, SC::TaperParameters::HAMMING);
        yHamming100 = waveform.getData();

        pct = 100*(0.1*2); // SAC to RTSeis; 20 pct of signal
        waveform.setData(x100);
        waveform.taper(pct, SC::TaperParameters::HANN);
        yHann100 = waveform.getData();

        pct = 100*(0.3*2); // SAC to RTSeis; 60 pct of signal
        waveform.setData(x100);
        waveform.taper(pct, SC::TaperParameters::SINE);
        ySine100 = waveform.getData(); 
    }
    catch (const std::invalid_argument &ia)
    {
        std::cerr << "Taper 100 failed " << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
    // Compare
    for (size_t i=0; i<ySine100.size(); i++)
    {
        if (std::abs(yHamming100[i] - yHamming100Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Hamming 100 failed %ld,%lf,%lf", 
                          i, yHamming100[i], yHamming100Ref[i])
            return EXIT_FAILURE;
        }
        if (std::abs(yHann100[i] - yHanning100Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Hann 100 failed %ld,%lf,%lf", 
                          i, yHann100[i], yHanning100Ref[i])
            return EXIT_FAILURE;
        }
        if (std::abs(ySine100[i] - ySine100Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Sine 100 failed %ld,%lf,%lf", 
                          i, ySine100[i], ySine100Ref[i])
            return EXIT_FAILURE;
        }
    }
    // Repeat for 101 points
    std::vector<double> x101(101, 1); 
    std::vector<double> yHamming101;
    std::vector<double> yHann101;
    std::vector<double> ySine101;
    try 
    {   
        double pct = 100*(0.05*2); // SAC to RTSeis; 10 pct of signal
        waveform.setData(x101);
        waveform.taper(pct, SC::TaperParameters::HAMMING);
        yHamming101 = waveform.getData();

        pct = 100*(0.1*2); // SAC to RTSeis; 20 pct of signal
        waveform.setData(x101);
        waveform.taper(pct, SC::TaperParameters::HANN);
        yHann101 = waveform.getData();

        pct = 100*(0.15*2); // SAC to RTSeis; 30 pct of signal
        waveform.setData(x101);
        waveform.taper(pct, SC::TaperParameters::SINE);
        ySine101 = waveform.getData(); 
    }   
    catch (const std::invalid_argument &ia)
    {   
        std::cerr << "Taper 101 failed " << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
    for (size_t i=0; i<101; i++)
    {
        if (std::abs(yHamming101[i] - yHamming101Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Hamming 101 failed %ld,%lf,%lf", 
                          i, yHamming101[i], yHamming101Ref[i])
            return EXIT_FAILURE;
        }
        if (std::abs(yHann101[i] - yHanning101Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Hann 101 failed %ld,%lf,%lf", 
                          i, yHann101[i], yHanning101Ref[i])
            return EXIT_FAILURE;
        }   
        if (std::abs(ySine101[i] - ySine101Ref[i]) > 1.e-6)
        {   
            RTSEIS_ERRMSG("Sine 101 failed %ld,%lf,%lf", 
                          i, ySine101[i], ySine101Ref[i])
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

void readData(const std::string &fname, std::vector<double> &x)
{
    x.reserve(12000);
    x.resize(0);
    std::string line;
    std::ifstream textFile(fname);
    while (std::getline(textFile, line))
    {
        double xval;
        std::sscanf(line.c_str(), "%lf\n", &xval);
        x.push_back(xval);
    }
    return;
}
