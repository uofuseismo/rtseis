#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/postProcessing/singleChannel/waveform.hpp"
#include "rtseis/utilities/design/fir.hpp"
#include "rtseis/utilities/design/iir.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterRepresentations/sos.hpp"
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/utilities/filterImplementations/sosFilter.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#include "rtseis/utilities/filterImplementations/iiriirFilter.hpp"
#include "rtseis/utilities/filterImplementations/iirFilter.hpp"

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
int testFilter(const std::vector<double> &x);
int testBandSpecificFilters(const std::vector<double> &x);
int testTaper(void);
void readData(const std::string &fname, std::vector<double> &x);

int main(void)
{
    rtseis_utils_verbosity_setLoggingLevel(RTSEIS_SHOW_ALL);
    std::vector<double> gse2;
    readData(gse2FileName, gse2);
    if (gse2.size() < 1)
    {
        RTSEIS_ERRMSG("%s", "No data");
        return EXIT_FAILURE;
    }
    int ierr;
    ierr = testDemean();
    if (ierr != EXIT_SUCCESS)
    { 
        RTSEIS_ERRMSG("%s", "Failed demean test");
        return EXIT_FAILURE;
    } 
    RTSEIS_INFOMSG("%s", "Passed demean test");

    ierr = testDetrend();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed detrend test");
        return EXIT_FAILURE;
    } 
    RTSEIS_INFOMSG("%s", "Passed detrend test");

    ierr = testFilter(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed detrend test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed generic filter test");

    ierr = testBandSpecificFilters(gse2);
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed band specific filter tests");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed band specific filter tests");

    ierr = testTaper();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed taper test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed window test");
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
    Utilities::FilterRepresentations::FIR fir;
    try
    {
        Utilities::FilterDesign::FIR::FIR1Lowpass(order, fc, fir,
                                  Utilities::FilterDesign::FIRWindow::HAMMING);
        PostProcessing::SingleChannel::Waveform waveform; 
        waveform.setData(x);
        waveform.firFilter(fir);
        waveform.getData(y);     
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
        RTSEIS_ERRMSG("%s", "Failed to print fir filter");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int testBandSpecificFilters(const std::vector<double> &x)
{
    double l1Norm;
    int npts = static_cast<int> (x.size());
    Utilities::FilterRepresentations::SOS sos;
    Utilities::FilterRepresentations::FIR fir;
    Utilities::FilterRepresentations::BA ba;
    const double fsamp = 200;      // Sampling rate
    const double dt = 1/fsamp;     // Sampling period
    const double fnyq = 1./(2*dt); // Nyquist frequency
    std::vector<double> y;
    y.reserve(npts);
    // By this point the core filters work.  By this point we are verifying
    // that the higher order functions correctly call the low level library. 
    double fcV[2] = {0, 0};
    fcV[0] = 10/fnyq;
    Utilities::FilterDesign::IIR::iirfilter(2, fcV, 5, 0, 
                      Utilities::FilterDesign::Bandtype::LOWPASS,
                      Utilities::FilterDesign::IIRPrototype::CHEBYSHEV1,
                      sos,
                      Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    int ns = sos.getNumberOfSections();
    std::vector<double> bs = sos.getNumeratorCoefficients();
    std::vector<double> as = sos.getDenominatorCoefficients();
    Utilities::FilterImplementations::SOSFilter sosFilt; 
    std::vector<double> ysosRef(npts);
    sosFilt.initialize(ns, bs.data(), as.data(),
                       RTSeis::ProcessingMode::POST_PROCESSING,
                       RTSeis::Precision::DOUBLE);
    sosFilt.apply(npts, x.data(), ysosRef.data());
    sosFilt.clear();
    // Zero-phase SOS filter
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
        waveform.getData(y);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
//! [ppSCSOSLowpassExample]
    // Compare
    ippsNormDiff_L1_64f(ysosRef.data(), y.data(), npts, &l1Norm);
    if (l1Norm > 1.e-12)
    {
        fprintf(stderr, "Failed sos filter test with error=%e\n", l1Norm);
        return EXIT_FAILURE;
    }
    // Highpass filter
    fcV[0] = 10/fnyq;
    Utilities::FilterDesign::IIR::iirfilter(2, fcV, 0, 0,
                      Utilities::FilterDesign::Bandtype::HIGHPASS,
                      Utilities::FilterDesign::IIRPrototype::BESSEL,
                      sos,
                      Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    ns = sos.getNumberOfSections();
    bs = sos.getNumeratorCoefficients();
    as = sos.getDenominatorCoefficients();
    std::vector<double> ytemp(npts);
    sosFilt.initialize(ns, bs.data(), as.data(),
                       RTSeis::ProcessingMode::POST_PROCESSING,
                       RTSeis::Precision::DOUBLE);
    sosFilt.apply(npts, x.data(), ytemp.data());
    std::reverse(ytemp.begin(), ytemp.end());
    sosFilt.apply(npts, ytemp.data(), ysosRef.data());
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
        waveform.getData(y);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ytemp.data(), y.data(), npts, &l1Norm);
    if (l1Norm > 1.e-12)
    {
        fprintf(stderr, "Failed sos filter test with error=%e\n", l1Norm);
        return EXIT_FAILURE;
    }
    // Bandpass filter
    fcV[0] = 1/fnyq;
    fcV[1] = 10/fnyq;
    Utilities::FilterDesign::IIR::iirfilter(2, fcV, 0, 0,
                      Utilities::FilterDesign::Bandtype::BANDPASS,
                      Utilities::FilterDesign::IIRPrototype::BUTTERWORTH,
                      sos,
                      Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    ns = sos.getNumberOfSections();
    bs = sos.getNumeratorCoefficients();
    as = sos.getDenominatorCoefficients();
    sosFilt.initialize(ns, bs.data(), as.data(),
                       RTSeis::ProcessingMode::POST_PROCESSING,
                       RTSeis::Precision::DOUBLE);
    sosFilt.apply(npts, x.data(), ytemp.data());
    std::reverse(ytemp.begin(), ytemp.end());
    sosFilt.apply(npts, ytemp.data(), ysosRef.data());
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
        waveform.getData(y);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ytemp.data(), y.data(), npts, &l1Norm);
    if (l1Norm > 1.e-12)
    {
        fprintf(stderr, "Failed sos bp filter test with error=%e\n", l1Norm);
        return EXIT_FAILURE;
    }
    // Bandpass filter
    fcV[0] = 1/fnyq;
    fcV[1] = 10/fnyq;
    Utilities::FilterDesign::IIR::iirfilter(3, fcV, 0, 0,
                          Utilities::FilterDesign::Bandtype::BANDSTOP,
                          Utilities::FilterDesign::IIRPrototype::BESSEL,
                          sos,
                          Utilities::FilterDesign::IIRFilterDomain::DIGITAL);
    ns = sos.getNumberOfSections();
    bs = sos.getNumeratorCoefficients();
    as = sos.getDenominatorCoefficients();
    sosFilt.initialize(ns, bs.data(), as.data(),
                       RTSeis::ProcessingMode::POST_PROCESSING,
                       RTSeis::Precision::DOUBLE);
    sosFilt.apply(npts, x.data(), ysosRef.data());
    sosFilt.clear();
    for (int j=0; j<1; j++)
    {
    try
    {
        // Design a 2nd order Bessel filter with 10 Hz cutoff
        int order = 3;                     // 2nd order (becomes 4 poles)
        std::pair<double,double> fc(1,10); // Passband is 1 to 10 Hz
        double ripple = 0;                 // N/A
        bool lzeroPhase = false;           // Not Zero-phase
        waveform.setData(x);
        waveform.sosBandstopFilter(order, fc,
                                   SingleChannel::IIRPrototype::BESSEL,
                                   ripple, lzeroPhase);
        waveform.getData(y);
    }
    catch (const std::invalid_argument &ia)
    {
        fprintf(stderr, "%s", ia.what());
        return EXIT_FAILURE;
    }
    }
    // Compare
    ippsNormDiff_L1_64f(ysosRef.data(), y.data(), npts, &l1Norm);
    if (l1Norm > 1.e-12)
    {
        fprintf(stderr, "Failed sos bs filter test with error=%e\n", l1Norm);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS; 
}

//============================================================================//

int testDemean(void)
{
    // Do two passes
    for (int j=0; j<2; j++)
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
        waveform.getData(y);    // Get the demeaned data in a vector, y
    }
    catch (std::invalid_argument &ia)
    {
        fprintf(stderr, "Demean failed %s\n", ia.what()); 
        return EXIT_FAILURE;
    }
//! [ppSCDemeanExample]
    // Verify
    if (y.size() != x.size())
    {
        RTSEIS_ERRMSG("%s", "Inconsistent sizes");
        return EXIT_FAILURE;
    }
    double maxAbs = 0;
    ippsMaxAbs_64f(y.data(), npts, &maxAbs);
    if (maxAbs > 1.e-13)
    {
        RTSEIS_ERRMSG("Demean failed %e", maxAbs);
    }
    } // Loop on passes

    return EXIT_SUCCESS;
}
//============================================================================//
int testDetrend(void)
{
    for (int j=0; j<2; j++)
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
        waveform.setData(x);    // Set the data to detrend
        waveform.detrend();     // Detrend
        waveform.getData(y);    // Get the detrended data in a vector, y
    }   
    catch (std::invalid_argument &ia)
    {
        fprintf(stderr, "Detrend failed %s\n", ia.what()); 
        return EXIT_FAILURE;
    }   
//! [ppSCDetrendExample]
    double maxAbs = 0;
    ippsMaxAbs_64f(y.data(), npts, &maxAbs);
    // ~50,000 points loses about 5 digits which is sensible 
    if (maxAbs > 1.e-9)
    {
        RTSEIS_ERRMSG("Demean failed %e", maxAbs);
        return EXIT_FAILURE;
    }
    }; // Loop on npts
    return EXIT_SUCCESS;
}
//============================================================================//
int testTaper(void)
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
            RTSEIS_ERRMSG("Premature end of file %s\n",
                          taperSolns100FileName.c_str());
            return EXIT_FAILURE;
        }
        std::sscanf(line100.c_str(), "%lf, %lf, %lf\n",
                    &yHamming100Ref[i], &yHanning100Ref[i], &ySine100Ref[i]);
        if (!std::getline(taper101File, line101))
        {
            RTSEIS_ERRMSG("Premature end of file %s\n",
                          taperSolns101FileName.c_str());
            return EXIT_FAILURE;
        }
        std::sscanf(line101.c_str(), "%lf, %lf, %lf\n",
                    &yHamming101Ref[i], &yHanning101Ref[i], &ySine101Ref[i]);
    }
    if (!std::getline(taper101File, line101))
    {
        RTSEIS_ERRMSG("Premature end of file %s\n",
                      taperSolns101FileName.c_str());
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
        double pct = 100*(0.2*2); // SAC to RTSeis; 20 pct of signal
        waveform.setData(x100);
        waveform.taper(pct, SC::TaperParameters::HAMMING);
        waveform.getData(yHamming100);

        pct = 100*(0.1*2); // SAC to RTSeis; 20 pct of signal
        waveform.setData(x100);
        waveform.taper(pct, SC::TaperParameters::HANN);
        waveform.getData(yHann100);

        pct = 100*(0.3*2); // SAC to RTSeis; 60 pct of signal
        waveform.setData(x100);
        waveform.taper(pct, SC::TaperParameters::SINE);
        waveform.getData(ySine100); 
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("Taper 100 failed %s", ia.what());
        return EXIT_FAILURE;
    }
    // Compare
    for (size_t i=0; i<ySine100.size(); i++)
    {
        if (std::abs(yHamming100[i] - yHamming100Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Sine 100 failed %ld,%lf,%lf", 
                          i, yHamming100[i], yHamming100Ref[i])
            return EXIT_FAILURE;
        }
        if (std::abs(yHann100[i] - yHanning100Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Sine 100 failed %ld,%lf,%lf", 
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
        waveform.getData(yHamming101);

        pct = 100*(0.1*2); // SAC to RTSeis; 20 pct of signal
        waveform.setData(x101);
        waveform.taper(pct, SC::TaperParameters::HANN);
        waveform.getData(yHann101);

        pct = 100*(0.15*2); // SAC to RTSeis; 30 pct of signal
        waveform.setData(x101);
        waveform.taper(pct, SC::TaperParameters::SINE);
        waveform.getData(ySine101); 
    }   
    catch (const std::invalid_argument &ia)
    {   
        RTSEIS_ERRMSG("Taper 101 failed %s", ia.what());
        return EXIT_FAILURE;
    }
    for (size_t i=0; i<101; i++)
    {
        if (std::abs(yHamming101[i] - yHamming101Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Sine 101 failed %ld,%lf,%lf", 
                          i, yHamming101[i], yHamming101Ref[i])
            return EXIT_FAILURE;
        }
        if (std::abs(yHann101[i] - yHanning101Ref[i]) > 1.e-6)
        {
            RTSEIS_ERRMSG("Sine 101 failed %ld,%lf,%lf", 
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
