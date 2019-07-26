#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <complex>
#include <vector>
#include <ipps.h>
#include "rtseis/utilities/filterDesign/fir.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterImplementations/decimate.hpp"
#include "rtseis/utilities/filterImplementations/downsample.hpp"
#include "rtseis/utilities/filterImplementations/iirFilter.hpp"
#include "rtseis/utilities/filterImplementations/iiriirFilter.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#include "rtseis/utilities/filterImplementations/multiRateFIRFilter.hpp"
#include "rtseis/utilities/filterImplementations/medianFilter.hpp"
#include "rtseis/utilities/filterImplementations/sosFilter.hpp"
#include "rtseis/utilities/filterImplementations/enums.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis::Utilities::FilterImplementations;
static void read_decimate(const int nq, std::vector<double> *xdecim);
static int readTextFile(int *npts, double *xPtr[],
                        const std::string fileName = "utils/data/gse2.txt");
static void lowpassFilterThenDownsample(const int npts, const int nfir,
                                        const int downFactor, const double x[], 
                                        std::vector<double> *y,
                                        const bool lremovePhase);
/*
static int filters_downsample_test(); //const int npts, const double x[]);
static int filters_firFilter_test(const int npts, const double x[],
                           const std::string fileName);
static int filters_firMRFilter_test(const int npts, const double x[],
                                    const std::string fileName);
static int filters_medianFilter_test(const int npts, const double x[],
                                     const std::string fileName);
static int filters_sosFilter_test(const int npts, const double x[],
                                  const std::string fileName);
int filters_iirFilter_test(const int npts, const double x[],
                           const std::string fileName1,
                           const std::string fileName2);
int filters_iiriirFilter_test(const int npts, const double x[],
                              const std::string fileName);

int rtseis_test_utils_filters(void)
{
//    double dt = 1.0/200.0;
    int npts;
    double *x = nullptr;
    std::string dataDir = "data/";
    // Load the data
    int ierr = readTextFile(&npts, &x, dataDir + "gse2.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed to read gse2 data");
        return EXIT_FAILURE;
    }
    // Apply the downsampler
    ierr = filters_downsample_test(); //npts, x);
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed downsampler test");
        return EXIT_FAILURE;
    }
    // Apply the IIR filter
    ierr = filters_iirFilter_test(npts, x,
                                  dataDir + "iirReference1.txt",
                                  dataDir + "iirReference2.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed iir filter test");
        return EXIT_FAILURE;
    } 
    // Apply the zero-phase IIR filter
    ierr = filters_iiriirFilter_test(npts, x,
                                     dataDir + "iiriirReference.txt"); 
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed iiriir filter test");
        return EXIT_FAILURE;
    }
    // Apply the median filter
    ierr = filters_medianFilter_test(npts, x,
                                     dataDir + "medianFilterReference.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed median filter test");
        return EXIT_FAILURE;
    }
    // Apply the second order section filter
    ierr = filters_sosFilter_test(npts, x,
                                  dataDir + "sosReference.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed sos filter test");
        return EXIT_FAILURE;
    }
    // Apply the FIR filter
    ierr = filters_firFilter_test(npts, x, dataDir + "firReference.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed fir filter test");
        return EXIT_FAILURE;
    }
    ierr = filters_firMRFilter_test(npts, x,
                                    dataDir + "firResampleReference.txt");
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed firmr filter test");
        return EXIT_FAILURE;
    }

    if (x != nullptr){free(x);}
    return EXIT_SUCCESS;
}
*/
//============================================================================//
//int filters_iirFilter_test(const int npts, const double x[],
//                           const std::string fileName1,
//                           const std::string fileName2)
TEST(UtilitiesFilterImplementations, iir)
{
    double *x = NULL;
    int npts;
    auto ierr = readTextFile(&npts, &x, "data/gse2.txt");
    EXPECT_EQ(ierr, 0);
    const std::string fileName1 = "data/iirReference1.txt";
    const std::string fileName2 = "data/iirReference2.txt";
    // Hardwire a high-order butterworth filter
    int na = 9;
    int nb = 9;
    const double b[9] = {0.000401587491686,
                         0.0,
                        -0.001606349966746,
                         0.0,
                         0.002409524950119,
                         0.0,
                        -0.001606349966746,
                         0.0,
                         0.000401587491686};
    const double a[9] = {1.000000000000000,
                        -7.185226122700763,
                        22.615376628798678,
                       -40.733465892344896,
                        45.926605646620146,
                       -33.196326377161412,
                        15.023103545324197,
                        -3.891997997268024,
                         0.441930568732716};
    int nb2 = 5;
    int na2 = 5;
    // [b,a] = signal.iirfilter(2, [0.5/200, 10/200], btype='band', ftype='butter')
    const double b2[5] = {0.005028154401050966,
                          0.0,
                         -0.010056308802101932,
                          0.0,
                          0.005028154401050966};
    const double a2[5] = {1.0,
                         -3.787291923503399,
                          5.384566634300907,
                         -3.407019620837378,
                          0.8097462844356603};

    // Load a reference solution
    double *yref1 = nullptr;
    double *yref2 = nullptr;
    int npref = 0; 
    ierr = readTextFile(&npref, &yref1, fileName1);
    EXPECT_EQ(ierr, 0);
    EXPECT_EQ(npts, npref);
    ierr = readTextFile(&npref, &yref2, fileName2);
    EXPECT_EQ(ierr, 0);
    EXPECT_EQ(npts, npref);
    // Compute the zero-phase IIR filter alternative
    IIRFilter iir;
    IIRFilter iir_slow;
    EXPECT_NO_THROW(iir.initialize(nb, b, na, a,
                                   RTSeis::ProcessingMode::POST_PROCESSING,
                                   RTSeis::Precision::DOUBLE,
                                   IIRDFImplementation::DF2_FAST));
    EXPECT_NO_THROW(iir_slow.initialize(nb, b, na, a,
                                        RTSeis::ProcessingMode::POST_PROCESSING,
                                        RTSeis::Precision::DOUBLE,
                                        IIRDFImplementation::DF2_SLOW));
    double *yref = new double[npts];
    double *yref_slow = new double[npts];
    double *y1 = new double[npts];
    double *y2 = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(iir.apply(npts, x, &y1));
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double error = 0;
    ippsNormDiff_L2_64f(y1, yref1, npts, &error);
    error = error/static_cast<double> (npts);
    EXPECT_LE(error, 3.e-2);
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Fast reference solution 1 computation time %.8lf (s)\n",
            tdif.count());
    std::copy(y1, y1+npts, yref);

    // Repeat for slow
    timeStart = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(iir_slow.apply(npts, x, &y2));
    timeEnd = std::chrono::high_resolution_clock::now();
    ippsNormDiff_Inf_64f(y2, yref1, npts, &error);
    EXPECT_LE(error, 1.e-3);
    tdif = timeEnd - timeStart;
    fprintf(stdout, "Slow reference solution 1 computation time %.8lf (s)\n",
            tdif.count());
    std::copy(y2, y2+npts, yref_slow);
    // Try a lower order filter - Note the difference in accuracy.
    // I think the filter is a little unstable.
    IIRFilter iir2;
    EXPECT_NO_THROW(iir2.initialize(nb2, b2, na2, a2,
                                    RTSeis::ProcessingMode::POST_PROCESSING,
                                    RTSeis::Precision::DOUBLE,
                                    IIRDFImplementation::DF2_FAST));
    timeStart = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(iir2.apply(npts, x, &y2));
    timeEnd = std::chrono::high_resolution_clock::now();
    error = 0;
    ippsNormDiff_Inf_64f(y2, yref2, npts, &error);
    EXPECT_LE(error, 5.e-7);
    tdif = timeEnd - timeStart;
    fprintf(stdout, "Fast reference solution 2 computation time %.8lf (s)\n",
            tdif.count());
    // Do a real-time test of high-order filter
    IIRFilter iirrt;
    EXPECT_NO_THROW(iirrt.initialize(nb, b, na, a,
                                     RTSeis::ProcessingMode::REAL_TIME,
                                     RTSeis::Precision::DOUBLE,
                                     IIRDFImplementation::DF2_FAST));
    iir = iirrt;
    std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                 1000, 1024, 1200, 2048, 4000, 4096, 5000});
    for (int job=0; job<2; job++)
    {
        for (size_t ip=0; ip<packetSize.size(); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0;
            int nptsPass = 0;
            while (nxloc < npts)
            {
                nptsPass = packetSize[ip];
                if (job == 1)
                {
                     nptsPass = std::max(1, nptsPass + rand()%50 - 25);
                }
                nptsPass = std::min(nptsPass, npts - nxloc);
                double *y1Temp = &y1[nxloc];
                EXPECT_NO_THROW(iir.apply(nptsPass, &x[nxloc], &y1Temp));//[nxloc]);
                nxloc = nxloc + nptsPass;
            }
            iir.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            ippsNormDiff_Inf_64f(yref, y1, npts, &error);
            EXPECT_LE(error, 1.e-12);
            if (job == 0)
            {
                fprintf(stdout,
                        "Passed IIR filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed IIR filter random in %.8e (s)\n", tdif.count());
            }
        }
    }
    // Retry this for the slow filter implementation
    EXPECT_NO_THROW(iir_slow.initialize(nb, b, na, a,
                                        RTSeis::ProcessingMode::REAL_TIME,
                                        RTSeis::Precision::DOUBLE,
                                        IIRDFImplementation::DF2_SLOW));
    iir = iir_slow;
    for (int job=0; job<2; job++)
    {    
        for (size_t ip=0; ip<packetSize.size(); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0; 
            int nptsPass = 0; 
            while (nxloc < npts)
            {
                nptsPass = packetSize[ip];
                if (job == 1)
                {
                     nptsPass = std::max(1, nptsPass + rand()%50 - 25); 
                }
                nptsPass = std::min(nptsPass, npts - nxloc);
                double *y1Temp = &y1[nxloc];
                EXPECT_NO_THROW(iir.apply(nptsPass, &x[nxloc], &y1Temp)); //y1[nxloc]);
                nxloc = nxloc + nptsPass;
            }
            iir.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            ippsNormDiff_Inf_64f(yref_slow, y1, npts, &error);
            EXPECT_LE(error, 1.e-12);
            if (job == 0)
            {
                fprintf(stdout,
                        "Passed IIR slow filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed IIR slow filter random in %.8e (s)\n", tdif.count());
            }
        }
    }

    free(yref1);
    free(yref2);
    delete[] yref;
    delete[] yref_slow;
    delete[] y1;
    delete[] y2;
    free(x);
}
//============================================================================//
//int filters_iiriirFilter_test(const int npts, const double x[],
//                              const std::string fileName)
TEST(UtilitiesFilterImplementations, iiriir)
{
    double *x = NULL;
    int npts;
    auto ierr = readTextFile(&npts, &x, "data/gse2.txt");
    EXPECT_EQ(ierr, 0);
    const std::string fileName = "data/iiriirReference.txt";
    // Hardwire a bandpass butterworth filter
    const int na = 9;
    const int nb = 9;
    const double b[9] = {0.000401587491686,
                         0.0,
                        -0.001606349966746,
                         0.0,
                         0.002409524950119,
                         0.0,
                        -0.001606349966746,
                         0.0,
                         0.000401587491686};
    const double a[9] = {1.000000000000000,
                        -7.185226122700763,
                        22.615376628798678,
                       -40.733465892344896,
                        45.926605646620146,
                       -33.196326377161412,
                        15.023103545324197,
                        -3.891997997268024,
                         0.441930568732716};
    // Load a reference solution
    double *yref = nullptr;
    int npref = 0;
    ierr = readTextFile(&npref, &yref, fileName);
    EXPECT_EQ(ierr, 0);
    EXPECT_EQ(npts, npref);
    // Compute the zero-phase IIR filter alternative
    IIRIIRFilter iiriir;
    EXPECT_NO_THROW(iiriir.initialize(nb, b, na, a,
                                      RTSeis::Precision::DOUBLE));
    double *y = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(iiriir.apply(npts, x, &y));
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double error = 0;
    ippsNormDiff_L2_64f(y, yref, npts, &error);
    error = error/static_cast<double> (npts);
    EXPECT_LE(error, 4.e-2);
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    delete[] y;
    free(yref);

    // Do the Intel test
    iiriir.clear();
    const int LEN = 256;
    const double SMP_RATE = 1000;
    const double F1 = 50;
    const double F2 = 120;
    const double AMPL = 10;
    const double CUT_OFF = 0.3;
    const int ORDER = 3;
    Ipp64f *f1 = ippsMalloc_64f(LEN);
    Ipp64f *f2 = ippsMalloc_64f(LEN);
    double *xi = new double[LEN];
    double *yi = new double[LEN];
    double *yiRef = new double[LEN];
    // Generate 50 Hz sine tone
    Ipp64f phase = 3*static_cast<Ipp64f> (IPP_PI2);
    ippsTone_64f(f1, LEN, AMPL, F1/SMP_RATE, &phase, ippAlgHintAccurate);
    phase =  3*static_cast<Ipp64f> (IPP_PI2); // Restore phase
    // Generate 120 Hz sine tone
    ippsTone_64f(f2, LEN, AMPL, F2/SMP_RATE, &phase, ippAlgHintAccurate);
    ippsAdd_64f(f1, f2, xi, LEN); 
    ippsAbs_64f_I(xi, LEN); // Make positive; 
    // Make a filter
    int bufferSizeGen, bufferSizeIIRIIR;
    ippsIIRIIRGetStateSize_64f(ORDER, &bufferSizeIIRIIR);
    ippsIIRGenGetBufferSize(ORDER, &bufferSizeGen);
    Ipp8u *pBufIIRIIR = ippsMalloc_8u(IPP_MAX(bufferSizeGen, bufferSizeIIRIIR));
    Ipp64f *pTaps = ippsMalloc_64f( 2 * ( ORDER + 1 ));
    ippsIIRGenLowpass_64f(CUT_OFF/2, 0.0, ORDER, pTaps,
                          ippButterworth, pBufIIRIIR);
    // Filter
    IppsIIRState_64f *pStateIIRIIR = NULL;
    ippsIIRIIRInit_64f(&pStateIIRIIR, pTaps, ORDER, NULL, pBufIIRIIR );
    ippsIIRIIR_64f(xi, yiRef, LEN, pStateIIRIIR);
    // Repeat with RTSeis
    EXPECT_NO_THROW(iiriir.initialize(ORDER+1, &pTaps[0],
                                      ORDER+1, &pTaps[ORDER+1]));
    EXPECT_NO_THROW(iiriir.apply(LEN, xi, &yi));
    ippsNormDiff_Inf_64f(yi, yiRef, LEN, &error);
    EXPECT_LE(error, 1.e-14);
    ippsFree(f1);
    ippsFree(f2);
    ippsFree(pTaps);
    ippsFree(pBufIIRIIR);
    delete[] xi;
    delete[] yi;
    delete[] yiRef;
    free(x);
}
//============================================================================//
//int filters_firMRFilter_test(const int npts, const double x[],
//                             const std::string fileName)
TEST(UtiltiesFilterImplementations, multirateFIR)
{
    double *x = NULL;
    int npts;
    auto ierr = readTextFile(&npts, &x, "data/gse2.txt");
    EXPECT_EQ(ierr, 0);
    const std::string fileName = "data/firResampleReference.txt";
    fprintf(stdout, "Testing FIR multirate filter...\n");
    const int nb = 51; 
    int upFactor = 3;
    int downFactor = 8; 
    //const int na = 1;
    const double b[51] = {-0.000000000000000, -0.001056235801065,
                          -0.000769341020163,  0.000956323223723,
                           0.001976082742122, -0.000000000000000,
                          -0.003265384800345, -0.002568519852901,
                           0.003234633130890,  0.006519908075213,
                          -0.000000000000000, -0.009825739114867,
                          -0.007365685405410,  0.008881348924986,
                           0.017256056989442, -0.000000000000000,
                          -0.024784271698734, -0.018417666768131,
                           0.022299534288278,  0.044222443880910,
                          -0.000000000000000, -0.071469809226860,
                          -0.060430328816090,  0.092317626953209,
                           0.302027315266443,  0.400523418058701,
                           0.302027315266443,  0.092317626953209,
                          -0.060430328816090, -0.071469809226860,
                          -0.000000000000000,  0.044222443880910,
                           0.022299534288278, -0.018417666768131,
                          -0.024784271698734, -0.000000000000000,
                           0.017256056989442,  0.008881348924986,
                          -0.007365685405410, -0.009825739114867,
                          -0.000000000000000,  0.006519908075213,
                           0.003234633130890, -0.002568519852901,
                          -0.003265384800345, -0.000000000000000,
                           0.001976082742122,  0.000956323223723,
                          -0.000769341020163, -0.001056235801065,
                          -0.000000000000000};
    // Load a reference solution
    double *yref = nullptr;
    int npref = 0;
    ierr = readTextFile(&npref, &yref, fileName);
    EXPECT_EQ(ierr, 0);
    // Initialize filter
    //double gain = static_cast<double> (upFactor);
    //ippsMulC_64f(b, gain, bgain, nb); // Fix gain
    MultiRateFIRFilter firmr;
    EXPECT_NO_THROW(firmr.initialize(upFactor, downFactor, nb, b,
                                     RTSeis::ProcessingMode::POST_PROCESSING,
                                     RTSeis::Precision::DOUBLE));
    // Estimate space in output
    int nest = firmr.estimateSpace(npts);
    EXPECT_EQ(nest, 4500);
    //if (nest != 4500)
    //{
    //    RTSEIS_WARNMSG("%s", "Anticipating 4500 points in comparison");
    //}
    int ncomp = std::min(nest, npref);
    // Filter 
    double *y = new double[npts];
    int ny;
    auto timeStart = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(firmr.apply(npts, x, npts, &ny, &y));
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double error = 0;
    ippsNormDiff_Inf_64f(y, yref, ncomp, &error); 
    EXPECT_LE(error, 1.e-10); 
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Repeat for real-time 
    MultiRateFIRFilter firmrRT;
    EXPECT_NO_THROW(firmrRT.initialize(upFactor, downFactor, nb, b,
                                       RTSeis::ProcessingMode::REAL_TIME,
                                       RTSeis::Precision::DOUBLE));
    firmr = firmrRT; 
    std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                 1000, 1024, 1200, 2048, 4000, 4096, 5000});
    for (auto job=0; job<2; job++)
    {   
        for (auto ip=0; ip<static_cast<int> (packetSize.size()); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0;
            int nyloc = 0;
            int nptsPass = 0;
            while (nxloc < npts)
            {
                nptsPass = packetSize[ip];
                if (job == 1)
                {
                     nptsPass = std::max(1, nptsPass + rand()%50 - 25);
                }
                nptsPass = std::min(nptsPass, npts - nxloc);
                int nwork = npts+1-nyloc;
                int nyDec = 0;
                double *yptr = &y[nyloc];
                EXPECT_NO_THROW(firmr.apply(nptsPass, &x[nxloc], nwork,
                                            &nyDec, &yptr)); //&y[nyloc]);
                nxloc = nxloc + nptsPass;
                nyloc = nyloc + nyDec;
            }
            firmr.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            EXPECT_EQ(nyloc, ny);
            ippsNormDiff_Inf_64f(yref, y, nyloc, &error);
            EXPECT_LE(error, 1.e-8);
            if (job == 0)
            {    
                fprintf(stdout,
                        "Passed FIR multirate filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else 
            {
                fprintf(stdout,
                        "Passed FIR multirate filter random in %.8e (s)\n", tdif.count());
            }
        }
    }
    delete[] y;
    free(yref);
    free(x);
}
//============================================================================//
//int filters_firFilter_test(const int npts, const double x[],
//                           const std::string fileName)
TEST(UtilitiesFilterImplementations, fir)
{
    double *x = NULL;
    int npts;
    auto ierr = readTextFile(&npts, &x, "data/gse2.txt");
    EXPECT_EQ(ierr, 0);
    const std::string fileName = "data/firReference.txt";
    const int nb = 51;
    //const int na = 1;
    const double b[51] = {-0.000000000000000, -0.001056235801065,
                          -0.000769341020163,  0.000956323223723,
                           0.001976082742122, -0.000000000000000,
                          -0.003265384800345, -0.002568519852901,
                           0.003234633130890,  0.006519908075213,
                          -0.000000000000000, -0.009825739114867,
                          -0.007365685405410,  0.008881348924986,
                           0.017256056989442, -0.000000000000000,
                          -0.024784271698734, -0.018417666768131,
                           0.022299534288278,  0.044222443880910,
                          -0.000000000000000, -0.071469809226860,
                          -0.060430328816090,  0.092317626953209,
                           0.302027315266443,  0.400523418058701,
                           0.302027315266443,  0.092317626953209,
                          -0.060430328816090, -0.071469809226860,
                          -0.000000000000000,  0.044222443880910,
                           0.022299534288278, -0.018417666768131,
                          -0.024784271698734, -0.000000000000000,
                           0.017256056989442,  0.008881348924986,
                          -0.007365685405410, -0.009825739114867,
                          -0.000000000000000,  0.006519908075213,
                           0.003234633130890, -0.002568519852901,
                          -0.003265384800345, -0.000000000000000,
                           0.001976082742122,  0.000956323223723,
                          -0.000769341020163, -0.001056235801065,
                          -0.000000000000000};
    // Load a reference solution
    double *yref = nullptr;
    int npref = 0;
    ierr = readTextFile(&npref, &yref, fileName);
    EXPECT_EQ(ierr, 0);
    EXPECT_EQ(npts, npref);
    // Make a post-processing solution
    FIRFilter fir;
    EXPECT_NO_THROW(fir.initialize(nb, b,
                                   RTSeis::ProcessingMode::POST_PROCESSING,
                                   RTSeis::Precision::DOUBLE,
                                   FIRImplementation::DIRECT));
    double *y = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(fir.apply(npts, x, &y));
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double error = 0;
    ippsNormDiff_Inf_64f(y, yref, npts, &error);
    EXPECT_LE(error, 1.e-10);
    std::copy(y, y+static_cast<size_t> (npts), yref);
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Do packetized tests
    FIRFilter firrt;
    EXPECT_NO_THROW(firrt.initialize(nb, b,
                                     RTSeis::ProcessingMode::REAL_TIME,
                                     RTSeis::Precision::DOUBLE,
                                     FIRImplementation::DIRECT));
    fir = firrt;
    std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                 1000, 1024, 1200, 2048, 4000, 4096, 5000});
    for (int job=0; job<2; job++)
    {
        for (size_t ip=0; ip<packetSize.size(); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0;
            int nptsPass = 0;
            while (nxloc < npts)
            {
                nptsPass = packetSize[ip];
                if (job == 1)
                {
                     nptsPass = std::max(1, nptsPass + rand()%50 - 25);
                }
                nptsPass = std::min(nptsPass, npts - nxloc);
                double *yptr = &y[nxloc];
                EXPECT_NO_THROW(fir.apply(nptsPass, &x[nxloc], &yptr));
                nxloc = nxloc + nptsPass;
            }
            fir.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            ippsNormDiff_Inf_64f(yref, y, npts, &error);
            EXPECT_LE(error, 1.e-10);
            if (job == 0)
            {
                fprintf(stdout,
                        "Passed FIR filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed FIR filter random in %.8e (s)\n", tdif.count());
            }
        }
    }
    delete[] y;
    free(yref);
    free(x);
}
//============================================================================//
//int filters_sosFilter_test(const int npts, const double x[],
//                           const std::string fileName)
TEST(UtilitiesFilterImplementations, sos)
{
    double *x = NULL;
    int npts;
    auto ierr = readTextFile(&npts, &x, "data/gse2.txt");
    EXPECT_EQ(ierr, 0);
    const std::string fileName = "data/sosReference.txt";
    int ns = 7; // number of sections 
    const double bs7[21] = {6.37835424e-05,  6.37835424e-05,  0.00000000e+00, 
                           1.00000000e+00, -1.78848938e+00,  1.00000000e+00,
                           1.00000000e+00, -1.93118487e+00,  1.00000000e+00,
                           1.00000000e+00, -1.95799864e+00,  1.00000000e+00,
                           1.00000000e+00, -1.96671846e+00,  1.00000000e+00,
                           1.00000000e+00, -1.97011885e+00,  1.00000000e+00,
                           1.00000000e+00, -1.97135784e+00,  1.00000000e+00};
    const double as7[21] = {1.00000000e+00, -9.27054679e-01,  0.00000000e+00,
                           1.00000000e+00, -1.87008942e+00,  8.78235919e-01,
                           1.00000000e+00, -1.90342568e+00,  9.17455718e-01,
                           1.00000000e+00, -1.93318668e+00,  9.52433552e-01,
                           1.00000000e+00, -1.95271141e+00,  9.75295685e-01,
                           1.00000000e+00, -1.96423610e+00,  9.88608056e-01,
                           1.00000000e+00, -1.97157693e+00,  9.96727086e-01};
    double yref40[40] = {6.37835424e-05,  1.23511272e-04,  1.34263690e-04,
                         1.78634911e-04,  2.50312740e-04,  3.46332848e-04,
                         4.66239952e-04,  6.11416691e-04,  7.84553129e-04,
                         9.89232232e-04,  1.22960924e-03,  1.51016546e-03,
                         1.83551947e-03,  2.21028135e-03,  2.63893773e-03,
                         3.12575784e-03,  3.67471270e-03,  4.28940130e-03,
                         4.97297977e-03,  5.72809028e-03,  6.55678845e-03,
                         7.46046851e-03,  8.43978671e-03,  9.49458408e-03,
                         1.06238101e-02,  1.18254496e-02,  1.30964547e-02,
                         1.44326848e-02,  1.58288573e-02,  1.72785101e-02,
                         1.87739799e-02,  2.03063976e-02,  2.18657022e-02,
                         2.34406756e-02,  2.50189979e-02,  2.65873261e-02,
                         2.81313940e-02,  2.96361349e-02,  3.10858256e-02,
                         3.24642512e-02};
    double *impulse = new double[40];
    double *y40 = new double[40]; 
    std::fill(impulse, impulse+40, 0);
    impulse[0] = 1;
    SOSFilter sos;
    EXPECT_NO_THROW(sos.initialize(ns, bs7, as7,
                                   RTSeis::ProcessingMode::POST_PROCESSING,
                                   RTSeis::Precision::DOUBLE));
    EXPECT_NO_THROW(sos.apply(40, impulse, &y40));
    double error;
    ippsNormDiff_Inf_64f(y40, yref40, 40, &error);
    EXPECT_LE(error, 1.e-8);
    delete[] y40;
    delete[] impulse;
    // Load a reference solution
    double *yref = nullptr;
    int npref = 0;
    ierr = readTextFile(&npref, &yref, fileName);
    EXPECT_EQ(ierr, 0);
    EXPECT_EQ(npts, npref);
    ns = 4;
    const double bs[12] = {0.000401587491686,  0.000803175141692,  0.000401587491549,
                           1.000000000000000, -2.000000394412897,  0.999999999730209,
                           1.000000000000000,  1.999999605765104,  1.000000000341065,
                           1.000000000000000, -1.999999605588274,  1.000000000269794};
    const double as[12] = {1.000000000000000, -1.488513049541281,  0.562472929601870,
                           1.000000000000000, -1.704970593447777,  0.792206889942566,
                           1.000000000000000, -1.994269533089365,  0.994278822534674,
                           1.000000000000000, -1.997472946622339,  0.997483252685326};
    EXPECT_NO_THROW(sos.initialize(ns, bs, as,
                                   RTSeis::ProcessingMode::POST_PROCESSING,
                                   RTSeis::Precision::DOUBLE));
    double *y = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    EXPECT_NO_THROW(sos.apply(npts, x, &y));
    auto timeEnd = std::chrono::high_resolution_clock::now();
    ippsNormDiff_Inf_64f(y, yref, npts, &error);
    EXPECT_LE(error, 1.e-8);
    std::copy(y, y+static_cast<size_t> (npts), yref);
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Do packetized tests 
    SOSFilter sosrt;
    EXPECT_NO_THROW(sosrt.initialize(ns, bs, as,
                                     RTSeis::ProcessingMode::REAL_TIME,
                                     RTSeis::Precision::DOUBLE));
    sos = sosrt;
    std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                 1000, 1024, 1200, 2048, 4000, 4096, 5000});
    for (auto job=0; job<2; job++)
    {
        for (auto ip=0; ip<static_cast<int> (packetSize.size()); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0;
            int nptsPass = 0;
            while (nxloc < npts)
            {
                nptsPass = packetSize[ip];
                if (job == 1)
                {
                     nptsPass = std::max(1, nptsPass + rand()%50 - 25);
                }
                nptsPass = std::min(nptsPass, npts - nxloc);
                double *yptr = &y[nxloc];
                EXPECT_NO_THROW(sos.apply(nptsPass, &x[nxloc], &yptr)); //&y[nxloc]);
                nxloc = nxloc + nptsPass;
            }
            sos.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            ippsNormDiff_Inf_64f(yref, y, npts, &error);
            EXPECT_LE(error, 1.e-8);
            if (job == 0)
            {
                fprintf(stdout,
                        "Passed SOS filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed SOS filter random in %.8e (s)\n", tdif.count());
            }
        }
    }
    delete[] y;
    free(yref);
    free(x);
}
//============================================================================//
//int filters_medianFilter_test(const int npts, const double x[],
//                              const std::string fileName)
TEST(UtilitiesFilterImplementations, medianFilter)
{
    double *x = NULL;
    int npts;
    auto ierr = readTextFile(&npts, &x, "data/gse2.txt");
    EXPECT_EQ(ierr, 0);
    const std::string fileName = "data/medianFilterReference.txt";
    double xin[8] = {1, 2, 127, 4, 5, 0, 7, 8};
    double y8[8];
    double yref3[8] = {1, 2, 4, 5, 4, 5, 7, 7}; // Matlab soln; IPP has edge effect
    double yref5[8] = {1, 2, 4, 4, 5, 5, 5, 0};
    MedianFilter median;
    EXPECT_NO_THROW(
        median.initialize(3,
                          RTSeis::ProcessingMode::POST_PROCESSING,
                          RTSeis::Precision::DOUBLE));
    double *yptr = &y8[0];
    EXPECT_NO_THROW(median.apply(8, xin, &yptr));
    double error;
    ippsNormDiff_Inf_64f(&y8[2], &yref3[1], 6, &error);
    EXPECT_LE(error, 1.e-14);

    EXPECT_NO_THROW(median.initialize(5,
                    RTSeis::ProcessingMode::POST_PROCESSING,
                    RTSeis::Precision::DOUBLE));
    yptr = &y8[0];
    median.apply(8, xin, &yptr);
    ippsNormDiff_Inf_64f(&y8[4], &yref5[2], 4, &error);
    EXPECT_LE(error, 1.e-14);
    // Load a reference solution
    double *yref = nullptr;
    int npref;
    ierr = readTextFile(&npref, &yref, fileName);
    EXPECT_EQ(ierr, 0);
    EXPECT_EQ(npref - 11/2, npts);
    EXPECT_NO_THROW(median.initialize(11,
                    RTSeis::ProcessingMode::POST_PROCESSING,
                    RTSeis::Precision::DOUBLE));
    auto timeStart = std::chrono::high_resolution_clock::now();
    double *y = new double[npts];
    EXPECT_NO_THROW(median.apply(npts, x, &y));
    auto timeEnd = std::chrono::high_resolution_clock::now();
    ippsNormDiff_Inf_64f(y, yref, npts, &error);
    EXPECT_LE(error, 1.e-14);
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Now do the packetized tests
    MedianFilter medianrt = median;
    EXPECT_NO_THROW(medianrt.initialize(11,
                                        RTSeis::ProcessingMode::REAL_TIME,
                                        RTSeis::Precision::DOUBLE));
    std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                 1000, 1024, 1200, 2048, 4000, 4096, 5000});
    for (auto job=0; job<2; job++)
    {
        for (auto ip=0; ip<static_cast<int> (packetSize.size()); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0;
            int nptsPass = 0;
            while (nxloc < npts)
            {
                nptsPass = packetSize[ip];
                if (job == 1)
                {
                     nptsPass = std::max(1, nptsPass + rand()%50 - 25);  
                }
                nptsPass = std::min(nptsPass, npts - nxloc);
                double *yptr = &y[nxloc];
                EXPECT_NO_THROW(medianrt.apply(nptsPass, &x[nxloc], &yptr)); //[nxloc]);
                nxloc = nxloc + nptsPass;
            }
            medianrt.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            ippsNormDiff_Inf_64f(yref, y, npts, &error);
            EXPECT_LE(error, 1.e-10);
            if (job == 0)
            {
                fprintf(stdout,
                        "Passed median filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed median filter random in %.8e (s)\n", tdif.count());
            }
        }
    }   
    free(yref);
    delete[] y;
    free(x);
}
//============================================================================//
//int filters_downsample_test() //const int npts, const double x[])
TEST(UtilitiesFilterImplementations, downsample)
{
    double *x = NULL;
    int npts;
    auto ierr = readTextFile(&npts, &x, "data/gse2.txt");
    EXPECT_EQ(ierr, 0);
    const int nq = 7;
    RTSeis::Precision precision = RTSeis::Precision::DOUBLE;
    // Call this in post-processing for a couple different decimation rates
    srand(10245);
    Downsample downsample;
    double *y = static_cast<double *>
                (calloc(static_cast<size_t> (npts), sizeof(double)));
    double *yref = static_cast<double *>
                   (calloc(static_cast<size_t> (npts), sizeof(double)));
    for (auto iq=1; iq<nq+1; iq++)
    {
        // Do a post-processing test
        memset(y, 0, static_cast<size_t> (npts)*sizeof(double));
        memset(yref, 0, static_cast<size_t> (npts)*sizeof(double));
        EXPECT_NO_THROW(
            downsample.initialize(iq,
                                  RTSeis::ProcessingMode::POST_PROCESSING,
                                  precision)
        );
        auto timeStart = std::chrono::high_resolution_clock::now();
        int ny;
        EXPECT_NO_THROW(downsample.apply(npts, x, npts, &ny, &y));
        auto timeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> tdif = timeEnd - timeStart;
        // Do a manual downsample
        int j = 0;
        double error = 0;
        for (int i=0; i<npts; i=i+iq)
        {
            error = std::max(error, y[j] - x[i]);
            j = j + 1;
        }
        EXPECT_EQ(j, ny);
        EXPECT_LE(error, 1.e-10);
        downsample.clear();
        fprintf(stdout,
                "Post-processing execution time for nq=%d is %.8lf (s)\n",
                iq, tdif.count());
        // Make a copy of the correct answer
        int nyref = ny;
        std::copy(y, y+ny, yref);
        //for (auto iy=0; iy<ny; iy++){yref[iy] = y[iy];}
        // Do a real-time test
        EXPECT_NO_THROW(
            downsample.initialize(iq,
                                  RTSeis::ProcessingMode::REAL_TIME,
                                  precision));
        std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                     1000, 1024, 1200, 2048, 4000, 4096, 5000});
        for (auto ip=0; ip<static_cast<int> (packetSize.size()); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0;
            int nyloc = 0;
            while (nxloc < npts)
            {
                int nptsPass = std::min(packetSize[ip], npts - nxloc);
                int nyDec = 0;
                double *yptr = &y[nyloc];
                EXPECT_NO_THROW(
                    downsample.apply(nptsPass, &x[nxloc],
                                     npts+1-nyloc, &nyDec, &yptr)
                );
                nxloc = nxloc + nptsPass;
                nyloc = nyloc + nyDec;
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            downsample.resetInitialConditions();
            EXPECT_EQ(nyloc, nyref);
            ippsNormDiff_Inf_64f(y, yref, nyref, &error);
            EXPECT_LE(error, 1.e-10);
            if (iq == 7)
            {
                fprintf(stdout,
                   "Passed downsampler fixed packet size %4d w/ nq=%d in %.8e (s)\n",
                    packetSize[ip], iq, tdif.count());
            }
        }
        // Random packet sizes
        timeStart = std::chrono::high_resolution_clock::now();
        int nxloc = 0;
        int nyloc = 0;
        int packetLen = 100;
        while (nxloc < npts)
        {
            int nptsPass = std::min(packetLen, npts - nxloc);
            int nyDec = 0;
            double *yptr = &y[nyloc];
            EXPECT_NO_THROW(downsample.apply(nptsPass, &x[nxloc],
                                             npts+1-nyloc, &nyDec, &yptr));
            nxloc = nxloc + nptsPass;
            nyloc = nyloc + nyDec;
            packetLen = std::max(1, packetLen + rand()%50 - 25);
        }
        timeEnd = std::chrono::high_resolution_clock::now();
        EXPECT_EQ(nyloc, nyref);
        ippsNormDiff_Inf_64f(y, yref, nyref, &error);
        EXPECT_LE(error, 1.e-10);
        tdif = timeEnd - timeStart;
        if (iq == 7)
        {
            fprintf(stdout,
                    "Passed downsampler random packet size w/ nq=%d in %.8e (s)\n",
                    iq, tdif.count());
        }
        // Loop 
        downsample.clear();
    }
    free(y);
    free(yref);
    free(x);
}
//============================================================================//
TEST(UtilitiesFilterImplementations, decimate)
{
    double *x = NULL;
    int npts;
    auto ierr = readTextFile(&npts, &x, "data/gse2.txt");
    EXPECT_EQ(ierr, 0);
    int nywork = npts;
    std::vector<double> y;
    std::vector<double> xpad;
    std::vector<double> ylpds;
    // Decimate by some default factors
    Decimate decimate;
    //int downFactor = 5;
    int filterLength = 95;
    bool lremovePhaseShift = true;
    for (auto downFactor=2; downFactor<=10; downFactor++)
    {
        // Read the reference solution
        std::vector<double> xdecim_ref;
        read_decimate(downFactor, &xdecim_ref);
        // Now do the decimation
        EXPECT_NO_THROW(decimate.initialize(downFactor,
                                            filterLength,
                                            lremovePhaseShift,
                                            RTSeis::ProcessingMode::POST_PROCESSING,
                                            RTSeis::Precision::DOUBLE));
        EXPECT_TRUE(decimate.isInitialized());
        int ndecim = decimate.estimateSpace(npts);
        int nyDecim = 0;
        y.resize(nywork);
        std::fill(y.begin(), y.end(), 0);
        auto timeStart = std::chrono::high_resolution_clock::now();
        auto yptr = y.data();
        EXPECT_NO_THROW(decimate.apply(npts, x, ndecim, &nyDecim, &yptr));
        EXPECT_EQ(nyDecim, ndecim);
        EXPECT_EQ(static_cast<int> (xdecim_ref.size()), nyDecim);
        auto timeEnd = std::chrono::high_resolution_clock::now();
        // Lowpass filter and downsample a reference signal
        int nfir = decimate.getFIRFilterLength();
        EXPECT_EQ(nfir%2, 1);
        double error = 0;
        ippsNormDiff_Inf_64f(y.data(), xdecim_ref.data(), nyDecim, &error);
        EXPECT_LE(error, 1.e-8);
        std::chrono::duration<double> tdif = timeEnd - timeStart;
        fprintf(stdout, "Reference solution time for nq=%d in %.8e (s)\n",
                downFactor, tdif.count());
        // Now do a real time solution
        Decimate rtDecim;
        EXPECT_NO_THROW(rtDecim.initialize(downFactor,
                                           nfir, //filterLength,
                                           false,
                                           RTSeis::ProcessingMode::REAL_TIME,
                                           RTSeis::Precision::DOUBLE));
        decimate = rtDecim; // Test copy assignent
        std::fill(y.begin(), y.end(), 0.0);
        std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512, 
                                     1000, 1024, 1200, 2048, 4000, 4096, 5000});
        for (auto job=0; job<2; job++)
        {
            for (auto ip=0; ip<static_cast<int> (packetSize.size()); ip++)
            {
                timeStart = std::chrono::high_resolution_clock::now();
                int nxloc = 0;
                int nyloc = 0;
                while (nxloc < npts)
                {
                    int nptsPass = std::min(packetSize[ip], npts - nxloc);
                    if (job == 1)
                    {
                        nptsPass = std::min(packetSize[ip] + rand()%50 - 25,
                                            npts - nxloc);
                    }
                    int nyDec = 0;
                    double *yptr = &y[nyloc];
                    EXPECT_NO_THROW(decimate.apply(nptsPass, &x[nxloc],
                                                   npts+1-nyloc, &nyDec, &yptr));
                    if (nptsPass <= 0){continue;}
                    nxloc = nxloc + nptsPass;
                    nyloc = nyloc + nyDec;
                }
                auto timeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> tdif = timeEnd - timeStart;
                decimate.resetInitialConditions();
                EXPECT_EQ(nyloc, nyDecim);
                int groupDelay = nfir/2;
                int ncomp = std::min(nyloc, nyDecim) - groupDelay/downFactor;
                ippsNormDiff_Inf_64f(y.data()+groupDelay/downFactor,
                                     xdecim_ref.data(), ncomp, &error);
                EXPECT_LE(error, 1.e-8);
                if (job == 0 && downFactor == 6)
                {
                    fprintf(stdout,
                            "Passed decimate fixed packet size %4d in %.8e (s)\n",
                            packetSize[ip], tdif.count());
                }
                else if (job == 1 && downFactor == 6)
                {
                    fprintf(stdout,
                            "Passed decimate random in %.8e (s)\n",
                            tdif.count());
                }
            } // Loop on packets 
        } // Loop on deterministic/random
    } // Loop on different downsampling factors
    free(x);
}
//============================================================================//
void read_decimate(const int nq, std::vector<double> *xdecim)
{
    xdecim->resize(0);
    std::string fileName = "data/decimate_" + std::to_string(nq) + ".txt";
    std::ifstream dfile;
    dfile.open(fileName);
    if (dfile.is_open())
    { 
        std::string line;
        int n = 0;
        while (std::getline(dfile, line))
        {
            n = n + 1;
        }
        dfile.close();
        dfile.open(fileName);
        xdecim->reserve(n);
        while (std::getline(dfile, line))
        {
            double xval;
            sscanf(line.c_str(), "%lf\n", &xval);
            xdecim->push_back(xval);
        }
        dfile.close();
    } 
}
void lowpassFilterThenDownsample(const int npts, const int nfir,
                                 const int downFactor, const double x[], 
                                 std::vector<double> *y,
                                 const bool lremovePhase)
{
    // Design an equivalent lowpass filter
    auto order = nfir - 1; 
    auto r = 1.0/static_cast<double> (downFactor);
    auto fir = RTSeis::Utilities::FilterDesign::FIR::FIR1Lowpass(order, r,
                           RTSeis::Utilities::FilterDesign::FIRWindow::HAMMING);
    // Initialize the lowpass single rate filter and filter the signal
    FIRFilter firFilter;
    auto taps = fir.getFilterTaps();
    firFilter.initialize(taps.size(), taps.data(),
                         RTSeis::ProcessingMode::POST_PROCESSING,
                         RTSeis::Precision::DOUBLE);
    // Create a downsampler
    Downsample downsample;
    downsample.initialize(downFactor,
                          RTSeis::ProcessingMode::POST_PROCESSING,
                          RTSeis::Precision::DOUBLE); 
    if (lremovePhase)
    {
        int groupDelay = static_cast<int> (taps.size()/2);
        int nptsPad = npts + groupDelay;
        std::vector<double> xpad(nptsPad, 0);
        std::vector<double> xfilt(nptsPad, 0);
        std::copy(x, x+npts, xpad.begin());
        auto xfiltPtr = xfilt.data();
        firFilter.apply(nptsPad, xpad.data(), &xfiltPtr);
        // Now downsample starting at the group delay
        auto nywork = downsample.estimateSpace(npts); // Filter last npts
        y->resize(nywork);
        std::fill(y->begin(), y->end(), 0.0);
        auto yPtr = y->data();
        int ny;
        downsample.apply(npts, xfiltPtr+groupDelay,
                         nywork, &ny, &yPtr); 
        assert(ny == nywork);
         
    }
    else
    {

    }

}

int readTextFile(int *npts, double *xPtr[], const std::string fileName)
{
    *xPtr = nullptr;
    *npts = 0;
    char line[64];
    double *x = nullptr;
    FILE *fl = fopen(fileName.c_str(), "r");
    *npts = 0;
    while (fscanf(fl, "%s", line) != EOF)
    {
        *npts = *npts + 1;
    }
    rewind(fl);
    if (*npts < 1)
    {
        fprintf(stderr, "%s", "No data points in file\n");
        return EXIT_FAILURE;
    }
    x = static_cast<double *> (calloc(static_cast<size_t> (*npts),
                              sizeof(double)));
    for (int i=0; i<*npts; i++)
    {
        memset(line, 0, 64*sizeof(char));
        fgets(line, 64, fl);
        sscanf(line, "%lf\n", &x[i]);
    }
    fclose(fl);
    *xPtr = x;
    return EXIT_SUCCESS;
} 

}
