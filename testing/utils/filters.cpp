#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <complex>
#include <vector>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/filterImplementations/downsample.hpp"
#include "rtseis/utilities/filterImplementations/iirFilter.hpp"
#include "rtseis/utilities/filterImplementations/iiriirFilter.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#include "rtseis/utilities/filterImplementations/multiRateFIRFilter.hpp"
#include "rtseis/utilities/filterImplementations/medianFilter.hpp"
#include "rtseis/utilities/filterImplementations/sosFilter.hpp"
#include "rtseis/utilities/filterImplementations/enums.hpp"
#include "rtseis/log.h"
#include "utils.hpp"
#include <ipps.h>

using namespace RTSeis::Utilities::FilterImplementations;
static int readTextFile(int *npts, double *xPtr[],
                        const std::string fileName = "utils/data/gse2.txt");
static int filters_downsample_test(const int npts, const double x[]);
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
    ierr = filters_downsample_test(npts, x);
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
//============================================================================//
int filters_iirFilter_test(const int npts, const double x[],
                           const std::string fileName1,
                           const std::string fileName2)
{
    fprintf(stdout, "Testing IIR filter...\n");
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
    int ierr = readTextFile(&npref, &yref1, fileName1);
    if (ierr != 0 || npts != npref)
    {
        RTSEIS_ERRMSG("%s", "Failed to load reference data");
        return EXIT_FAILURE;
    }
    ierr = readTextFile(&npref, &yref2, fileName2);
    if (ierr != 0 || npts != npref)
    {
        RTSEIS_ERRMSG("%s", "Failed to load reference data 2");
        return EXIT_FAILURE;
    }
    // Compute the zero-phase IIR filter alternative
    IIRFilter iir;
    IIRFilter iir_slow;
    ierr = iir.initialize(nb, b, na, a,
                          RTSeis::ProcessingMode::POST_PROCESSING,
                          RTSeis::Precision::DOUBLE,
                          IIRDFImplementation::DF2_FAST);
    if (ierr != 0)
    {    
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        return EXIT_FAILURE;
    }
    ierr = iir_slow.initialize(nb, b, na, a,
                               RTSeis::ProcessingMode::POST_PROCESSING,
                               RTSeis::Precision::DOUBLE,
                               IIRDFImplementation::DF2_SLOW);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize slow filter");
        return EXIT_FAILURE;
    }
    double *yref = new double[npts];
    double *yref_slow = new double[npts];
    double *y1 = new double[npts];
    double *y2 = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    ierr = iir.apply(npts, x, y1);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to filter signal");
        return EXIT_FAILURE;
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double error = 0;
    for (int i=0; i<npts; i++)
    {
        error = std::max(std::pow(y1[i] - yref1[i], 2), error);
    }
    error = error/static_cast<double> (npts);
    if (error > 3.e-2)
    {
        RTSEIS_ERRMSG("Failed matlab test1: %e\n", error);
        return EXIT_FAILURE;
    }
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Fast reference solution 1 computation time %.8lf (s)\n",
            tdif.count());
    std::copy(y1, y1+npts, yref);
    // Repeat for slow
    timeStart = std::chrono::high_resolution_clock::now();
    ierr = iir_slow.apply(npts, x, y2);
    timeEnd = std::chrono::high_resolution_clock::now();
    error = 0; 
    for (int i=0; i<npts; i++) 
    {
        error = std::max(std::pow(y2[i] - yref1[i], 2), error);
    }
    error = error/static_cast<double> (npts);
    if (error > 1.e-10)
    {
        RTSEIS_ERRMSG("Failed matlab test1 slow: %e\n", error);
        return EXIT_FAILURE;
    }
    tdif = timeEnd - timeStart;
    fprintf(stdout, "Slow reference solution 1 computation time %.8lf (s)\n",
            tdif.count());
    std::copy(y2, y2+npts, yref_slow);
    // Try a lower order filter - Note the difference in accuracy.
    // I think the filter is a little unstable.
    IIRFilter iir2;
    ierr = iir2.initialize(nb2, b2, na2, a2,
                           RTSeis::ProcessingMode::POST_PROCESSING,
                           RTSeis::Precision::DOUBLE,
                           IIRDFImplementation::DF2_FAST);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        return EXIT_FAILURE;
    }
    timeStart = std::chrono::high_resolution_clock::now();
    ierr = iir2.apply(npts, x, y2);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to filter signal");
        return EXIT_FAILURE;
    }
    timeEnd = std::chrono::high_resolution_clock::now();
    error = 0;
    for (int i=0; i<npts; i++)
    {
        error = std::max(std::pow(y2[i] - yref2[i], 2), error); 
    }
    if (error > 1.e-12)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute filter 2");
        return EXIT_FAILURE;
    }
    tdif = timeEnd - timeStart;
    fprintf(stdout, "Fast reference solution 2 computation time %.8lf (s)\n",
            tdif.count());
    // Do a real-time test of high-order filter
    IIRFilter iirrt;
    ierr = iirrt.initialize(nb, b, na, a,
                            RTSeis::ProcessingMode::REAL_TIME,
                            RTSeis::Precision::DOUBLE,
                            IIRDFImplementation::DF2_FAST);
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
                ierr = iir.apply(nptsPass, &x[nxloc], &y1[nxloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply iir filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
            }
            iir.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            for (int i=0; i<npts; i++)
            {
                if (std::abs(yref[i] - y1[i]) > 1.e-12)
                {
                    RTSEIS_ERRMSG("Failed to compute reference soln %d %lf %lf",
                                  i, yref[i], y1[i]);
                    return EXIT_FAILURE;
                }
            }
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
    ierr = iir_slow.initialize(nb, b, na, a,
                               RTSeis::ProcessingMode::REAL_TIME,
                               RTSeis::Precision::DOUBLE,
                               IIRDFImplementation::DF2_SLOW);
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
                ierr = iir.apply(nptsPass, &x[nxloc], &y1[nxloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply iir filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
            }
            iir.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            for (int i=0; i<npts; i++) 
            {
                if (std::abs(yref_slow[i] - y1[i]) > 1.e-12)
                {
                    RTSEIS_ERRMSG("Failed to compute reference soln %d %lf %lf",
                                  i, yref_slow[i], y1[i]);
                    return EXIT_FAILURE;
                }
            }
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
    return EXIT_SUCCESS;
}
//============================================================================//
int filters_iiriirFilter_test(const int npts, const double x[],
                              const std::string fileName)
{
    fprintf(stdout, "Testing IIRIIR filter...\n");
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
    int ierr = readTextFile(&npref, &yref, fileName);
    if (ierr != 0 || npts != npref)
    {   
        RTSEIS_ERRMSG("%s", "Failed to load reference data");
        return EXIT_FAILURE;
    } 
    // Compute the zero-phase IIR filter alternative
    IIRIIRFilter iiriir;
    ierr = iiriir.initialize(nb, b, na, a, RTSeis::Precision::DOUBLE);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        return EXIT_FAILURE;
    }
    double *y = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    ierr = iiriir.apply(npts, x, y); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply zero-phase IIR filter");
        return EXIT_FAILURE;
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double error = 0;
    for (int i=0; i<npts; i++)
    {
        error = error + std::pow(y[i] - yref[i], 2);
    }
    error = std::sqrt(error)/static_cast<double> (npts);
    if (error > 3.e-2)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute reference solution");
        return EXIT_FAILURE;
    }
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    delete[] y;
    free(yref);
    // Do the Intel test
#ifdef IPPS_H__
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
    ierr = iiriir.initialize(ORDER+1, &pTaps[0], ORDER+1, &pTaps[ORDER+1]);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        return EXIT_FAILURE;
    }
    ierr = iiriir.apply(LEN, xi, yi);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return EXIT_FAILURE;
    }
    for (int i=0; i<LEN; i++)
    {
        if (std::abs(yiRef[i] - yi[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Failed ref test %lf %lf", yiRef[i], yi[i]);
        }
    }
    ippsFree(f1);
    ippsFree(f2);
    ippsFree(pTaps);
    ippsFree(pBufIIRIIR);
    delete[] xi;
    delete[] yi;
    delete[] yiRef;
#endif 
    return EXIT_SUCCESS;
}
//============================================================================//
int filters_firMRFilter_test(const int npts, const double x[],
                             const std::string fileName)
{
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
    int ierr = readTextFile(&npref, &yref, fileName);
    if (ierr != 0)
    {   
        RTSEIS_ERRMSG("%s", "Failed to load reference data");
        return EXIT_FAILURE;
    }
    // Initialize filter
    //double gain = static_cast<double> (upFactor);
    //ippsMulC_64f(b, gain, bgain, nb); // Fix gain
    MultiRateFIRFilter firmr;
    ierr = firmr.initialize(upFactor, downFactor, nb, b,
                            RTSeis::ProcessingMode::POST_PROCESSING,
                            RTSeis::Precision::DOUBLE);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize fir filter");
        return EXIT_FAILURE;
    }
    // Estimate space in output
    int nest = firmr.estimateSpace(npts);
    if (nest != 4500)
    {
        RTSEIS_WARNMSG("%s", "Anticipating 4500 points in comparison");
    }
    int ncomp = std::min(nest, npref);
    // Filter 
    double *y = new double[npts];
    int ny;
    auto timeStart = std::chrono::high_resolution_clock::now();
    ierr = firmr.apply(npts, x, npts, &ny, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return EXIT_FAILURE;
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double error = 0;
    for (int i=0; i<ncomp; i++)
    {
        error = std::max(error, std::abs(y[i] - yref[i]));
    }
    if (error > 1.e-10)
    {
        RTSEIS_ERRMSG("Failed to apply upfirndn post-proc %.10e\n", error); 
        return EXIT_FAILURE;
    }
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Repeat for real-time 
    MultiRateFIRFilter firmrRT;
    ierr = firmrRT.initialize(upFactor, downFactor, nb, b,
                              RTSeis::ProcessingMode::REAL_TIME,
                              RTSeis::Precision::DOUBLE); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize real-time filter");
        return EXIT_FAILURE;
    }
    firmr = firmrRT; 
    std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                 1000, 1024, 1200, 2048, 4000, 4096, 5000});
    for (int job=0; job<2; job++)
    {   
        for (size_t ip=0; ip<packetSize.size(); ip++)
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
                ierr = firmr.apply(nptsPass, &x[nxloc], nwork, &nyDec, &y[nyloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply firmr filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
                nyloc = nyloc + nyDec;
            }
            firmr.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            if (nyloc != ny)
            {
                RTSEIS_ERRMSG("Failed resampling %d %d",  nyloc, ny);
                return EXIT_FAILURE;
            }
            for (int i=0; i<nyloc; i++) 
            {
                if (std::abs(yref[i] - y[i]) > 1.e-8)
                {
                    RTSEIS_ERRMSG("Failed to compute reference soln %d %lf %lf",
                                  i, y[i], yref[i]);
                    return EXIT_FAILURE;
                }
            }
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
    return EXIT_SUCCESS;
}
//============================================================================//
int filters_firFilter_test(const int npts, const double x[],
                           const std::string fileName)
{
    fprintf(stdout, "Testing FIR filter...\n");
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
    int ierr = readTextFile(&npref, &yref, fileName);
    if (ierr != 0 || npts != npref)
    {
        RTSEIS_ERRMSG("%s", "Failed to load reference data");
        return EXIT_FAILURE;
    }
    // Make a post-processing solution
    FIRFilter fir;
    ierr = fir.initialize(nb, b,
                          RTSeis::ProcessingMode::POST_PROCESSING,
                          RTSeis::Precision::DOUBLE,
                          FIRImplementation::DIRECT);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize FIR filter");
        return EXIT_FAILURE;
    }
    double *y = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    ierr = fir.apply(npts, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply FIR filter");
        return EXIT_FAILURE;
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    for (int i=0; i<npts; i++)
    {
        if (std::abs(y[i] - yref[i]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to apply reference FIR filter %lf %lf",
                          yref[i], y[i]);
            return EXIT_FAILURE;
        }
    }
    std::copy(y, y+static_cast<size_t> (npts), yref);
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Do packetized tests
    FIRFilter firrt;
    ierr = firrt.initialize(nb, b,
                          RTSeis::ProcessingMode::REAL_TIME,
                          RTSeis::Precision::DOUBLE,
                          FIRImplementation::DIRECT);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        return EXIT_FAILURE;
    }
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
                ierr = fir.apply(nptsPass, &x[nxloc], &y[nxloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply fir filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
            }
            fir.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            for (int i=0; i<npts; i++)
            {
                if (std::abs(yref[i] - y[i]) > 1.e-8)
                {
                    RTSEIS_ERRMSG("Failed to compute reference soln %d %lf %lf",
                                  i, y[i], yref[i]);
                    return EXIT_FAILURE;
                }
            }
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
    return EXIT_SUCCESS;
}
//============================================================================//
int filters_sosFilter_test(const int npts, const double x[],
                           const std::string fileName)
{
    fprintf(stdout, "Testing SOS filter...\n");
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
    int ierr = sos.initialize(ns, bs7, as7,
                              RTSeis::ProcessingMode::POST_PROCESSING,
                              RTSeis::Precision::DOUBLE);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize sos");
        return EXIT_FAILURE;
    }
    ierr = sos.apply(40, impulse, y40);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return EXIT_FAILURE;
    }
    for (int i=0; i<40; i++)
    {
        if (std::abs(y40[i] - yref40[i]) > 1.e-8)
        {
            RTSEIS_ERRMSG("Impulse response failed %lf %lf", yref40[i], y40[i]);
            return EXIT_FAILURE;
        }
    }
    delete[] y40;
    delete[] impulse;
    // Load a reference solution
    double *yref = nullptr;
    int npref = 0;
    ierr = readTextFile(&npref, &yref, fileName);
    if (ierr != 0 || npts != npref)
    {
        RTSEIS_ERRMSG("%s", "Failed to load reference data");
        return EXIT_FAILURE;
    }
    ns = 4;
    const double bs[12] = {0.000401587491686,  0.000803175141692,  0.000401587491549,
                           1.000000000000000, -2.000000394412897,  0.999999999730209,
                           1.000000000000000,  1.999999605765104,  1.000000000341065,
                           1.000000000000000, -1.999999605588274,  1.000000000269794};
    const double as[12] = {1.000000000000000, -1.488513049541281,  0.562472929601870,
                           1.000000000000000, -1.704970593447777,  0.792206889942566,
                           1.000000000000000, -1.994269533089365,  0.994278822534674,
                           1.000000000000000, -1.997472946622339,  0.997483252685326};
    ierr = sos.initialize(ns, bs, as,
                          RTSeis::ProcessingMode::POST_PROCESSING,
                          RTSeis::Precision::DOUBLE);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        return EXIT_FAILURE;
    } 
    double *y = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    ierr = sos.apply(npts, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute post-processing soln");
        return EXIT_FAILURE;
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    for (int i=0; i<npts; i++)
    {
        if (std::abs(y[i] - yref[i]) > 1.e-8)
        {
            RTSEIS_ERRMSG("Failed to compute sample %d %lf %lf",
                          i, yref[i], y[i]);
            return EXIT_FAILURE;
        }
    } 
    std::copy(y, y+static_cast<size_t> (npts), yref);
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Do packetized tests 
    SOSFilter sosrt;
    ierr = sosrt.initialize(ns, bs, as,
                            RTSeis::ProcessingMode::REAL_TIME,
                            RTSeis::Precision::DOUBLE);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        return EXIT_FAILURE;
    }
    sos = sosrt;
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
                ierr = sos.apply(nptsPass, &x[nxloc], &y[nxloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply sos filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
            }
            sos.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            for (int i=0; i<npts; i++)
            {
                if (std::abs(yref[i] - y[i]) > 1.e-8)
                {
                    RTSEIS_ERRMSG("Failed to compute reference soln %d %lf %lf",
                                  i, y[i], yref[i]);
                    return EXIT_FAILURE;
                }
            }
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
    return EXIT_SUCCESS;
    
}
//============================================================================//
int filters_medianFilter_test(const int npts, const double x[],
                              const std::string fileName)
{
    fprintf(stdout, "Testing median filter...\n");
    double xin[8] = {1, 2, 127, 4, 5, 0, 7, 8};
    double y8[8];
    double yref3[8] = {1, 2, 4, 5, 4, 5, 7, 7}; // Matlab soln; IPP has edge effect
    double yref5[8] = {1, 2, 4, 4, 5, 5, 5, 0};
    int ierr;
    MedianFilter median;
    ierr = median.initialize(3,
                             RTSeis::ProcessingMode::POST_PROCESSING,
                             RTSeis::Precision::DOUBLE);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize filter");
        return EXIT_FAILURE;
    }
    ierr = median.apply(8, xin, y8); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to apply filter");
        return EXIT_FAILURE;
    }
    for (int i=1; i<8-1; i++)
    {
        if (std::abs(y8[i+1] - yref3[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Failed test %lf %lf", y8[i+1], yref3[i]);
            return EXIT_FAILURE;
        }
    }
    ierr = median.initialize(5,
                             RTSeis::ProcessingMode::POST_PROCESSING,
                             RTSeis::Precision::DOUBLE);
    ierr = median.apply(8, xin, y8);
    for (int i=2; i<8-2; i++)
    {
        if (std::abs(y8[i+2] - yref5[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Failed test %lf %lf", y8[i+2], yref5[i]);
            return EXIT_FAILURE;
        }
    }
    // Load a reference solution
    double *yref = nullptr;
    int npref;
    ierr = readTextFile(&npref, &yref, fileName);
    if (ierr != 0 || npref - 11/2 != npts)
    {
        RTSEIS_ERRMSG("%s", "Failed to load reference data");
        return EXIT_FAILURE;
    }
    median.initialize(11,
                      RTSeis::ProcessingMode::POST_PROCESSING,
                      RTSeis::Precision::DOUBLE);
    auto timeStart = std::chrono::high_resolution_clock::now();
    double *y = new double[npts];
    ierr = median.apply(npts, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute reference solution");
        return EXIT_FAILURE;
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    for (int i=0; i<npts; i++)
    {
        if (std::abs(y[i] - yref[i]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute reference soln %d %lf %lf",
                          i, y[i], yref[i]);
            return EXIT_FAILURE;
        }
    }
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Now do the packetized tests
    MedianFilter medianrt = median;
    medianrt.initialize(11,
                        RTSeis::ProcessingMode::REAL_TIME,
                        RTSeis::Precision::DOUBLE);
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
                ierr = medianrt.apply(nptsPass, &x[nxloc], &y[nxloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply median filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
            }
            medianrt.resetInitialConditions();
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            for (int i=0; i<npts; i++)
            {
                if (std::abs(yref[i] - y[i]) > 1.e-10)
                {
                    RTSEIS_ERRMSG("Failed to compute reference soln %d %lf %lf",
                                  i, y[i], yref[i]);
                    return EXIT_FAILURE;
                }
            }
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
    return EXIT_SUCCESS;
}
//============================================================================//
int filters_downsample_test(const int npts, const double x[])
{
    const int nq = 7;
    RTSeis::Precision precision = RTSeis::Precision::DOUBLE;
    // Call this in post-processing for a couple different decimation rates
    srand(10245);
    Downsample downsample;
    double *y = static_cast<double *>
                (calloc(static_cast<size_t> (npts), sizeof(double)));
    double *yref = static_cast<double *>
                   (calloc(static_cast<size_t> (npts), sizeof(double)));
    int ierr = 0;
    fprintf(stdout, "Testing downsampler...\n");
    for (int iq=1; iq<nq+1; iq++)
    {
        // Do a post-processing test
        memset(y, 0, static_cast<size_t> (npts)*sizeof(double));
        memset(yref, 0, static_cast<size_t> (npts)*sizeof(double));
        ierr = downsample.initialize(iq,
                                     RTSeis::ProcessingMode::POST_PROCESSING,
                                     precision); 
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to intiialized downsample");
            return EXIT_FAILURE;
        }
        auto timeStart = std::chrono::high_resolution_clock::now();
        int ny;
        ierr = downsample.apply(npts, x, npts, &ny, y);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to call downsampler");
            return EXIT_FAILURE;
        }
        auto timeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> tdif = timeEnd - timeStart;
        // Do a manual downsample
        int j = 0;
        for (int i=0; i<npts; i=i+iq)
        {
            if (std::abs(y[j] - x[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("%s", "Post-processing downsample failed");
                return EXIT_FAILURE;
            }
            j = j + 1;
        }
        if (j != ny)
        {
            RTSEIS_ERRMSG("%s", "Incorrect number of output points");
            return EXIT_FAILURE;
        } 
        downsample.clear();
        fprintf(stdout,
                "Post-processing execution time for nq=%d is %.8lf (s)\n",
                iq, tdif.count());
        // Make a copy of the correct answer
        int nyref = ny;
        for (int iy=0; iy<ny; iy++){yref[iy] = y[iy];}
        // Do a real-time test
        ierr = downsample.initialize(iq,
                                     RTSeis::ProcessingMode::REAL_TIME,
                                     precision);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to intiialized downsample");
            return EXIT_FAILURE;
        }
        std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512,
                                     1000, 1024, 1200, 2048, 4000, 4096, 5000});
        for (size_t ip=0; ip<packetSize.size(); ip++)
        {
            timeStart = std::chrono::high_resolution_clock::now();
            int nxloc = 0;
            int nyloc = 0;
            while (nxloc < npts)
            {
                int nptsPass = std::min(packetSize[ip], npts - nxloc);
                int nyDec = 0;
                ierr = downsample.apply(nptsPass, &x[nxloc],
                                        npts+1-nyloc, &nyDec, &y[nyloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("Failed to apply downsampler for iq=%d", iq);
                    return EXIT_FAILURE; 
                }
                nxloc = nxloc + nptsPass;
                nyloc = nyloc + nyDec;
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            downsample.resetInitialConditions();
            if (nyloc != nyref)
            {
                RTSEIS_ERRMSG("Failed fixed packet size test %d %d",
                              nyloc, nyref);
                return EXIT_FAILURE;
            }
            for (int iy=0; iy<nyref; iy++)
            {
                if (std::abs(y[iy] - yref[iy]) > 1.e-10)
                {
                    RTSEIS_ERRMSG("Failed fixed packet size test %lf %lf",
                                  y[iy], yref[iy]);
                    return EXIT_FAILURE;
                }
            }
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
            ierr = downsample.apply(nptsPass, &x[nxloc],
                                    npts+1-nyloc, &nyDec, &y[nyloc]);
            if (ierr != 0)
            {
                RTSEIS_ERRMSG("Failed to apply downsampler for iq=%d", iq);
                return EXIT_FAILURE;
            }
            nxloc = nxloc + nptsPass;
            nyloc = nyloc + nyDec;
            packetLen = std::max(1, packetLen + rand()%50 - 25);
        }
        if (nyloc != nyref)
        {
            RTSEIS_ERRMSG("Failed fixed packet size test %d %d", nyloc, nyref);
            return EXIT_FAILURE;
        }
        for (int iy=0; iy<nyref; iy++)
        {
            if (std::abs(y[iy] - yref[iy]) > 1.e-10)
            {
                RTSEIS_ERRMSG("%s", "Failed fixed packet size test");
                return EXIT_FAILURE;
            }
        }
        timeEnd = std::chrono::high_resolution_clock::now();
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
    fprintf(stdout, "Passed downsampler test\n");
    return EXIT_SUCCESS;
}
//============================================================================//
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
        RTSEIS_ERRMSG("%s", "No data points in file\n");
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
