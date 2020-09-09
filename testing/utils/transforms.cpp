#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <cstring>
#include <cassert>
#include <string>
#include <cfloat>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <complex>
#include <vector>
#include <ipps.h>
#include <fftw/fftw3.h>
#include "rtseis/utilities/transforms/enums.hpp"
#include "rtseis/utilities/transforms/dftRealToComplex.hpp"
#include "rtseis/utilities/transforms/dft.hpp"
#include "rtseis/utilities/transforms/hilbert.hpp"
#include "rtseis/utilities/transforms/envelope.hpp"
#include "rtseis/utilities/transforms/firEnvelope.hpp"
#include "rtseis/utilities/transforms/welch.hpp"
#include "rtseis/utilities/transforms/slidingWindowRealDFTParameters.hpp"
#include "rtseis/utilities/transforms/slidingWindowRealDFT.hpp"
#include "rtseis/utilities/transforms/utilities.hpp"
#include "rtseis/utilities/transforms/wavelets/morlet.hpp"
#include "rtseis/utilities/transforms/continuousWavelet.hpp"
#include "rtseis/utilities/windowFunctions.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis::Utilities::Transforms;

void readEnvelopeFile(const std::string &fileName,
                      std::vector<double> *x,
                      std::vector<double> *upRef,
                      std::vector<double> *loRef);
void readTextFile(const std::string &fileName,
                  std::vector<double> *x);
int fft(const int nx, const std::complex<double> *x, 
        const int ny, std::complex<double> *y);
int ifft(const int nx, const std::complex<double> *x, 
         const int ny, std::complex<double> *y);
int rfft(const int nx, double x[], const int n,
         const int ny, std::complex<double> y[]);
//int irfft(const int nx, const std::complex<double> x[],
//          const int n, double y[]);

/*
int rtseis_test_utils_transforms(void)
{
    const std::string dataDir = "data/"; 
    int ierr = transforms_nextPowerOfTwo_test();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed nextpow2 test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed nextPowerOfTwo test");

    ierr = transforms_phase_test();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed phase test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed phase test");

    ierr = transforms_unwrap_test();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed to unwrap phase");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed unwrap test");

    ierr = transforms_test_dftr2c();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute rdft");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed rdft test");

    ierr = transforms_test_dft();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute dft");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed dft test");

    ierr = transforms_test_hilbert();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hilbert");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed hilbert test");

    ierr = transforms_test_envelope(dataDir + "/envelopeChirpReference.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute envelope");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed envelope test");

    ierr = transforms_test_firEnvelope(dataDir + "/envelopeChirpReference300.txt",
                                       dataDir + "/envelopeChirpReference301.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute FIR envelope");
    }
    RTSEIS_INFOMSG("%s", "Passed FIR envelope test");

    return EXIT_SUCCESS;
}
*/

//int transforms_nextPowerOfTwo_test(void)
TEST(UtilitiesTransforms, NextPow2)
{
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(0), 1);
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(1), 1);
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(2), 2);
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(3), 4);
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(4), 4);
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(5), 8);
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(1200), 2048);
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(120000), 131072);
    EXPECT_EQ(DFTUtilities::nextPowerOfTwo(131072), 131072);
}

//int transforms_unwrap_test(void)
TEST(UtilitiesTransforms, Unwrap)
{
    int n = 23;
    std::vector<double> p{0, -1.5728, -1.5747, -1.5772, -1.5790,
                          -1.5816, -1.5852, -1.5877, -1.5922,
                          -1.5976, -1.6044, -1.6129, -1.6269,
                          -1.6512, -1.6998, -1.8621,  1.7252,
                           1.6124,  1.5930,  1.5916,  1.5708,
                           1.5708,  1.5708};
    std::vector<double> qref{0.00000,  -1.57280,  -1.57470,  -1.57720,
                            -1.57900,  -1.58160,  -1.58520,  -1.58770,
                            -1.59220,  -1.59760,  -1.60440,  -1.61290,
                            -1.62690,  -1.65120,  -1.69980,  -1.86210,
                            -4.557985307179586,  -4.670785307179586,
                            -4.690185307179586,  -4.691585307179587,
                            -4.712385307179586,  -4.712385307179586,
                            -4.712385307179586};
    std::vector<double> qref2{1.000000000000000,  -3.710385307179586,
                             -3.708485307179586,  -3.705985307179586,
                             -3.704185307179586,  -3.701585307179586,
                             -3.697985307179586,  -3.695485307179586,
                             -3.690985307179586,  -3.685585307179586,
                             -3.678785307179586,  -3.670285307179586,
                             -3.656285307179586,  -3.631985307179586,
                             -3.583385307179586,  -3.421085307179586,
                             -0.725200000000000,  -0.612400000000000,
                             -0.593000000000000,  -0.591600000000000,
                             -0.570800000000000,  -0.570800000000000,
                             -0.570800000000000};
    // Should run with default tol = M_PI
    std::vector<double> q;
    EXPECT_NO_THROW(q = DFTUtilities::unwrap(p));
    double emax;
    ippsNormDiff_Inf_64f(q.data(), qref.data(), n, &emax);
    ASSERT_LE(emax, 1.e-10);
    // Switch tolerance to 90 degrees and adjust p
    for (auto i=0; i<n; ++i){p[i] =-p[i] + 1;}
    double tol = M_PI/2;
    EXPECT_NO_THROW(q = DFTUtilities::unwrap(p, tol));
    ippsNormDiff_Inf_64f(q.data(), qref2.data(), n, &emax);
    ASSERT_LE(emax, 1.e-10);
}

//int transforms_phase_test(void)
TEST(UtilitiesTransforms, Phase)
{
    const int n = 7;
    std::vector<double> tr{-0.785398163397448, 0.463647609000806,
                          -0.3217505543966, 0.244978663126864,
                           0, 1.570796326794897, 0};
    std::vector<double> tr2 = tr;
    std::vector<std::complex<double>> z(n);
    z[0] = std::complex<double> (1, -1);
    z[1] = std::complex<double> (2, +1);
    z[2] = std::complex<double> (3, -1);
    z[3] = std::complex<double> (4, +1);
    z[4] = std::complex<double> (0, +0);
    z[5] = std::complex<double> (0, +1);
    z[6] = std::complex<double> (1, +0);
    std::vector<double> angle;
    EXPECT_NO_THROW(angle = DFTUtilities::phase(z));
    double emax;
    ippsNormDiff_Inf_64f(angle.data(), tr.data(), n, &emax);
    ASSERT_LE(emax, 1.e-10);
    // Repeat same test but get result in degrees
    bool lwantDeg = true;
    EXPECT_NO_THROW(angle = DFTUtilities::phase(z, lwantDeg));
    for (auto i=0; i<n; i++){tr2[i] = tr[i]*180.0/M_PI;}
    ippsNormDiff_Inf_64f(angle.data(), tr2.data(), n, &emax);
    ASSERT_LE(emax, 1.e-10);
}

TEST(UtilitiesTransforms, rfftFreqs)
{
    // Edge case
    auto rfft1 = DFTUtilities::realToComplexDFTFrequencies(1, 1.0);
    EXPECT_EQ(rfft1[0], 0.0);
    // Length 8 which is even
    auto rfft2 = DFTUtilities::realToComplexDFTFrequencies(8, 5.0);
    EXPECT_EQ(static_cast<int> (rfft2.size()), 5);
    for (auto i=0; i<5; ++i)
    {
        EXPECT_NEAR(rfft2[i], static_cast<double> (i)*0.025, 1.e-14);
    }
    // Length 13 which is odd
    auto rfft3 = DFTUtilities::realToComplexDFTFrequencies(13, 4.0);
    EXPECT_EQ(static_cast<int> (rfft3.size()), 7);
    for (auto i=0; i<7; ++i)
    {
        EXPECT_NEAR(rfft3[i], i*0.019230769230769232, 1.e-13);
    } 
}

TEST(UtilitiesTransforms, fftShift)
{
    std::vector<double> f1({1});
    auto r1 = DFTUtilities::fftShift(f1);
    EXPECT_EQ(r1[0], f1[0]);

    std::vector<double> f2({1, 2});
    auto r2 = DFTUtilities::fftShift(f2);
    EXPECT_EQ(r2[0], f2[1]);
    EXPECT_EQ(r2[1], f1[0]);

    std::vector<double> f3({2, 3, 1});
    auto r3 = DFTUtilities::fftShift(f3);
    for (int i=0; i<3; ++i)
    {
        EXPECT_EQ(static_cast<double> (i+1), r3[i]);
    }

    std::vector<double> f9({5, 6, 7, 8, 9, 1, 2, 3, 4});
    auto r9 = DFTUtilities::fftShift(f9);
    for (int i=0; i<9; ++i)
    {
        EXPECT_EQ(static_cast<double> (i+1), r9[i]);
    }

    std::vector<double> f10({6, 7, 8, 9, 10, 1, 2, 3, 4, 5});
    auto r10 = DFTUtilities::fftShift(f10);
    for (int i=0; i<10; ++i)
    {
        EXPECT_EQ(static_cast<double> (i+1), r10[i]);
    }
}

//int transforms_test_dft()
TEST(UtilitiesTransforms, dft)
{
    std::vector<std::complex<double>> x5(5);

    x5[0] = std::complex<double> (0.293848340517710, 0.543040331839914);
    x5[1] = std::complex<double> (0.432658290043139, 0.507949811338692);
    x5[2] = std::complex<double> (0.638136660660825, 0.014141443357630);
    x5[3] = std::complex<double> (0.523377526876359, 0.431524418557386);
    x5[4] = std::complex<double> (0.146439036982003, 0.982433436513920);
    std::vector<std::complex<double>> y5ref(5);
    y5ref[0] = std::complex<double> ( 2.034459855080035, 2.479089441607542);
    y5ref[1] = std::complex<double> (-1.163477761938224, 0.303378415338802);
    y5ref[2] = std::complex<double> ( 0.302336705735982,-0.584079752541479);
    y5ref[3] = std::complex<double> ( 0.066216063691746,-0.465893684778681);
    y5ref[4] = std::complex<double> ( 0.229706840019013, 0.982707239573387);
    std::vector<std::complex<double>> y8ref(8);
    y8ref[0] = std::complex<double> ( 2.034459855080035,+2.479089441607542);
    y8ref[1] = std::complex<double> ( 0.761711158054026,-1.699508261045320);
    y8ref[2] = std::complex<double> (-0.121423890379805, 1.602051561829424);
    y8ref[3] = std::complex<double> ( 0.861724646436443,-0.531316766704685);
    y8ref[4] = std::complex<double> ( 0.122388221241041,+0.600140981815385);
    y8ref[5] = std::complex<double> (-0.438609664267351,-0.455551269624340);
    y8ref[6] = std::complex<double> (-0.274274675942418,+1.420613088162984);
    y8ref[7] = std::complex<double> (-0.595188926080287,+0.928803878678324);

    std::vector<std::complex<double>> y5(5), y8(8), y5fftw(5),
                                      x5inv(5), x5invref(5);
    fft(x5.size(), x5.data(), y5fftw.size(), y5fftw.data());
    ifft(y5ref.size(), y5ref.data(), x5invref.size(), x5invref.data());

    DFT<double> dft;
    EXPECT_NO_THROW(dft.initialize(x5.size()));
    EXPECT_TRUE(dft.isInitialized());
    ASSERT_EQ(dft.getInverseTransformLength(), 5);
    ASSERT_EQ(dft.getTransformLength(), 5);
    std::complex<double> *yPtr = y5.data();
    EXPECT_NO_THROW(dft.forwardTransform(5, x5.data(), 5, &yPtr)); //y5.data()));
    std::complex<double> *xPtr = x5inv.data();
    EXPECT_NO_THROW(dft.inverseTransform(5, y5.data(), 5, &xPtr)); //x5inv.data()));
    // Check forward transform
    double emax = 0;
    for (auto i=0; i<static_cast<int> (y5ref.size()); i++)
    {
        emax = std::max(emax, std::abs(y5ref[i] - y5[i]));
        if (std::abs(y5ref[i] - y5[i]) > 1.e-12)
        {
            fprintf(stderr, "Failed to compute dft %d %e",
                    i, std::abs(y5ref[i] - y5[i]));
        }
    }
    ASSERT_LE(emax, 1.e-12);
    // Check that inverse transform is covered
    emax = 0;
    for (auto i=0; i<static_cast<int> (x5inv.size()); i++)
    {
        emax = std::max(emax, std::abs(x5inv[i] - x5[i]));
        emax = std::max(emax, std::abs(x5invref[i] - x5inv[i]));
        if (std::abs(x5inv[i] - x5[i]) > 1.e-12 ||
            std::abs(x5invref[i] - x5inv[i]) > 1.e-12)
        {
           fprintf(stderr, "Failed to compute idft %d %e %e",
                   i, std::abs(x5[i] - x5inv[i]),
                   std::abs(x5invref[i] - x5inv[i]));
        }
    }
    ASSERT_LE(emax, 1.e-12);

    // Try padding
    EXPECT_NO_THROW(dft.initialize(8));
    ASSERT_EQ(dft.getInverseTransformLength(), 8);
    ASSERT_EQ(dft.getTransformLength(), 8);
    yPtr = y8.data();
    EXPECT_NO_THROW(dft.forwardTransform(5, x5.data(), 8, &yPtr)); //y8.data()));
    emax = 0;
    for (auto i=0; i<static_cast<int> (y8.size()); i++)
    {
        emax = std::max(emax, std::abs(y8[i] - y8ref[i]));
        if (std::abs(y8[i] - y8ref[i]) > 1.e-12)
        {
            fprintf(stderr, "Failed to compute dft %d %e",
                    i, std::abs(y8[i] - y8ref[i]));
        }
    }
    ASSERT_LE(emax, 1.e-12);
    // Do a larger test
    int niter = 50;
    int np0 = 12001;
    std::vector<std::complex<double>> x(np0+1), y(np0+1),
                                      yref(np0+1), xinv(np0+1);
    for (auto i=0; i<np0+1; i++)
    {
        x[i] = std::complex<double> (static_cast<double> (rand())/RAND_MAX,
                                     static_cast<double> (rand())/RAND_MAX);
    }
    for (auto j=0; j<2; j++)
    { 
        // Compute reference solution with FFTw
        auto npts = np0 + j;
        fft(npts, x.data(), npts, yref.data()); 
        // Try again with IPP
        dft.clear();
        EXPECT_NO_THROW(dft.initialize(npts));
        yPtr = y.data();
        EXPECT_NO_THROW(dft.forwardTransform(npts, x.data(), npts, &yPtr)); //y.data()));
        xPtr = xinv.data();
        EXPECT_NO_THROW(dft.inverseTransform(npts, y.data(), npts, &xPtr)); //xinv.data()));
        emax  = 0;
        double emaxi = 0;
        for (auto i=0; i<npts; i++)
        {
            emax = std::max(emax, std::abs(y[i] - yref[i]));
            emaxi= std::max(emaxi,std::abs(x[i] - xinv[i]));
        }
        ASSERT_LE(emax, 5.e-11);
        ASSERT_LE(emaxi, 1.e-12);
        // Stress test
        if (j == 1)
        {
            DFT<double> dftStress;
            dftStress.initialize(npts);
            auto timeStart = std::chrono::high_resolution_clock::now();
            for (auto i=0; i<niter; i++)
            {
                yPtr = y.data();
                dftStress.forwardTransform(npts, x.data(), npts, &yPtr); //y.data());
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            emax = 0;
            for (auto i=0; i<npts; i++)
            {
                emax  = std::max(emax,  std::abs(y[i] - yref[i]));
            }
            ASSERT_LE(emax, 5.e-11);
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            fprintf(stdout, "Average DFT time %.8lf (s)\n",
                    tdif.count()/static_cast<double>(niter));
 
            timeStart = std::chrono::high_resolution_clock::now();
            for (auto i=0; i<niter; i++)
            {
                xPtr = x.data();
                dftStress.inverseTransform(npts, y.data(), npts, &xPtr); //xinv.data());
            }
            timeEnd = std::chrono::high_resolution_clock::now();
            emaxi = 0;
            for (auto i=0; i<npts; i++)
            {   
                emaxi = std::max(emaxi, std::abs(x[i] - xinv[i]));
            }
            ASSERT_LE(emaxi, 1.e-12);
            tdif = timeEnd - timeStart;
            fprintf(stdout, "Average inverse DFT time %.8lf (s)\n",
                    tdif.count()/static_cast<double>(niter));
        }
    }
}


//int transforms_test_dftr2c(void)
TEST(UtilitiesTransforms, dftr2c)
{
    int niter = 50;
    int np0 = 12001;
    double *x = new double[np0+1];
    for (auto i=0; i<np0+1; i++)
    {
        x[i] = static_cast<double> (rand())/RAND_MAX;
    }

    for (auto j=0; j<2; j++)
    {
        auto npts = np0 + j;
        // Compute a DFT w/ FFTw
        auto lendft = npts/2 + 1;
        std::complex<double> *zrefDFT = new std::complex<double>[lendft];
        ASSERT_EQ(rfft(npts, x, npts, lendft, zrefDFT), 0);
        // Compute an FFT w/ FFTw
        auto np2 = DFTUtilities::nextPowerOfTwo(npts);
        auto lenfft = np2/2 + 1;
        std::complex<double> *zrefFFT = new std::complex<double>[lenfft];
        ASSERT_EQ(rfft(npts, x, np2, lenfft, zrefFFT), 0);
        // Initialize the DFT
        DFTRealToComplex<double> dft; 
        EXPECT_NO_THROW(dft.initialize(npts,
                        FourierTransformImplementation::DFT));
        ASSERT_EQ(dft.getTransformLength(), lendft);
        std::complex<double> *z = new std::complex<double>[lendft];
        EXPECT_NO_THROW(dft.forwardTransform(npts, x, lendft, &z));
        double error = 0; 
        for (auto i=0; i<lendft; i++)
        {
            error = std::max(error, std::abs(z[i] - zrefDFT[i]));
        }
        ASSERT_LE(error, 1.e-11);
        // Inverse DFT
        auto npout = dft.getInverseTransformLength(); 
        ASSERT_EQ(npout, npts);
        double *xinv = new double[npout];
        EXPECT_NO_THROW(dft.inverseTransform(lendft, z, npout, &xinv));
        error = 0;
        for (auto i=0; i<npts; i++)
        {
            error = std::max(error, std::abs(x[i] - xinv[i]));
        }
        ASSERT_LE(error, 1.e-10);
        delete[] xinv;
        // Stress test it
        if (j == 1)
        {
            auto timeStart = std::chrono::high_resolution_clock::now();
            for (auto kiter=0; kiter<niter; kiter++)
            {
                EXPECT_NO_THROW(dft.forwardTransform(npts, x, lendft, &z));
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            error = 0;
            for (auto i=0; i<lendft; i++)
            {
                error = std::max(error, std::abs(z[i] - zrefDFT[i]));
            }
            ASSERT_LE(error, 1.e-11);
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            fprintf(stdout, "Average real DFT time %.8lf (s)\n",
                    tdif.count()/static_cast<double>(niter)); 
 
        }
        delete[] z;
        // Redo this for an FFT
        EXPECT_NO_THROW(dft.initialize(npts,
                        FourierTransformImplementation::FFT));
        ASSERT_EQ(dft.getTransformLength(), lenfft);
        z = new std::complex<double>[lenfft];
        EXPECT_NO_THROW(dft.forwardTransform(npts, x, lenfft, &z));
        error = 0;  
        for (auto i=0; i<lenfft; i++)
        {
            error = std::max(error, std::abs(z[i] - zrefFFT[i]));
        }
        ASSERT_LE(error, 1.e-12);
        npout = dft.getInverseTransformLength(); 
        ASSERT_EQ(npout, np2);
        xinv = new double[npout];
        EXPECT_NO_THROW(dft.inverseTransform(lenfft, z, npout, &xinv));
        error = 0;
        for (auto i=0; i<npts; i++)
        {
            error = std::max(error, std::abs(x[i] - xinv[i]));
        }
        ASSERT_LE(error, 1.e-10);
        delete[] xinv;
        // Stress test it
        if (j == 1)
        {
            auto timeStart = std::chrono::high_resolution_clock::now();
            for (auto kiter=0; kiter<niter; kiter++)
            {
                EXPECT_NO_THROW(dft.forwardTransform(npts, x, lenfft, &z));
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            error = 0;  
            for (auto i=0; i<lenfft; i++)
            {
                error = std::max(error, std::abs(z[i] - zrefFFT[i]));
            }
            ASSERT_LE(error, 1.e-12);
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            fprintf(stdout, "Average real FFT time %.8lf (s)\n",
                    tdif.count()/static_cast<double>(niter)); 
 
        }
        delete[] z;
        // Free reference space 
        delete[] zrefFFT;
        delete[] zrefDFT;
    }
    delete[] x;
}

TEST(UtilitiesTransforms, Hilbert)
{
    std::vector<std::complex<double>> h10(10), h11(11);
    h10[0] = std::complex<double> (0, 5.50552768);
    h10[1] = std::complex<double> (1,-0.64983939);
    h10[2] = std::complex<double> (2,-0.64983939);
    h10[3] = std::complex<double> (3,-2.10292445);
    h10[4] = std::complex<double> (4,-2.10292445);
    h10[5] = std::complex<double> (5,-2.10292445);
    h10[6] = std::complex<double> (6,-2.10292445);
    h10[7] = std::complex<double> (7,-0.64983939);
    h10[8] = std::complex<double> (8,-0.64983939);
    h10[9] = std::complex<double> (9, 5.50552768);

    h11[0]  = std::complex<double> ( 0, 6.42868554);
    h11[1]  = std::complex<double> ( 1,-0.52646724);
    h11[2]  = std::complex<double> ( 2,-0.23284074);
    h11[3]  = std::complex<double> ( 3,-2.42253531);
    h11[4]  = std::complex<double> ( 4,-1.77987433);
    h11[5]  = std::complex<double> ( 5,-2.93393585);
    h11[6]  = std::complex<double> ( 6,-1.77987433);
    h11[7]  = std::complex<double> ( 7,-2.42253531);
    h11[8]  = std::complex<double> ( 8,-0.23284074);
    h11[9]  = std::complex<double> ( 9,-0.52646724);
    h11[10] = std::complex<double> (10, 6.42868554);
    Hilbert<double> hilbert;    
    std::vector<std::complex<double>> h;
    std::vector<double> x;
    int n = 10;
    x.reserve(n+1);
    x.resize(n);
    for (auto i=0; i<n; i++){x[i] = static_cast<double> (i);}
    EXPECT_NO_THROW(hilbert.initialize(n));
    EXPECT_EQ(hilbert.getTransformLength(), n);
    h.resize(n);
    std::complex<double> *hPtr = h.data();
    EXPECT_NO_THROW(hilbert.transform(n, x.data(), &hPtr));
    double emax = 0;
    for (auto i=0; i<static_cast<int> (x.size()); i++)
    {
        emax = std::max(emax, std::abs(h[i] - h10[i]));
        if (std::abs(h[i] - h10[i]) > 1.e-8)
        {
            fprintf(stderr, "Failed hilbert %d %lf\n", 
                    i, std::abs(h[i] - h10[i]));
        }
    }
    ASSERT_LE(emax, 1.e-8);
    // Test copy constructor with n = 11
    n = 11;
    Hilbert<double> hilbert11;
    EXPECT_NO_THROW(hilbert11.initialize(n));
    hilbert = hilbert11;
    EXPECT_EQ(hilbert.getTransformLength(), n);
    x.resize(n);
    h.resize(n);
    for (auto i=0; i<static_cast<int> (x.size()); i++)
    {
        x[i] = static_cast<double> (i);
    }
    hPtr = h.data();
    EXPECT_NO_THROW(hilbert.transform(n, x.data(), &hPtr));
    emax = 0;
    for (auto i=0; i<static_cast<int>(x.size()); i++)
    {
        emax = std::max(emax, std::abs(h[i] - h11[i]));
        if (std::abs(h[i] - h11[i]) > 1.e-8)
        {
            fprintf(stderr, "Failed hilbert %d %lf\n", 
                    i, std::abs(h[i] - h11[i]));
        }
    }
    ASSERT_LE(emax, 1.e-8);
}

TEST(UtilitiesTransforms, Envelope)
{
    const std::string fileName = "data/envelopeChirpReference.txt";
    Envelope<double> envelope;
    // Do an edge case analytically
    int n = 1;
    double x1[1] = {9};
    double upRef1[1] = {9};
    double loRef1[1] = {9};
    double up1[1], lo1[1];
    EXPECT_NO_THROW(envelope.initialize(n));
    double *upPtr = up1; 
    double *loPtr = lo1;
    EXPECT_NO_THROW(envelope.transform(1, x1, &upPtr, &loPtr));
    ASSERT_LE(std::abs(upRef1[0] - up1[0]), 1.e-15);
    ASSERT_LE(std::abs(loRef1[0] - lo1[0]), 1.e-15);
    // Load the data
    std::vector<double> x, upRef, loRef;
    readEnvelopeFile(fileName, &x, &upRef, &loRef);
    auto npts = static_cast<int> (x.size());
    EXPECT_EQ(npts, 4000);
    // Do a real test
    std::vector<double> yupper(npts), ylower(npts);
    Envelope<double> envelopeChirp;
    EXPECT_NO_THROW(envelopeChirp.initialize(npts));
    envelope = envelopeChirp; // Test copy assignment operator
    EXPECT_TRUE(envelope.isInitialized()); // Verify it's initialized
    double *yPtrUp = yupper.data();
    double *yPtrLo = ylower.data();
    EXPECT_NO_THROW(envelope.transform(npts, x.data(), &yPtrUp, &yPtrLo));
    double errorLower;
    double errorUpper;
    ippsNormDiff_L1_64f(upRef.data(), yupper.data(), npts, &errorUpper);
    ippsNormDiff_L1_64f(loRef.data(), ylower.data(), npts, &errorLower);
    ASSERT_LE(errorUpper, 1.e-9);
    ASSERT_LE(errorLower, 1.e-9);
}

TEST(UtilitiesTransforms, firEnvelope)
{
    const std::string fileName1 = "data/envelopeChirpReference300.txt";
    const std::string fileName2 = "data/envelopeChirpReference301.txt";
    // Try a simple case to test the post-processing filtering logic
    std::vector<double> xSimple{2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 0, -1,
                               -2, -3, -4, -3, -2, -1, 0};
    std::vector<double> ySimpleRef1{2.104854409188811, 3.054435023431979,
                                    4.036561092009176, 5.027493422155269,
                                    6.000000000000000, 5.027493422155269,
                                    4.036561092009176, 3.054435023431979,
                                    2.104854409188811, 1.469790661363077,
                                    2.104854409188811, 3.054435023431979,
                                    4.036561092009176, 5.027493422155269,
                                    6.000000000000000, 5.027493422155269,
                                    4.036561092009176, 3.054435023431979,
                                    2.104854409188811};
    std::vector<double> ySimpleRef2{2.652854154751004, 3.594503058524673,
                                    4.567414567447291, 5.633676421237045,
                                    5.633676421237045, 4.567414567447291,
                                    3.594503058524673, 2.652708722278031,
                                    1.855246233954807, 1.855246233954807,
                                    2.652708722278031, 3.594503058524673,
                                    4.567414567447291, 5.633676421237045,
                                    5.633676421237045, 4.567414567447291,
                                    3.594503058524673, 2.652854154751004,
                                    1.787936920537643};
    std::vector<double> ySimple(xSimple.size(), 0);
    FIREnvelope<double> env;
    EXPECT_NO_THROW(env.initialize(5, RTSeis::ProcessingMode::POST_PROCESSING));
    double *yPtr = ySimple.data();
    EXPECT_NO_THROW(env.transform(xSimple.size(),
                    xSimple.data(), &yPtr));
#ifdef __STDCPP_MATH_SPEC_FUNCS__ // __cplusplus > 201402L
    double tol = 1.e-13;
#else
    double tol = 1.e-5;
#endif
    double error;
    ippsNormDiff_L1_64f(ySimpleRef1.data(), ySimple.data(),
                        ySimple.size(), &error);
    ASSERT_LE(error, tol);

    EXPECT_NO_THROW(env.initialize(6, RTSeis::ProcessingMode::POST_PROCESSING));
    yPtr = ySimple.data();
    EXPECT_NO_THROW(env.transform(xSimple.size(), xSimple.data(), &yPtr));
    ippsNormDiff_L1_64f(ySimpleRef2.data(), ySimple.data(), 
                        ySimple.size(), &error);
    ASSERT_LE(error, tol);
    // Do a more substantial post-processing test
    std::vector<double> x, upRef300, loRef300, upRef301, loRef301;
    readEnvelopeFile(fileName1, &x, &upRef300, &loRef300);
    readEnvelopeFile(fileName2, &x, &upRef301, &loRef301);
    ASSERT_EQ(static_cast<int> (upRef300.size()), 4000);
    ASSERT_EQ(static_cast<int> (upRef301.size()), 4000);
    // Create with copy constructor
    FIREnvelope<double> env300;
    EXPECT_NO_THROW(env300.initialize(300,
                    RTSeis::ProcessingMode::POST_PROCESSING));
    env = env300;
    std::vector<double> up(x.size());
    yPtr = up.data();
    env.transform(x.size(), x.data(), &yPtr);
    ippsNormDiff_L1_64f(upRef300.data(), up.data(),  upRef300.size(), &error);
    ASSERT_LE(error/upRef300.size(), 1.e-8);
    // Test move constructor
    FIREnvelope<double> env301;
    EXPECT_NO_THROW(env301.initialize(301,
                    RTSeis::ProcessingMode::POST_PROCESSING));
    env = std::move(env301);
    yPtr = up.data();
    env.transform(x.size(), x.data(), &yPtr);
    ippsNormDiff_L1_64f(upRef301.data(), up.data(),  upRef301.size(), &error);
    ASSERT_LE(error/upRef301.size(), 1.e-8);
    // Remove the mean to make comparison easier
    double mean;
    ippsMean_64f(x.data(), x.size(), &mean);
    ippsSubC_64f_I(mean, x.data(), x.size());
    auto timeStart = std::chrono::high_resolution_clock::now();
    env.initialize(300, RTSeis::ProcessingMode::POST_PROCESSING);
    ippsZero_64f(upRef300.data(), upRef300.size()); 
    yPtr = upRef300.data();
    env.transform(x.size(), x.data(), &yPtr);
    env.initialize(301, RTSeis::ProcessingMode::POST_PROCESSING);
    ippsZero_64f(upRef301.data(), upRef301.size());
    yPtr = upRef301.data();
    env.transform(x.size(), x.data(), &yPtr);
    // Test the real-time component
    FIREnvelope<double> envrt300, envrt301;
    std::vector<double> up300(x.size()), up301(x.size());
    envrt300.initialize(300, RTSeis::ProcessingMode::REAL_TIME);
    envrt301.initialize(301, RTSeis::ProcessingMode::REAL_TIME);
    auto timeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference time: %.8e (s)\n", tdif.count());
    std::vector<int> packetSize({1, 2, 3, 16, 64, 100, 200, 512, 
                                 1000, 1024, 1200, 2048, 4000});
    int npts = x.size();
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
                const double *xptr = x.data() + nxloc;
                double *yptr300 = up300.data() + nxloc;
                double *yptr301 = up301.data() + nxloc;
                EXPECT_NO_THROW(envrt300.transform(nptsPass, xptr, &yptr300));
                EXPECT_NO_THROW(envrt301.transform(nptsPass, xptr, &yptr301));
                nxloc = nxloc + nptsPass;
            } // Loop on acquisition loop
            envrt300.resetInitialConditions();
            envrt301.resetInitialConditions();
            timeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            // Need to deal with the impulse response
            int groupDelay = 301/2;
            int ncomp = npts - groupDelay;
            ippsNormDiff_L1_64f(upRef300.data(), up300.data()+groupDelay,
                                ncomp, &error);
            ASSERT_LE(error, 1.e-10);
            ippsNormDiff_L1_64f(upRef301.data(), up301.data()+groupDelay,
                                ncomp, &error); 
            ASSERT_LE(error, 1.e-10);
            if (job == 0)
            {
                fprintf(stdout,
                        "Passed envfir filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed envfir random in %.8e (s)\n", tdif.count());
            }
        }
    }
}

TEST(UtilitiesTransforms, SlidingWindowRealDFTParameters)
{
    SlidingWindowRealDFTParameters parameters;
    int nSamples = 1040;
    int windowLength = 88;
    int overlapLength = 22;
    int dftLength = 124;
    std::vector<double> windowRef(88);
    std::vector<double> window;
    EXPECT_NO_THROW(parameters.setNumberOfSamples(nSamples));
    EXPECT_FALSE(parameters.isValid());
    EXPECT_NO_THROW(parameters.setWindow(windowLength,
                                         SlidingWindowType::HANN));
    EXPECT_TRUE(parameters.isValid());
    EXPECT_NO_THROW(parameters.setDFTLength(dftLength));
    EXPECT_NO_THROW(parameters.setNumberOfSamplesInOverlap(overlapLength));
    EXPECT_NO_THROW(parameters.setDetrendType(SlidingWindowDetrendType::REMOVE_MEAN));
    EXPECT_NO_THROW(parameters.setPrecision(RTSeis::Precision::FLOAT));
 
    EXPECT_EQ(parameters.getNumberOfSamples(), nSamples);
    EXPECT_EQ(parameters.getWindowType(), SlidingWindowType::HANN);
    EXPECT_EQ(parameters.getWindowLength(), windowLength); 
    EXPECT_EQ(parameters.getNumberOfSamplesInOverlap(), overlapLength);
    EXPECT_EQ(parameters.getDetrendType(), SlidingWindowDetrendType::REMOVE_MEAN);
    EXPECT_EQ(parameters.getPrecision(), RTSeis::Precision::FLOAT);
    window = parameters.getWindow();
    double *windowData = windowRef.data();
    RTSeis::Utilities::WindowFunctions::hann(windowLength, &windowData);
    double emax = 0;
    ippsNormDiff_Inf_64f(window.data(), windowRef.data(),
                         windowLength, &emax);
    EXPECT_LE(emax, 1.e-15);
    // set a custom window
    EXPECT_NO_THROW(parameters.setWindow(windowLength, window.data()));
    EXPECT_EQ(parameters.getWindowType(), SlidingWindowType::CUSTOM);
    window = parameters.getWindow();
    EXPECT_EQ(parameters.getWindowLength(), windowLength);
    EXPECT_EQ(parameters.getDFTLength(), windowLength); // Override DFT length
    EXPECT_EQ(parameters.getNumberOfSamplesInOverlap(), 0); // Override overlap
    ippsNormDiff_Inf_64f(window.data(), windowRef.data(),
                         windowLength, &emax);
    EXPECT_LE(emax, 1.e-15);

    // Test copy operator
    parameters.setNumberOfSamplesInOverlap(overlapLength);
    parameters.setDFTLength(dftLength);
    SlidingWindowRealDFTParameters pcopy(parameters); 
    EXPECT_EQ(pcopy.getNumberOfSamples(), nSamples);
    EXPECT_EQ(pcopy.getWindowType(), SlidingWindowType::CUSTOM);
    EXPECT_EQ(pcopy.getWindowLength(), windowLength);
    EXPECT_EQ(pcopy.getNumberOfSamplesInOverlap(), overlapLength);
    EXPECT_EQ(pcopy.getDetrendType(), SlidingWindowDetrendType::REMOVE_MEAN);
    EXPECT_EQ(pcopy.getPrecision(), RTSeis::Precision::FLOAT);
    window = pcopy.getWindow();
    ippsNormDiff_Inf_64f(window.data(), windowRef.data(),
                         windowLength, &emax);
    EXPECT_LE(emax, 1.e-15);
}

TEST(UtilitiesTransforms, Spectrogram)
{
    SlidingWindowRealDFT sdft;
    // Load the chirp
    const std::string signalFileName = "data/spectrogramNoisyChirpSignal.txt";
    const std::string specFileName = "data/spectrogramNoisyChirp.txt";
    std::ifstream signalFile(signalFileName);
    std::string line;
    std::vector<double> sig;
    sig.reserve(2048);
    while (std::getline(signalFile, line)) 
    {
        double sample;
        std::sscanf(line.c_str(), "%lf\n", &sample);
        sig.push_back(sample);
    }
    EXPECT_EQ(static_cast<int> (sig.size()), 2048);
    // Load the spectrogram
    std::ifstream specFile(specFileName);
    std::vector<std::complex<double>> spec;
    spec.reserve(129*199);
    int ic = 0;
    int jc = 0;
    while (std::getline(specFile, line))
    {
        ic = ic + 1;
        if (ic == 200)
        {
            jc = jc + 1;
            ic = 0;
            continue;
        }
        double re, im;
        std::sscanf(line.c_str(), "%lf %lf\n", &re, &im);
        spec.push_back(std::complex<double> (re, im));
    } 
    EXPECT_EQ(jc, 129);
    EXPECT_EQ(static_cast<int> (spec.size()), 129*199);
    // Create a Kaiser window
    int nSamples = 2048;
    int nSamplesPerSegment = 63;
    int windowLength = nSamplesPerSegment;
    int dftLength = 256;
    int nSamplesInOverlap = windowLength - 10;
    std::vector<double> window(windowLength);
    const double beta = 17;
    double *kdata = window.data();
    EXPECT_NO_THROW(
       RTSeis::Utilities::WindowFunctions::kaiser(windowLength, &kdata, beta));
    // Set the parameters
    SlidingWindowRealDFTParameters parameters;
    EXPECT_NO_THROW(parameters.setNumberOfSamples(nSamples));
    EXPECT_NO_THROW(parameters.setWindow(windowLength, window.data()));
    EXPECT_NO_THROW(parameters.setDFTLength(dftLength));
    EXPECT_NO_THROW(parameters.setNumberOfSamplesInOverlap(nSamplesInOverlap));
    EXPECT_NO_THROW(
       parameters.setDetrendType(SlidingWindowDetrendType::REMOVE_NONE));
    EXPECT_NO_THROW(parameters.setPrecision(RTSeis::Precision::DOUBLE));
    EXPECT_TRUE(parameters.isValid());
    // Check the outputs
    EXPECT_EQ(parameters.getNumberOfSamples(), nSamples);
    EXPECT_EQ(parameters.getWindowLength(), windowLength);
    std::vector<double> windowBack = parameters.getWindow();
    double resmax = 0;
    ippsNormDiff_Inf_64f(window.data(), windowBack.data(),
                         windowLength, &resmax);
    EXPECT_LE(resmax, 1.e-15);
    EXPECT_EQ(parameters.getDFTLength(), dftLength);
    EXPECT_EQ(parameters.getNumberOfSamplesInOverlap(), nSamplesInOverlap);
    EXPECT_EQ(parameters.getDetrendType(), SlidingWindowDetrendType::REMOVE_NONE);
    EXPECT_EQ(parameters.getPrecision(), RTSeis::Precision::DOUBLE);
    // Initialize
    EXPECT_NO_THROW(sdft.initialize(parameters));
/*
    EXPECT_NO_THROW(sdft.initialize(nSamples,
                    nSamplesPerSegment,
                    dftLength,
                    nSamplesInOverlap,
                    windowLength,
                    window.data(),
                    SlidingWindowDetrendType::REMOVE_NONE,
                    RTSeis::Precision::DOUBLE));
*/
    EXPECT_EQ(sdft.getNumberOfSamples(), nSamples);
    EXPECT_EQ(sdft.getNumberOfFrequencies(), 129); // Number of rows
    EXPECT_EQ(sdft.getNumberOfTransformWindows(), 199); // Number of columns
    EXPECT_NO_THROW(sdft.transform(sig.size(), sig.data()));
    // Loop on the windows
    resmax = 0;
    for (auto i=0; i<sdft.getNumberOfTransformWindows(); ++i)
    {
        const std::complex<double> *cptr = nullptr;
        EXPECT_NO_THROW(cptr = sdft.getTransform64f(i));
        for (auto j=0; j<sdft.getNumberOfFrequencies(); j++)
        {
            int indx = sdft.getNumberOfTransformWindows()*j + i;
            double res = std::abs(cptr[j] - spec[indx]);
            resmax = std::max(resmax, res); 
            if (res > 1.e-7)
            {
                fprintf(stderr, "%d %+lf %+lfi, %+lf %+lfi, %lf\n",
                        indx, std::real(cptr[j]), std::imag(cptr[j]),
                        std::real(spec[indx]), std::imag(spec[indx]), res);
            }
        }
    }
    ASSERT_LE(resmax, 1.e-7);
}

TEST(UtilitiesTransforms, Welch)
{
    // Dirty trick - I need to read a 3 column text file so I can use envelope
    const std::string fileName1 = "data/welchTest1.txt";
    const std::string fileName2 = "data/welchTest2.txt";
    std::vector<double> fRef1, fRef2, powerRef1, powerRef2, psdRef1, psdRef2;
    readEnvelopeFile(fileName1, &fRef1, &powerRef1, &psdRef1);
    readEnvelopeFile(fileName2, &fRef2, &powerRef2, &psdRef2);
    // Setup the parameters 
    SlidingWindowRealDFTParameters parameters;
    double samplingRate = 10e3;
    int nSamples = 10000;
    double freq = 1234.0;
    int nWindowLength = 1024;
    SlidingWindowType windowType = SlidingWindowType::HANN;
    int fftLength = nWindowLength; // Default
    int nSamplesInOverlap = nWindowLength/2; // Good choice for Welch
    
    // Create the sine wave whose frequency is 1234 Hz
    std::vector<double> xSineSignal(nSamples);
    for (auto i=0; i<nSamples; ++i)
    {
        auto time = static_cast<double> (i)/samplingRate;
        xSineSignal[i] = std::sin(2.0*M_PI*time*freq);
    }
    // Put the parameters into the parameters structure 
    EXPECT_NO_THROW(parameters.setNumberOfSamples(nSamples));
    EXPECT_NO_THROW(parameters.setWindow(nWindowLength, windowType));
    EXPECT_NO_THROW(parameters.setNumberOfSamplesInOverlap(nSamplesInOverlap));
    EXPECT_NO_THROW(parameters.setDetrendType(SlidingWindowDetrendType::REMOVE_MEAN));
    EXPECT_NO_THROW(parameters.setDFTLength(fftLength));
    EXPECT_TRUE(parameters.isValid());
    Welch welch;
    EXPECT_NO_THROW(welch.initialize(parameters, samplingRate));
    EXPECT_EQ(welch.getNumberOfSamples(), nSamples);
    int nFrequencies = welch.getNumberOfFrequencies(); 
    EXPECT_EQ(nFrequencies, fftLength/2 + 1); //nWindowLength/2 + 1);
    // Check the frequencies
    std::vector<double> frequencies(nFrequencies);
    double *freqPtr = frequencies.data();
    welch.getFrequencies(nFrequencies, &freqPtr);
    EXPECT_EQ(fRef1.size(), frequencies.size());
    double error = 0;
    ippsNormDiff_Inf_64f(fRef1.data(), frequencies.data(), nFrequencies,
                         &error);
    EXPECT_LE(error, 1.e-7);
    // Transform
    EXPECT_NO_THROW(welch.transform(nSamples, xSineSignal.data()));
    // Get and compare results
    std::vector<double> spectrum(nFrequencies); 
    double *sPtr = spectrum.data();
    EXPECT_NO_THROW(welch.getPowerSpectralDensity(nFrequencies, &sPtr));
    ippsNormDiff_Inf_64f(powerRef1.data(), spectrum.data(),
                         nFrequencies, &error);
    EXPECT_LE(error, 1.e-5);
    EXPECT_NO_THROW(welch.getPowerSpectrum(nFrequencies, &sPtr));
    ippsNormDiff_Inf_64f(psdRef1.data(), spectrum.data(),
                         nFrequencies, &error);
    EXPECT_LE(error, 1.e-5);
    // Now feed the program some awkward numbers
    nWindowLength = 1013;
    windowType = SlidingWindowType::HAMMING;
    fftLength = 1091;
    nSamplesInOverlap = 501;
    EXPECT_NO_THROW(parameters.setNumberOfSamples(nSamples));
    EXPECT_NO_THROW(parameters.setWindow(nWindowLength, windowType));
    EXPECT_NO_THROW(parameters.setNumberOfSamplesInOverlap(nSamplesInOverlap));
    EXPECT_NO_THROW(parameters.setDetrendType(SlidingWindowDetrendType::REMOVE_NONE));
    EXPECT_NO_THROW(parameters.setDFTLength(fftLength));
    EXPECT_TRUE(parameters.isValid());
    EXPECT_NO_THROW(welch.initialize(parameters, samplingRate));
    nFrequencies = welch.getNumberOfFrequencies();
    frequencies.resize(nFrequencies);
    EXPECT_EQ(fRef2.size(), frequencies.size());
    freqPtr = frequencies.data();
    welch.getFrequencies(nFrequencies, &freqPtr);
    ippsNormDiff_Inf_64f(fRef2.data(), frequencies.data(), nFrequencies,
                         &error);
    EXPECT_LE(error, 1.e-7);
    // Transform
    EXPECT_NO_THROW(welch.transform(nSamples, xSineSignal.data()));
    spectrum.resize(nFrequencies);
    sPtr = spectrum.data();
    EXPECT_NO_THROW(welch.getPowerSpectralDensity(nFrequencies, &sPtr));
    ippsNormDiff_Inf_64f(powerRef2.data(), spectrum.data(),
                         nFrequencies, &error);
    EXPECT_LE(error, 1.e-5);
    EXPECT_NO_THROW(welch.getPowerSpectrum(nFrequencies, &sPtr));
    ippsNormDiff_Inf_64f(psdRef2.data(), spectrum.data(),
                         nFrequencies, &error);
    EXPECT_LE(error, 1.e-5);
}

TEST(UtilitiesTransforms, CWT)
{
    // Read the signal
    std::vector<double> x;
    readTextFile("data/zwave_cwt_example.txt", &x);
    // Create morlet wavelet
    double omega0 = 6;
    Wavelets::Morlet morlet;
    EXPECT_NO_THROW(morlet.setParameter(omega0));
    // Initialize cwt
    ContinuousWavelet<double> cwt;
    int nSamples = 1008;
    EXPECT_EQ(nSamples, static_cast<int> (x.size()));
    double samplingRate = 100;
    double fmin = 1;
    double fmax = 40;
    int nf = 200;
    auto nScales = nf;
    double df = (fmax - fmin)/static_cast<int> (nf - 1);
    std::vector<double> scales(nScales, 0);
    for (int i=0; i<nf; ++i)
    {
        auto f = fmin + i*df;
        scales[i] = (omega0*samplingRate)/(2*M_PI*f); // Inversely propto f
    }
    //morlet.initialize(nSamples, morlet, samplingRate); 
    EXPECT_NO_THROW(cwt.initialize(nSamples, nScales, scales.data(),
                                   morlet, samplingRate));
    EXPECT_TRUE(cwt.isInitialized());
    EXPECT_EQ(cwt.getNumberOfSamples(), nSamples);
    EXPECT_EQ(cwt.getNumberOfScales(), nScales);
    EXPECT_NO_THROW(cwt.transform(x.size(), x.data()));
}

//============================================================================//
//                              Private functions                             //
//============================================================================//

// Reads and envelope solution file
void readEnvelopeFile(const std::string &fileName,
                      std::vector<double> *x,
                      std::vector<double> *upRef,
                      std::vector<double> *loRef)
{
    // Load the data
    x->resize(0);
    upRef->resize(0);
    loRef->resize(0);
    x->reserve(4000); 
    upRef->reserve(4000); 
    loRef->reserve(4000);
    std::ifstream textFile(fileName);
    std::string line;
    auto i = 0;
    while (std::getline(textFile, line))
    {
        double xVal, uVal, lVal;
        std::sscanf(line.c_str(), "%lf %lf %lf\n", &xVal, &uVal, &lVal);
        x->push_back(xVal);
        upRef->push_back(uVal);
        loRef->push_back(lVal);
        i = i + 1;
    }
#ifdef DEBUG
    auto npts = static_cast<int> (x->size());
    assert(npts == 4000);
#endif
}

// Reads a text file
void readTextFile(const std::string &fileName,
                  std::vector<double> *x)
{
    // Load the data
    x->resize(0);
    x->reserve(5000);
    std::ifstream textFile(fileName);
    std::string line;
    auto i = 0;
    while (std::getline(textFile, line))
    {
        double xVal;
        std::sscanf(line.c_str(), "%lf\n", &xVal);
        x->push_back(xVal);
        i = i + 1;
    }
}

int fft(const int nx, const std::complex<double> *x,
        const int ny, std::complex<double> *y)
{
    if (nx < 1 || ny < 1)
    {
        if (nx < 1){fprintf(stderr, "%s", "nx must be positive");}
        if (ny < 1){fprintf(stderr, "%s", "ny must be positive");}
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){fprintf(stderr, "%s", "x is NULL");}
        if (y == nullptr){fprintf(stderr, "%s", "y is NULL");}
        return -1; 
    }
    // Set space and make plan
    auto n = ny;
    size_t nbytes = sizeof(fftw_complex)*static_cast<size_t> (n);
    fftw_complex *in  = static_cast<fftw_complex *> (fftw_malloc(nbytes));
    auto inCopy = reinterpret_cast<std::complex<double> *> (in);
    memset(in, 0, nbytes);
    fftw_complex *out = reinterpret_cast<fftw_complex *> (y);
    fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    // Equal size transforms
    if (nx == ny) 
    {
        std::copy(x, x+n, inCopy); //std::memcpy(in, x, nbytes);
    }   
    // Truncate x to length of output array y
    else if (nx > ny) 
    {
        //size_t ncopy = sizeof(fftw_complex)*static_cast<size_t> (ny);
        std::copy(x, x+ny, inCopy); ////std::memcpy(in, x, ncopy);
    }   
    // Pad x to length of output array y
    else //if (nx < ny) 
    {
        //size_t ncopy = sizeof(fftw_complex)*static_cast<size_t> (nx);
        std::copy(x, x+nx, inCopy); //std::memcpy(in, x, ncopy);
    }
    // Transform
    fftw_execute(p);
    // Free plan and data
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_cleanup();
    return 0;
}

int ifft(const int nx, const std::complex<double> *x,
         const int ny, std::complex<double> *y)
{
    if (nx < 1 || ny < 1)
    {
        if (nx < 1){fprintf(stderr, "%s", "nx must be positive");}
        if (ny < 1){fprintf(stderr, "%s", "ny must be positive");}
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){fprintf(stderr, "%s", "x is NULL");}
        if (y == nullptr){fprintf(stderr, "%s", "y is NULL");}
        return -1; 
    }
    // Set space and make plan
    auto n = ny; 
    size_t nbytes = sizeof(fftw_complex)*static_cast<size_t> (n);
    fftw_complex *in  = static_cast<fftw_complex *> (fftw_malloc(nbytes));
    auto inCopy = reinterpret_cast<std::complex<double> *> (in);
    memset(in, 0, nbytes);
    fftw_complex *out = reinterpret_cast<fftw_complex *> (y); //(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(size_t) n);
    fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    // Equal size transforms
    if (nx == ny) 
    {
        std::copy(x, x+nx, inCopy); //std::memcpy(in, x, nbytes);
    }   
    // Truncate x to length of output array y
    else if (nx > ny) 
    {
        //size_t ncopy = sizeof(fftw_complex)*static_cast<size_t> (ny);
        std::copy(x, x+ny, inCopy); //std::memcpy(in, x, ncopy);
    }
    // Pad x to length of output array y
    else //if (nx < ny) 
    {
        //size_t ncopy = sizeof(fftw_complex)*static_cast<size_t> (nx);
        std::copy(x, x+nx, inCopy); //std::memcpy(in, x, ncopy);
    }
    // Transform
    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_cleanup();
    double xnorm = 1.0/static_cast<double> (ny);
    for (int i=0; i<ny; i++)
    {
        y[i] = xnorm*y[i];
    }

    return 0;
}

/*
int irfft(const int nx, const std::complex<double> x[],
          const int n, double y[])
{
    if (n < 1 || nx < 1)
    {
        if (n < 1){fprintf(stderr, "%s", "n must be positive");}
        if (nx < 1){fprintf(stderr, "%s", "nx must be positive");}
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){fprtinf(stderr, "%s", "x is NULL");}
        if (y == nullptr){fprintf(stderr, "%s", "y is NULL");}
        return -1;
    }
    int ntf = n/2 + 1;
    size_t nbytes = sizeof(fftw_complex)*static_cast<size_t> (ntf);
    fftw_complex *in = static_cast<fftw_complex *> (fftw_malloc(nbytes));
    double *out = y;
    fftw_plan p = fftw_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
    // Equal size transforms
    if (nx == ntf)
    {
        for (int i=0; i<nx; i++)
        {
            in[i][0] = std::real(x[i]);
            in[i][1] = std::imag(x[i]);
        }
    }
    // Truncate x to length of output array y
    else if (nx > ntf)
    {
        for (int i=0; i<ntf; i++)
        {
            in[i][0] = std::real(x[i]);
            in[i][1] = std::imag(x[i]);
        }
    }
    // Pad to length of transform
    {
        for (int i=0; i<nx; i++)
        {
            in[i][0] = std::real(x[i]);
            in[i][1] = std::imag(x[i]);
        }
        for (int i=nx; i<ntf; i++)
        {
            in[i][0] = 0;
            in[i][1] = 0;
        }
    }
    // Transform
    fftw_execute(p);
    // Copy and normalize
    double xnorm = 1/static_cast<double> (n);
    #pragma omp simd
    for (int i=0; i<n; i++){y[i] = xnorm*y[i];}
    // Clean up
    fftw_destroy_plan(p);  
    fftw_free(in);
    fftw_cleanup();
    return 0;
} 
*/

int rfft(const int nx, double x[], const int n,
         const int ny, std::complex<double> y[])
{
    if (n < 1 || nx < 1)
    {
        if (n < 1){fprintf(stderr, "%s", "n must be positive");}
        if (nx < 1){fprintf(stderr, "%s", "nx must be positive");}
        return -1;
    }
    int ntf = n/2 + 1;
    if (ny < ntf)
    {
        fprintf(stderr, "ny = %d must be at least %d", ny, ntf);
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){fprintf(stderr, "%s", "x is NULL");}
        fprintf(stderr, "%s", "y is NULL");
        return -1;
    }
    size_t nbytes = static_cast<size_t> (n)*sizeof(double);
    double *in = static_cast<double *> (fftw_malloc(nbytes));

    nbytes = sizeof(fftw_complex)*static_cast<size_t> (ntf);
    fftw_complex *out = reinterpret_cast<fftw_complex *> (y);
    fftw_plan p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
    // Equal size transforms
    if (nx == n)
    {
        std::copy(x, x+n, in); 
    }
    // Truncate x to length of output array y
    else if (nx > n)
    {
        std::copy(x, x+n, in);
    }
    // Pad x to length of output array y
    else // if (nx < n)
    {
        std::copy(x, x+nx, in);
        #pragma omp simd 
        for (int i=nx; i<n; i++){in[i] = 0;} 
    }
    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_cleanup();
    out = nullptr;
    return 0;
}

}
