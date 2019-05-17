#include <fstream>
#include <cstdio>
#include <cstdlib>
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
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/transforms/enums.hpp"
#include "rtseis/utilities/transforms/dftRealToComplex.hpp"
#include "rtseis/utilities/transforms/dft.hpp"
#include "rtseis/utilities/transforms/hilbert.hpp"
#include "rtseis/utilities/transforms/envelope.hpp"
#include "rtseis/utilities/transforms/firEnvelope.hpp"
#include "rtseis/utilities/transforms/utilities.hpp"
#include "rtseis/log.h"
#include "utils.hpp"

using namespace RTSeis::Utilities::Transforms;

int transforms_nextPow2_test();
int transforms_phase_test();
int transforms_unwrap_test();
int transforms_test_dftr2c();
int transforms_test_dft();
int transforms_test_hilbert();
int transforms_test_envelope(const std::string fileName = "data/envelopeChirp.txt");
int transforms_test_firEnvelope(const std::string fileName1 = "data/envelopeChirp300.txt",
                                const std::string fileName2 = "data/envelopeChirp301.txt");
void readEnvelopeFile(const std::string &fileName,
                      std::vector<double> *x,
                      std::vector<double> *upRef,
                      std::vector<double> *loRef);
int fft(const int nx, const std::complex<double> *x, 
        const int ny, std::complex<double> *y);
int ifft(const int nx, const std::complex<double> *x, 
         const int ny, std::complex<double> *y);
int rfft(const int nx, double x[], const int n,
         const int ny, std::complex<double> y[]);
int irfft(const int nx, const std::complex<double> x[],
          const int n, double y[]);

int rtseis_test_utils_transforms(void)
{
    const std::string dataDir = "data/"; 
    int ierr = transforms_nextPow2_test();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed nextpow2 test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed nextPow2 test");

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

int transforms_nextPow2_test(void)
{
    if (DFTUtilities::nextPow2(0) != 1)
    {   
        RTSEIS_ERRMSG("%s", "Failed 0 test");
        return EXIT_FAILURE;
    }   
    if (DFTUtilities::nextPow2(1) != 1)
    {
        RTSEIS_ERRMSG("%s", "Failed 1 test");
        return EXIT_FAILURE;
    }
    if (DFTUtilities::nextPow2(2) != 2)
    {
        RTSEIS_ERRMSG("%s", "Failed 2 test");
        return EXIT_FAILURE;
    }
    if (DFTUtilities::nextPow2(3) != 4)
    {
        RTSEIS_ERRMSG("%s", "Failed 3 test");
        return EXIT_FAILURE;
    }
    if (DFTUtilities::nextPow2(4) != 4)
    {
        RTSEIS_ERRMSG("%s", "Failed 4 test");
        return EXIT_FAILURE;
    }
    if (DFTUtilities::nextPow2(5) != 8)
    {
        RTSEIS_ERRMSG("%s", "Failed 8 test");
        return EXIT_FAILURE;
    }
    if (DFTUtilities::nextPow2(1200) != 2048)
    {
        RTSEIS_ERRMSG("%s", "Failed 1200 test");
        return EXIT_FAILURE;
    }
    if (DFTUtilities::nextPow2(120000) != 131072)
    {
        RTSEIS_ERRMSG("%s", "Failed 120000 test");
        return EXIT_FAILURE;
    }
    if (DFTUtilities::nextPow2(131072) != 131072)
    {
        RTSEIS_ERRMSG("%s", "Failed 131072 test");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int transforms_unwrap_test(void)
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
    try
    {
         q = DFTUtilities::unwrap(p);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    for (int i=0; i<n; i++)
    {
        if (std::abs(q[i] - qref[i]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to unwrap %lf %lf", q[i], qref[i]);
            return EXIT_FAILURE;
        }
    }
    // Switch tolerance to 90 degrees and adjust p
    for (int i=0; i<n; i++){p[i] =-p[i] + 1;}
    double tol = M_PI/2;
    try
    {
        q = DFTUtilities::unwrap(p, tol);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    for (int i=0; i<n; i++)
    {
        if (std::abs(q[i] - qref2[i]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to unwrap %lf %lf", q[i], qref2[i]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int transforms_phase_test(void)
{
    const int n = 7;
    std::vector<double> tr{-0.785398163397448, 0.463647609000806,
                          -0.3217505543966, 0.244978663126864,
                           0, 1.570796326794897, 0};
    std::vector<std::complex<double>> z(7);
    z[0] = std::complex<double> (1, -1);
    z[1] = std::complex<double> (2, +1);
    z[2] = std::complex<double> (3, -1);
    z[3] = std::complex<double> (4, +1);
    z[4] = std::complex<double> (0, +0);
    z[5] = std::complex<double> (0, +1);
    z[6] = std::complex<double> (1, +0);
    std::vector<double> angle;
    try
    {
        angle = DFTUtilities::phase(z);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    for (int i=0; i<n; i++)
    {
        if (std::abs(tr[i] - angle[i]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute angle %lf %lf", tr[i], angle[i]);
            return EXIT_FAILURE;
        } 
    } 

    bool lwantDeg = true;
    try
    {
        angle = DFTUtilities::phase(z, lwantDeg);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    for (int i=0; i<n; i++)
    {
        if (std::abs(tr[i]*180.0/M_PI - angle[i]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute angle %lf %lf", tr[i], angle[i]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS; 
}

int transforms_test_hilbert()
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

    Hilbert hilbert;    
    std::vector<std::complex<double>> h;
    std::vector<double> x;
    int n = 10;
    x.reserve(n+1);
    x.resize(n);
    for (int i=0; i<n; i++){x[i] = static_cast<double> (i);}
    try
    {
        hilbert.initialize(n);
    }
    catch (const std::invalid_argument &ia) 
    { 
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    assert(hilbert.getTransformLength() == n);
    h.resize(n);
    hilbert.transform(n, x.data(), h.data());
    for (auto i=0; i<x.size(); i++)
    {
        if (std::abs(h[i] - h10[i]) > 1.e-8)
        {
            RTSEIS_ERRMSG("Failed hilbert %d %lf\n", 
                          i, std::abs(h[i] - h10[i]));
            return EXIT_FAILURE;
        }
    }
    // Test copy constructor with n = 11
    n = 11;
    Hilbert hilbert11;
    try
    {
        hilbert11.initialize(n);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    hilbert = hilbert11;
    assert(hilbert.getTransformLength() == n);
    x.resize(n);
    h.resize(n);
    for (int i=0; i<x.size(); i++){x[i] = static_cast<double> (i);}
    hilbert.transform(n, x.data(), h.data());
    for (auto i=0; i<x.size(); i++)
    {
        if (std::abs(h[i] - h11[i]) > 1.e-8)
        {
            RTSEIS_ERRMSG("Failed hilbert %d %lf\n", 
                          i, std::abs(h[i] - h11[i]));
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int transforms_test_dft()
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
    try
    {
        DFT dft;
        dft.initialize(x5.size());
        assert(dft.getInverseTransformLength() == 5);
        assert(dft.isInitialized());
        //ASSERT_EQUAL(dft.getInverseTransformLength(), 5);
        assert(dft.getTransformLength() == 5);
        //ASSERT_EQUAL(dft.getTransformLengt(), 5);
        dft.forwardTransform(5, x5.data(), 5, y5.data()); 
        dft.inverseTransform(5, y5.data(), 5, x5inv.data());
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    // Check forward transform
    for (size_t i=0; i<y5ref.size(); i++)
    {
        if (std::abs(y5ref[i] - y5[i]) > 1.e-12)
        {
            RTSEIS_ERRMSG("Failed to compute dft %ld %e",
                           i, std::abs(y5ref[i] - y5[i]));
            return EXIT_FAILURE;
        }
    }
    // Check that inverse transform is covered
    for (size_t i=0; i<x5inv.size(); i++)
    {
        if (std::abs(x5inv[i] - x5[i]) > 1.e-12 ||
            std::abs(x5invref[i] - x5inv[i]) > 1.e-12)
        {
           RTSEIS_ERRMSG("Failed to compute idft %ld %e %e",
                          i, std::abs(x5[i] - x5inv[i]),
                          std::abs(x5invref[i] - x5inv[i]));
           return EXIT_FAILURE; 
        }
    }
    // Try padding
    try
    {
        DFT dft;
        dft.initialize(8);
        assert(dft.getInverseTransformLength() == 8); 
        //ASSERT_EQUAL(dft.getInverseTransformLength(), 8);
        assert(dft.getTransformLength() == 8); 
        //ASSERT_EQUAL(dft.getTransformLengt(), 8);
        dft.forwardTransform(5, x5.data(), 8, y8.data());
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    for (auto i=0; i<y8.size(); i++)
    {
        if (std::abs(y8[i] - y8ref[i]) > 1.e-12)
        {
            RTSEIS_ERRMSG("Failed to compute dft %d %e",
                           i, std::abs(y8[i] - y8ref[i]));
        }
    }
    // Do a larger test
    int niter = 50;
    int np0 = 12001;
    std::vector<std::complex<double>> x(np0+1), y(np0+1),
                                      yref(np0+1), xinv(np0+1);
    for (int i=0; i<np0+1; i++)
    {
        x[i] = std::complex<double> (static_cast<double> (rand())/RAND_MAX,
                                     static_cast<double> (rand())/RAND_MAX);
    }
    for (int j=0; j<2; j++)
    { 
        // Compute reference solution with FFTw
        int npts = np0 + j;
        fft(npts, x.data(), npts, yref.data()); 
        // Try again with IPP
        try
        {
            DFT dft;
            dft.initialize(npts);
            dft.forwardTransform(npts, x.data(), npts, y.data());
            dft.inverseTransform(npts, y.data(), npts, xinv.data()); 
        }
        catch (const std::invalid_argument &ia)
        {
            RTSEIS_ERRMSG("%s", ia.what());
            return EXIT_FAILURE;
        }
        double emax  = 0;
        double emaxi = 0;
        for (int i=0; i<npts; i++)
        { 
            emax  = std::max(emax,  std::abs(y[i] - yref[i]));
            emaxi = std::max(emaxi, std::abs(x[i] - xinv[i]));
        }
        if (emax > 1.e-12)
        {
            RTSEIS_ERRMSG("%s", "Forward transform failed");
        }
        if (emaxi > 1.e-12)
        {
            RTSEIS_ERRMSG("%s", "Inverse transform failed");
        }
        // Stress test
        if (j == 1)
        {
            DFT dft;
            dft.initialize(npts);
            auto timeStart = std::chrono::high_resolution_clock::now();
            for (int i=0; i<niter; i++)
            {
                dft.forwardTransform(npts, x.data(), npts, y.data());
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            emax = 0;
            for (int i=0; i<npts; i++)
            {
                emax  = std::max(emax,  std::abs(y[i] - yref[i]));
            }
            if (emax > 1.e-12)
            {
                RTSEIS_ERRMSG("%s", "Forward transform failed");
            }
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            fprintf(stdout, "Average DFT time %.8lf (s)\n",
                    tdif.count()/static_cast<double>(niter));
 
            timeStart = std::chrono::high_resolution_clock::now();
            for (int i=0; i<niter; i++)
            {
                dft.inverseTransform(npts, y.data(), npts, xinv.data());
            }
            timeEnd = std::chrono::high_resolution_clock::now();
            emaxi = 0;
            for (int i=0; i<npts; i++)
            {   
                emaxi = std::max(emaxi, std::abs(x[i] - xinv[i]));
            }
            if (emaxi > 1.e-12)
            {
                RTSEIS_ERRMSG("%s", "Inverse transform failed");
            }
            tdif = timeEnd - timeStart;
            fprintf(stdout, "Average inverse DFT time %.8lf (s)\n",
                    tdif.count()/static_cast<double>(niter));
        }
    }
    return EXIT_SUCCESS;
}

int transforms_test_dftr2c(void)
{
    int niter = 50;
    int np0 = 12001;
    double *x = new double[np0+1];
    for (int i=0; i<np0+1; i++)
    {
        x[i] = static_cast<double> (rand())/RAND_MAX;
    }

    for (int j=0; j<2; j++)
    {
        int npts = np0 + j;
        // Compute a DFT w/ FFTw
        int lendft = npts/2 + 1;
        std::complex<double> *zrefDFT = new std::complex<double>[lendft];
        int ierr = rfft(npts, x, npts, lendft, zrefDFT);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to compute rdft");
            return EXIT_FAILURE;
        }
        // Compute an FFT w/ FFTw
        int np2 = DFTUtilities::nextPow2(npts);
        int lenfft = np2/2 + 1;
        std::complex<double> *zrefFFT = new std::complex<double>[lenfft];
        ierr = rfft(npts, x, np2, lenfft, zrefFFT); 
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to compute rfft");
            return EXIT_FAILURE;
        }
        // Initialize the DFT
        DFTRealToComplex dft; 
        try
        {
            dft.initialize(npts,
                           FourierTransformImplementation::DFT,
                           RTSeis::Precision::DOUBLE); 
        }
        catch (const std::invalid_argument &ia)
        {
            RTSEIS_ERRMSG("%s", ia.what());
            return EXIT_FAILURE;
        }
        if (dft.getTransformLength() != lendft)
        {
            RTSEIS_ERRMSG("%s", "Inconsistent sizes");
            return EXIT_FAILURE;
        }
        std::complex<double> *z = new std::complex<double>[lendft];
        try
        {
            dft.forwardTransform(npts, x, lendft, z);
        }
        catch (const std::invalid_argument &ia)
        {
            RTSEIS_ERRMSG("%s", ia.what());
            return EXIT_FAILURE;
        }
        double error = 0; 
        for (int i=0; i<lendft; i++)
        {
            error = std::max(error, std::abs(z[i] - zrefDFT[i]));
        }
        if (error > 1.e-11)
        {
            RTSEIS_ERRMSG("Failed to compute dft %.10e", error);
            return EXIT_FAILURE;
        }
        // Inverse DFT
        int npout = dft.getInverseTransformLength(); 
        if (npout != npts)
        {
            RTSEIS_ERRMSG("%s", "Size inconsitency");
            return EXIT_FAILURE;
        }
        double *xinv = new double[npout];
        try
        {
            dft.inverseTransform(lendft, z, npout, xinv); 
        }
        catch (const std::invalid_argument &ia)
        {
            RTSEIS_ERRMSG("%s", ia.what());
            return EXIT_FAILURE;
        }
        error = 0;
        for (int i=0; i<npts; i++)
        {
            error = std::max(error, std::abs(x[i] - xinv[i]));
        }
        if (error > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute idft %.10e", error);
            return EXIT_FAILURE;
        } 
        delete[] xinv;
        // Stress test it
        if (j == 1)
        {
            auto timeStart = std::chrono::high_resolution_clock::now();
            for (int kiter=0; kiter<niter; kiter++)
            {
                try
                {
                    dft.forwardTransform(npts, x, lendft, z); 
                }
                catch (const std::exception &e)
                {
                    RTSEIS_ERRMSG("%s", e.what());
                    return EXIT_FAILURE;
                }
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            error = 0;
            for (int i=0; i<lendft; i++)
            {
                error = std::max(error, std::abs(z[i] - zrefDFT[i]));
            }
            if (error > 1.e-11)
            {
                RTSEIS_ERRMSG("Failed to compute dft %.10e", error);
                return EXIT_FAILURE;
            }
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            fprintf(stdout, "Average real DFT time %.8lf (s)\n",
                    tdif.count()/static_cast<double>(niter)); 
 
        }
        delete[] z;
        // Redo this for an FFT
        dft.initialize(npts,
                       FourierTransformImplementation::FFT,
                       RTSeis::Precision::DOUBLE);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to initialize fft");
            return EXIT_FAILURE;
        }
        if (dft.getTransformLength() != lenfft)
        {
            RTSEIS_ERRMSG("Inconsistent sizes %d %d",
                          dft.getTransformLength(), lenfft);
            return EXIT_FAILURE;
        }
        z = new std::complex<double>[lenfft];
        try
        {
            dft.forwardTransform(npts, x, lenfft, z); 
        }
        catch (const std::invalid_argument &ia)
        {
            RTSEIS_ERRMSG("%s", ia.what());
            return EXIT_FAILURE;
        }
        error = 0;  
        for (int i=0; i<lenfft; i++)
        {
            error = std::max(error, std::abs(z[i] - zrefFFT[i]));
        }
        if (error > 1.e-12)
        {
            RTSEIS_ERRMSG("Failed to compute fft %.10e", error);
            return EXIT_FAILURE;
        }
        npout = dft.getInverseTransformLength(); 
        if (npout != np2)
        {
            RTSEIS_ERRMSG("%s", "Size inconsitency");
            return EXIT_FAILURE;
        }
        xinv = new double[npout];
        try
        {
            dft.inverseTransform(lenfft, z, npout, xinv); 
        }
        catch (const std::invalid_argument &ia)
        {
            RTSEIS_ERRMSG("%s", ia.what());
            return EXIT_FAILURE;
        }
        error = 0;
        for (int i=0; i<npts; i++)
        {
            error = std::max(error, std::abs(x[i] - xinv[i]));
        }
        if (error > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute ifft %.10e", error);
            return EXIT_FAILURE;
        }
        delete[] xinv;
        // Stress test it
        if (j == 1)
        {
            auto timeStart = std::chrono::high_resolution_clock::now();
            for (int kiter=0; kiter<niter; kiter++)
            {
                try
                {
                    dft.forwardTransform(npts, x, lenfft, z);
                }
                catch (const std::invalid_argument &ia)
                {
                    RTSEIS_ERRMSG("%s", ia.what());
                    return EXIT_FAILURE;
                }
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            error = 0;  
            for (int i=0; i<lenfft; i++)
            {
                error = std::max(error, std::abs(z[i] - zrefFFT[i]));
            }
            if (error > 1.e-12)
            {
                RTSEIS_ERRMSG("Failed to compute fft %.10e", error);
                return EXIT_FAILURE;
            }
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
    return EXIT_SUCCESS;
}

int transforms_test_envelope(const std::string fileName)
{
    Envelope envelope;
    try
    {
        int n = 1;
        double x[1] = {9};
        double upRef[1] = {9};
        double loRef[1] = {9};
        double up[1], lo[1];
        envelope.initialize(n, RTSeis::Precision::DOUBLE);
        envelope.transform(1, x, up, lo);
        if (std::abs(upRef[0] - up[0]) > 1.e-15 ||
            std::abs(loRef[0] - lo[0]) > 1.e-15)
        {
            RTSEIS_ERRMSG("%s", "Failed envelope length 1");
            return EXIT_FAILURE;
        }
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("%s", e.what());
        return EXIT_FAILURE;
    }
    // Load the data
    std::vector<double> x, upRef, loRef;
    readEnvelopeFile(fileName, &x, &upRef, &loRef);
    int npts = static_cast<int> (x.size());
    if (npts != 4000)
    {
        RTSEIS_ERRMSG("Failed to load %s", fileName.c_str());
        return EXIT_FAILURE;
    }
    assert(npts == 4000);
    // Do a real test
    std::vector<double> yupper(npts), ylower(npts);
    Envelope envelopeChirp;
    try
    {
        envelopeChirp.initialize(npts, RTSeis::Precision::DOUBLE);
    }
    catch (std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return EXIT_FAILURE;
    }
    envelope = envelopeChirp; // Test copy assignment operator
    assert(envelope.isInitialized()); // Verify it's initialized
    try
    {
        envelope.transform(npts, x.data(), yupper.data(), ylower.data()); 
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("%s", e.what());
        return EXIT_FAILURE;
    }
    double errorLower;
    double errorUpper;
    ippsNormDiff_L1_64f(upRef.data(), yupper.data(), npts, &errorUpper);
    ippsNormDiff_L1_64f(loRef.data(), ylower.data(), npts, &errorLower);
    if (errorUpper > 1.e-9)
    { 
        RTSEIS_ERRMSG("Failed upper envelope %e", errorUpper);
        return EXIT_FAILURE;
    }
    if (errorLower > 1.e-9)
    {
        RTSEIS_ERRMSG("Failed lower envelope %e", errorLower);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int transforms_test_firEnvelope(const std::string fileName1,
                                const std::string fileName2)
{
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
    FIREnvelope env;
    try
    {
        env.initialize(5, RTSeis::ProcessingMode::POST_PROCESSING,
                       RTSeis::Precision::DOUBLE);
        env.transform(xSimple.size(), xSimple.data(), ySimple.data()); 
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("%s: Failed firEnv test 1", e.what());
        return EXIT_FAILURE;
    }
#ifdef __STDCPP_MATH_SPEC_FUNCS__ // __cplusplus > 201402L
    double tol = 1.e-13;
#else
    double tol = 1.e-5;
#endif
    double error;
    ippsNormDiff_L1_64f(ySimpleRef1.data(), ySimple.data(),
                        ySimple.size(), &error);
    if (error > tol)
    {
        RTSEIS_ERRMSG("Simple test 1 failed with error = %e", error);
        return EXIT_FAILURE;
    }
    try
    {
        env.initialize(6, RTSeis::ProcessingMode::POST_PROCESSING,
                       RTSeis::Precision::DOUBLE);
        env.transform(xSimple.size(), xSimple.data(), ySimple.data()); 
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("%s: Failed firEnv test 1", e.what());
        return EXIT_FAILURE;
    }
    ippsNormDiff_L1_64f(ySimpleRef2.data(), ySimple.data(), 
                        ySimple.size(), &error);
    if (error > tol)
    {
        RTSEIS_ERRMSG("Simple test 2 failed with error = %e", error);
        return EXIT_FAILURE;
    }
    // Do a more substantial post-processing test
    std::vector<double> x, upRef300, loRef300, upRef301, loRef301;
    readEnvelopeFile(fileName1, &x, &upRef300, &loRef300);
    readEnvelopeFile(fileName2, &x, &upRef301, &loRef301);
    if (upRef300.size() != 4000 || upRef301.size() != 4000)
    {
        RTSEIS_ERRMSG("%s", "Failed to load data file");
        return EXIT_FAILURE;
    }
    // Create with copy constructor
    FIREnvelope env300;
    try
    {
        env300.initialize(300,
                          RTSeis::ProcessingMode::POST_PROCESSING,
                          RTSeis::Precision::DOUBLE); 
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("%s; design n=300 failed", e.what());
        return EXIT_FAILURE;
    }
    env = env300;
    std::vector<double> up(x.size());
    env.transform(x.size(), x.data(), up.data());
    ippsNormDiff_L1_64f(upRef300.data(), up.data(),  upRef300.size(), &error);
    if (error/upRef300.size() > 1.e-8)
    {
        RTSEIS_ERRMSG("failed 300 with error = %e", error/upRef300.size());
        return EXIT_FAILURE;
    }
    // Test move construtor
    FIREnvelope env301;
    try  
    {
        env301.initialize(301,
                          RTSeis::ProcessingMode::POST_PROCESSING,
                          RTSeis::Precision::DOUBLE);
    }
    catch (const std::exception &e)
    {
        RTSEIS_ERRMSG("%s; design n=301 failed", e.what());
        return EXIT_FAILURE;
    }
    env = std::move(env301);
    env.transform(x.size(), x.data(), up.data());
    ippsNormDiff_L1_64f(upRef301.data(), up.data(),  upRef301.size(), &error);
    if (error/upRef301.size() > 1.e-8)
    {    
        RTSEIS_ERRMSG("failed 301 with error = %e", error/upRef301.size());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
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

int fft(const int nx, const std::complex<double> *x,
        const int ny, std::complex<double> *y)
{
    if (nx < 1 || ny < 1)
    {
        if (nx < 1){RTSEIS_ERRMSG("%s", "nx must be positive");}
        if (ny < 1){RTSEIS_ERRMSG("%s", "ny must be positive");}
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1; 
    }
    // Set space and make plan
    auto n = ny;
    size_t nbytes = sizeof(fftw_complex)*static_cast<size_t> (n);
    fftw_complex *in  = static_cast<fftw_complex *> (fftw_malloc(nbytes));
    memset(in, 0, nbytes);
    fftw_complex *out = reinterpret_cast<fftw_complex *> (y);
    fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    // Equal size transforms
    if (nx == ny) 
    {   
        std::memcpy(in, x, nbytes);
    }   
    // Truncate x to length of output array y
    else if (nx > ny) 
    {
        size_t ncopy = sizeof(fftw_complex)*static_cast<size_t> (ny);
        std::memcpy(in, x, ncopy);
    }   
    // Pad x to length of output array y
    else //if (nx < ny) 
    {
        size_t ncopy = sizeof(fftw_complex)*static_cast<size_t> (nx);
        std::memcpy(in, x, ncopy);
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
        if (nx < 1){RTSEIS_ERRMSG("%s", "nx must be positive");}
        if (ny < 1){RTSEIS_ERRMSG("%s", "ny must be positive");}
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
        return -1; 
    }
    // Set space and make plan
    auto n = ny; 
    size_t nbytes = sizeof(fftw_complex)*static_cast<size_t> (n);
    fftw_complex *in  = static_cast<fftw_complex *> (fftw_malloc(nbytes));
    memset(in, 0, nbytes);
    fftw_complex *out = reinterpret_cast<fftw_complex *> (y); //(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(size_t) n);
    fftw_plan p = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    // Equal size transforms
    if (nx == ny) 
    {
        std::memcpy(in, x, nbytes);
    }   
    // Truncate x to length of output array y
    else if (nx > ny) 
    {
        size_t ncopy = sizeof(fftw_complex)*static_cast<size_t> (ny);
        std::memcpy(in, x, ncopy);
    }
    // Pad x to length of output array y
    else //if (nx < ny) 
    {
        size_t ncopy = sizeof(fftw_complex)*static_cast<size_t> (nx);
        std::memcpy(in, x, ncopy);
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

int irfft(const int nx, const std::complex<double> x[],
          const int n, double y[])
{
    if (n < 1 || nx < 1)
    {
        if (n < 1){RTSEIS_ERRMSG("%s", "n must be positive");}
        if (nx < 1){RTSEIS_ERRMSG("%s", "nx must be positive");}
        return -1; 
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
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

int rfft(const int nx, double x[], const int n,
         const int ny, std::complex<double> y[])
{
    if (n < 1 || nx < 1)
    {
        if (n < 1){RTSEIS_ERRMSG("%s", "n must be positive");}
        if (nx < 1){RTSEIS_ERRMSG("%s", "nx must be positive");}
        return -1;
    }
    int ntf = n/2 + 1;
    if (ny < ntf)
    {
        RTSEIS_ERRMSG("ny = %d must be at least %d", ny, ntf);
        return -1;
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_ERRMSG("%s", "x is NULL");}
        if (y == nullptr){RTSEIS_ERRMSG("%s", "y is NULL");}
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
