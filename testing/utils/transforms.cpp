#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <float.h>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <complex>
#include <vector>
#define RTSEIS_LOGGING 1
#include "rtseis/utils/transforms.hpp"
#include "rtseis/log.h"
#include "utils.hpp"
#include <ipps.h>
#include <fftw/fftw3.h>

using namespace RTSeis::Utils::Transforms;

int transforms_nextPow2_test(void);
int transforms_phase_test(void);
int transforms_unwrap_test(void);
int transforms_test_dft(void);
int rfft(const int nx, double x[], const int n,
         const int ny, std::complex<double> y[]);
int irfft(const int nx, const std::complex<double> x[],
          const int n, double y[]);

int rtseis_test_utils_transforms(void)
{
    int ierr = transforms_nextPow2_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed nextpow2 test");
        return EXIT_FAILURE;
    }
    ierr = transforms_phase_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed phase test");
        return EXIT_FAILURE;
    }
    ierr = transforms_unwrap_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to unwrap phase");
        return EXIT_FAILURE;
    }
    ierr = transforms_test_dft();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute rdft");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int transforms_nextPow2_test(void)
{
    if (DFTUtils::nextPow2(0) != 1)
    {   
        RTSEIS_ERRMSG("%s", "Failed 0 test");
        return EXIT_FAILURE;
    }   
    if (DFTUtils::nextPow2(1) != 1)
    {
        RTSEIS_ERRMSG("%s", "Failed 1 test");
        return EXIT_FAILURE;
    }
    if (DFTUtils::nextPow2(2) != 2)
    {
        RTSEIS_ERRMSG("%s", "Failed 2 test");
        return EXIT_FAILURE;
    }
    if (DFTUtils::nextPow2(3) != 4)
    {
        RTSEIS_ERRMSG("%s", "Failed 3 test");
        return EXIT_FAILURE;
    }
    if (DFTUtils::nextPow2(4) != 4)
    {
        RTSEIS_ERRMSG("%s", "Failed 4 test");
        return EXIT_FAILURE;
    }
    if (DFTUtils::nextPow2(5) != 8)
    {
        RTSEIS_ERRMSG("%s", "Failed 8 test");
        return EXIT_FAILURE;
    }
    if (DFTUtils::nextPow2(1200) != 2048)
    {
        RTSEIS_ERRMSG("%s", "Failed 1200 test");
        return EXIT_FAILURE;
    }
    if (DFTUtils::nextPow2(120000) != 131072)
    {
        RTSEIS_ERRMSG("%s", "Failed 120000 test");
        return EXIT_FAILURE;
    }
    if (DFTUtils::nextPow2(131072) != 131072)
    {
        RTSEIS_ERRMSG("%s", "Failed 131072 test");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int transforms_unwrap_test(void)
{
    int n = 23;
    double p[23] = {0, -1.5728, -1.5747, -1.5772, -1.5790,
                    -1.5816, -1.5852, -1.5877, -1.5922,
                    -1.5976, -1.6044, -1.6129, -1.6269,
                    -1.6512, -1.6998, -1.8621,  1.7252,
                     1.6124,  1.5930,  1.5916,  1.5708,
                     1.5708,  1.5708};
    double qref[23] = {0.00000,  -1.57280,  -1.57470,  -1.57720,
                       -1.57900,  -1.58160,  -1.58520,  -1.58770,
                       -1.59220,  -1.59760,  -1.60440,  -1.61290,
                       -1.62690,  -1.65120,  -1.69980,  -1.86210,
                       -4.557985307179586,  -4.670785307179586,
                       -4.690185307179586,  -4.691585307179587,
                       -4.712385307179586,  -4.712385307179586,
                       -4.712385307179586};
    double qref2[23] = {1.000000000000000,  -3.710385307179586,
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
    double q[23];
    int ierr = DFTUtils::unwrap(n, p, q);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to called unwrap");
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
    ierr = DFTUtils::unwrap(n, p, q, tol);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to called unwrap");
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
    const double tr[7] = {-0.785398163397448, 0.463647609000806,
                          -0.3217505543966, 0.244978663126864,
                           0, 1.570796326794897, 0};
    std::complex<double> z[7];
    z[0] = std::complex<double> (1, -1);
    z[1] = std::complex<double> (2, +1);
    z[2] = std::complex<double> (3, -1);
    z[3] = std::complex<double> (4, +1);
    z[4] = std::complex<double> (0, +0);
    z[5] = std::complex<double> (0, +1);
    z[6] = std::complex<double> (1, +0);
    double angle[7];
    int ierr = DFTUtils::phase(n, z, angle);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute phase");
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
    ierr = DFTUtils::phase(n, z, angle, lwantDeg);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute phase");
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

int transforms_test_dft(void)
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
            RTSEIS_ERRMSG("%s", "Failed to compute rfft");
            return EXIT_FAILURE;
        }
        // Compute an FFT w/ FFTw
        int np2 = DFTUtils::nextPow2(npts);
        int lenfft = np2/2 + 1;
        std::complex<double> *zrefFFT = new std::complex<double>[lenfft];
        ierr = rfft(npts, x, np2, lenfft, zrefFFT); 

        // Initialize the DFT
        DFTR2C dft; 
        bool ldoFFT = false;
        ierr = dft.initialize(npts, ldoFFT, RTSeis::Precision::DOUBLE); 
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to initialize dft");
            return EXIT_FAILURE;
        }
        if (dft.getTransformLength() != lendft)
        {
            RTSEIS_ERRMSG("%s", "Inconsistent sizes");
            return EXIT_FAILURE;
        }
        std::complex<double> *z = new std::complex<double>[lendft];
        ierr = dft.forwardTransform(npts, x, lendft, z);
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to forward transform");
            return EXIT_FAILURE;
        }
        double error = 0; 
        for (int i=0; i<lendft; i++)
        {
            error = error + std::abs(z[i] - zrefDFT[i]);
        }
        if (error > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute dft %.10e", error);
            return EXIT_FAILURE;
        }
        // Stress test it
        if (j == 1)
        {
            auto timeStart = std::chrono::high_resolution_clock::now();
            for (int kiter=0; kiter<niter; kiter++)
            {
                ierr = dft.forwardTransform(npts, x, lendft, z); 
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed DFT");
                    return EXIT_FAILURE;
                }
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            error = 0;
            for (int i=0; i<lendft; i++)
            {
                error = error + std::abs(z[i] - zrefDFT[i]);
            }
            if (error > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed to compute dft %.10e", error);
                return EXIT_FAILURE;
            }
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            fprintf(stdout, "Average DFT time %.8lf (s)\n",
                    tdif.count()/static_cast<double>(niter)); 
 
        }
        delete[] z;
        // Redo this for an FFT
        ldoFFT = true;
        dft.initialize(npts, ldoFFT, RTSeis::Precision::DOUBLE);
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
        ierr = dft.forwardTransform(npts, x, lenfft, z); 
        if (ierr != 0)
        {
            RTSEIS_ERRMSG("%s", "Failed to forward transform");
            return EXIT_FAILURE;
        }
        error = 0;  
        for (int i=0; i<lenfft; i++)
        {
            error = error + std::abs(z[i] - zrefFFT[i]);
        }
        if (error > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute fft %.10e", error);
            return EXIT_FAILURE;
        }
        // Stress test it
        if (j == 1)
        {
            auto timeStart = std::chrono::high_resolution_clock::now();
            for (int kiter=0; kiter<niter; kiter++)
            {
                ierr = dft.forwardTransform(npts, x, lenfft, z);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed FFT");
                    return EXIT_FAILURE;
                }
            }
            auto timeEnd = std::chrono::high_resolution_clock::now();
            error = 0;  
            for (int i=0; i<lenfft; i++)
            {
                error = error + std::abs(z[i] - zrefFFT[i]);
            }
            if (error > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed to compute fft %.10e", error);
                return EXIT_FAILURE;
            }
            std::chrono::duration<double> tdif = timeEnd - timeStart;
            fprintf(stdout, "Average FFT time %.8lf (s)\n",
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
