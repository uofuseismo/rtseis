#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <chrono>
#define RTSEIS_LOGGING 1
#include "utils.hpp"
#include "rtseis/utilities/normalization/signBit.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::Normalization;

static int test_normalization_signBit(void);

int rtseis_test_utils_normalization(void)
{
    int ierr;
    ierr = test_normalization_signBit();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed sign bit normalization");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed signbit");
    return EXIT_SUCCESS;
}

int test_normalization_signBit(void)
{
    fprintf(stdout, "%s: Testing signBit...\n", __func__);
    SignBit signBit;
    signBit.initialize();
    // Create a random signal
    int npts = 6500;
    double *x = new double[npts];
    for (int i=0; i<npts; i++)
    {   
        x[i] = (static_cast<double> (rand())/RAND_MAX - 0.5)*100;
    }   
    double *y = new double[npts];
    double *yref = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    try
    {
        signBit.apply(npts, x, y); 
    }
    catch (std::exception &e)
    {
        RTSEIS_ERRMSG("%s", e.what());
        return EXIT_FAILURE;
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    for (int i=0; i<npts; i++)
    {
        double sign =-1;
        if (x[i] >= +0){sign =+1;}
        yref[i] = sign;
        if (std::abs(y[i] - sign) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute sign %lf %lf %lf",
                          x[i], sign, y[i]);
            return EXIT_FAILURE;
        }   
    }   
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    SignBit signBitRT = signBit; 
    // Do a real-time test
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
                try
                {
                    signBitRT.apply(nptsPass, &x[nxloc], &y[nxloc]);
                }
                catch (std::exception &e)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply onebit filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
            }
            signBitRT.resetInitialConditions();
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
                        "Passed signBit filter fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed onebit filter random in %.8e (s)\n", tdif.count());
            }
        }
    }
    delete[] x;
    delete[] y;
    delete[] yref;
    return EXIT_SUCCESS;
}
