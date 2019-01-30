#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <chrono>
#define RTSEIS_LOGGING 1
#include "rtseis/modules/classicSTALTA.hpp"
#include "rtseis/log.h"
#include "modules.hpp"

using namespace RTSeis::Modules;


int rtseis_test_modules_classicSTALTA(const int npts, const double x[],
                                      const std::string fileName)
{
    fprintf(stdout, "Testing classic STA/LTA...\n");
    srand(40234);
    // Load the data file
    double *yref = nullptr;
    int nptsIn;
    int ierr = readTextFile(&nptsIn, &yref, fileName);
    if (npts != nptsIn)
    {
        RTSEIS_ERRMSG("%s", "Data load failed");
        return EXIT_FAILURE;
    }
    // Set the post-processing and real-time parameters 
    double dt = 1.0/200;
    int nlta = static_cast<int> (10/dt);
    int nsta = static_cast<int> (5/dt);
    ClassicSTALTAParameters ppParms(nsta, nlta,
                                    RTSeis::ProcessingMode::POST_PROCESSING,
                                    RTSeis::Precision::DOUBLE);
    ClassicSTALTAParameters rtParms = ppParms;
    rtParms.setProcessingMode(RTSeis::ProcessingMode::REAL_TIME);
    rtParms.setChunkSize(512);
    // Initialize the stalta
    ClassicSTALTA stalta(ppParms);
    if (!stalta.isInitialized())
    {
        RTSEIS_ERRMSG("%s", "Failed to initialize stalta");
        return EXIT_FAILURE;
    }
    double *y = new double[npts];
    auto timeStart = std::chrono::high_resolution_clock::now();
    ierr = stalta.apply(npts, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute STA/LTA");
        return EXIT_FAILURE;
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    double error = 0;
    for (int i=0; i<npts; i++)
    {
        error = std::max(error, std::abs(y[i] - yref[i]));
    }
    if (error > 1.e-10)
    {
        RTSEIS_ERRMSG("Failed post-processing test w/ error %lf", error);
        return EXIT_FAILURE;
    }
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Passed post-processing test in %.8e (s)\n", tdif.count());
    // Do real-time test
    stalta = ClassicSTALTA(rtParms);
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
                ierr = stalta.apply(nptsPass, &x[nxloc], &y[nxloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply median filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
            }
            stalta.resetInitialConditions();
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
                        "Passed STA/LTA fixed packet size %4d in %.8e (s)\n",
                        packetSize[ip], tdif.count());
            }
            else
            {
                fprintf(stdout,
                        "Passed STA/LTA random in %.8e (s)\n", tdif.count());
            }
        }
    }

    free(yref);
    delete[] y;
    return EXIT_SUCCESS;
}

