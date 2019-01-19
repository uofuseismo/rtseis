#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <chrono>
#define RTSEIS_LOGGING 1
#include "rtseis/utils/filters.hpp"
#include "rtseis/log.h"
#include "utils.hpp"

using namespace RTSeis::Utils::Filters;
static int readTextFile(int *npts, double *xPtr[],
                        const std::string fileName = "utils/data/gse2.txt");
static int filters_downsample_test(const int npts, const double x[]);
static int filters_medianFilter_test(const int npts, const double x[],
                                     const std::string fileName);

int rtseis_test_utils_filters(void)
{
    double dt = 1.0/200.0;
    int npts;
    double *x = nullptr;
    std::string dataDir = "utils/data/";
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
    // Apply the median filter
    ierr = filters_medianFilter_test(npts, x,
                                     dataDir + "medianFilterReference.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed median filter test");
        return EXIT_FAILURE;
    }

    if (x != nullptr){free(x);}
    return EXIT_SUCCESS;
}
int filters_medianFilter_test(const int npts, const double x[], const std::string fileName)
{
    fprintf(stdout, "Testing median filter...\n");
    double xin[8] = {1, 2, 127, 4, 5, 0, 7, 8};
    double y8[8];
    double yref3[8] = {1, 2, 4, 5, 4, 5, 7, 7}; // Matlab soln; IPP has edge effect
    double yref5[8] = {1, 2, 4, 4, 5, 5, 5, 0};
    int ierr;
    MedianFilter median;
    bool lrt = false;
    ierr = median.initialize(3, lrt, RTSEIS_DOUBLE);
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
    ierr = median.initialize(5, lrt, RTSEIS_DOUBLE);
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
    median.initialize(11, lrt, RTSEIS_DOUBLE);
    auto timeStart = std::chrono::high_resolution_clock::now();
    double *y = new double[npts];
    ierr = median.apply(npts, x, y);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute reference solution");
        return EXIT_FAILURE;
    }
    for (int i=0; i<npts; i++)
    {
        if (std::abs(y[i] - yref[i]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed to compute reference soln %d %lf %lf",
                          i, y[i], yref[i]);
            return EXIT_FAILURE;
        }
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tdif = timeEnd - timeStart;
    fprintf(stdout, "Reference solution computation time %.8lf (s)\n",
            tdif.count());
    // Now do the packetized tests
    lrt = true;
    median.initialize(11, lrt, RTSEIS_DOUBLE);
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
                ierr = median.apply(nptsPass, &x[nxloc], &y[nxloc]);
                if (ierr != 0)
                {
                    RTSEIS_ERRMSG("%s", "Failed to apply median filter");
                    return EXIT_FAILURE;
                }
                nxloc = nxloc + nptsPass;
            }
            median.resetInitialConditions();
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
    const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE;
    // Call this in post-processing for a couple different decimation rates
    bool lrt = false;
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
        lrt = false;
        ierr = downsample.initialize(iq, lrt, precision); 
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
        lrt = true; 
        ierr = downsample.initialize(iq, lrt, precision);
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
                                        npts+1-nxloc, &nyDec, &y[nyloc]);
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
                RTSEIS_ERRMSG("%s", "Failed fixed packet size test");
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
                                    npts+1-nxloc, &nyDec, &y[nyloc]);
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
            RTSEIS_ERRMSG("%s", "Failed fixed packet size test");
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
