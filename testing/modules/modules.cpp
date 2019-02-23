#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#define RTSEIS_LOGGING 1
#include "modules.hpp"
#include "rtseis/log.h"

int main()
{
    int ierr;

    std::string dataDir = "data/";
    rtseis_utils_verbosity_setLoggingLevel(RTSEIS_SHOW_ALL);
    // Load data
    double *x = nullptr;
    int npts;
    ierr = readTextFile(&npts, &x, dataDir + "gse2.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed to load data");
        return EXIT_FAILURE;
    }

    ierr = rtseis_test_modules_demean();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed demean");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed demean");

    ierr = rtseis_test_modules_detrend();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed detrend");     
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed detrend");

    ierr = rtseis_test_modules_classicSTALTA(npts, x,
                                             dataDir + "classicSTALTA_ref.txt");
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed Classic STA/LTA module");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed classic STA/LTA");

    RTSEIS_INFOMSG("%s", "Passed all tets");
    free(x);
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
