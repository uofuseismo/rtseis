#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/modules/demean.hpp"
#include "rtseis/modules/detrend.hpp"
#include "rtseis/log.h"
#include "modules.hpp"

using namespace RTSeis::Modules;


int rtseis_test_modules_detrend(void)
{
    DetrendParameters parms;
    Detrend detrend(parms);
    int npts = 360001;
    for (int j=0; j<2; j++)
    {
        npts = npts + j;
        double *x = new double[npts];
        double *y = new double[npts];
        // Detrend should remove a line
        for (int i=0; i<npts; i++)
        {
            x[i] = 1.1 + 0.3*static_cast<double> (i);
        }
        try 
        {
            detrend.apply(npts, x, y);
        }
        catch (std::exception(&e))
        {
            fprintf(stderr, "%s", e.what());
            return EXIT_FAILURE;
        }
        double dmax = 0;
        for (int i=0; i<npts; i++)
        {
            dmax = std::max(dmax, std::abs(y[i]));
        }
        // For this many numbers this is about numerical precision
        if (dmax > 1.e-6)
        {
            RTSEIS_ERRMSG("%s", "Failed to compute detrend");
            return -1;
        }
        delete[] x;
        delete[] y;
    }
    return EXIT_SUCCESS;
}

int rtseis_test_modules_demean(void)
{
    DemeanParameters parms;
    Demean demean(parms);
    int npts = 3601;
    for (int j=0; j<2; j++)
    {
        npts = npts + j;
        double *x = new double[npts];
        double *y = new double[npts];
        // Detrend should remove a line
        for (int i=0; i<npts; i++)
        {
            x[i] = 1.1 + static_cast<double> (j);
        }
        try
        {
            demean.apply(npts, x, y); 
        }
        catch (std::exception(&e))
        {
            fprintf(stderr, "%s", e.what());
            return EXIT_FAILURE;
        }
        double dmax = 0;
        for (int i=0; i<npts; i++)
        {
            dmax = std::max(dmax, std::abs(y[i]));
        }
        // For this many numbers this is about numerical precision
        if (dmax > 1.e-10)
        {
            RTSEIS_ERRMSG("%s", "Failed to compute demean");
            return -1;
        }
        delete[] x;
        delete[] y;
    }
    return EXIT_SUCCESS;
}

