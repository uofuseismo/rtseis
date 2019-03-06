#include <cstdio>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/postProcessing/singleChannel/waveform.hpp"


using namespace RTSeis::PostProcessing;

int testDemean(void);
int testDetrend(void);

int main(void)
{
    rtseis_utils_verbosity_setLoggingLevel(RTSEIS_SHOW_ALL);
    int ierr;
    ierr = testDemean();
    if (ierr != EXIT_SUCCESS)
    { 
        RTSEIS_ERRMSG("%s", "Failed demean test");
        return EXIT_FAILURE;
    } 
    RTSEIS_INFOMSG("%s", "Passed demean test");

    ierr = testDetrend();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed detrend test");
        return EXIT_FAILURE;
    } 
    RTSEIS_INFOMSG("%s", "Passed detrend test");
    return EXIT_SUCCESS; 
}


int testDemean(void)
{
    // Do two passes
    for (int j=0; j<2; j++)
    {
    // Create a reference signal from which to remove the mean
    std::vector<double> x;
    int npts = 5000 + j;
    x.resize(npts);
    double xmean = 5.1 + static_cast<double> (j);
    std::fill(x.begin(), x.end(), xmean);
    std::vector<double> y;
//! [ppSCDemeanExample]
    RTSeis::PostProcessing::SingleChannel::Waveform waveform;
    try
    {
        waveform.setData(x);    // Set waveform data to demean
        waveform.demean();      // Demean the data
        waveform.getData(y);    // Get the demeaned data in a vector, y
    }
    catch (std::invalid_argument &ia)
    {
        fprintf(stderr, "Demean failed %s\n", ia.what()); 
        return EXIT_FAILURE;
    }
//! [ppSCDemeanExample]
    // Verify
    if (y.size() != x.size())
    {
        RTSEIS_ERRMSG("%s", "Inconsistent sizes");
        return EXIT_FAILURE;
    }
    double maxAbs = 0;
    ippsMaxAbs_64f(y.data(), npts, &maxAbs);
    if (maxAbs > 1.e-13)
    {
        RTSEIS_ERRMSG("Demean failed %e", maxAbs);
    }
    } // Loop on passes

    return EXIT_SUCCESS;
}
//============================================================================//
int testDetrend(void)
{
    for (int j=0; j<2; j++)
    {
    std::vector<double> x;
    std::vector<double> y;
    int npts = 50001 + j;
    x.resize(npts);
    #pragma omp simd
    for (int i=0; i<npts; i++){x[i] = 1.1 + 0.3*static_cast<double> (i);}
//! [ppSCDetrendExample]
    RTSeis::PostProcessing::SingleChannel::Waveform waveform;
    try 
    {
        waveform.setData(x);    // Set the data to detrend
        waveform.detrend();     // Detrend
        waveform.getData(y);    // Get the detrended data in a vector, y
    }   
    catch (std::invalid_argument &ia)
    {
        fprintf(stderr, "Detrend failed %s\n", ia.what()); 
        return EXIT_FAILURE;
    }   
//! [ppSCDetrendExample]
    double maxAbs = 0;
    ippsMaxAbs_64f(y.data(), npts, &maxAbs);
    // ~50,000 points loses about 5 digits which is sensible 
    if (maxAbs > 1.e-9)
    {
        RTSEIS_ERRMSG("Demean failed %e", maxAbs);
        return EXIT_FAILURE;
    }
    }; // Loop on npts
    return EXIT_SUCCESS;
}
