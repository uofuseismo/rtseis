#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define RTSEIS_LOGGING 1
#include "modules.hpp"
#include "rtseis/log.h"

int main()
{
    int ierr;

    rtseis_utils_verbosity_setLoggingLevel(RTSEIS_SHOW_ALL);
    ierr = rtseis_test_modules_oneBitNormalization();
 
    RTSEIS_INFOMSG("%s", "Passed one-bit normalization test");

    ierr = rtseis_test_modules_detrend();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed detrend");     
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed detrend");


    RTSEIS_INFOMSG("%s", "Passed all tets");
    return EXIT_SUCCESS;
}
