#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include "rtseis/config.h"
#include "rtseis/log.h"
#include "utils.hpp"
#include <mkl.h>

int main()
{
    int ierr;
    rtseis_utils_verbosity_setLoggingLevel(RTSEIS_SHOW_ALL);

    ierr = rtseis_test_utils_windowFunctions();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed window test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed window test");

    ierr = rtseis_test_utils_normalization();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed normalization");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed normalization test");

    ierr = rtseis_test_utils_interpolation();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed interpolation");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed interpolation test");

    ierr = rtseis_test_utils_convolve();
    if (ierr != EXIT_SUCCESS)
    {   
        RTSEIS_ERRMSG("%s", "Failed convolve test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed convolve test");

    ierr = rtseis_test_utils_polynomial();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed polynomial tests");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed polynomial test");

    ierr = rtseis_test_utils_transforms();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed transform tests");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed transform test");
 
    ierr = rtseis_test_utils_design_iir_ap();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed analog prototype design");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed analog prototype design");

    ierr = rtseis_test_utils_design_iir();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed protype-based digital filter design");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed prototype-based digital filter design");

    ierr = rtseis_test_utils_design_zpk2sos();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed sos filter design");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed sos filter design");

    ierr = rtseis_test_utils_design_fir_fir1();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed FIR window design");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed FIR window design");

    ierr = rtseis_test_utils_design_freqs();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed freqs/freqz test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed freqs/freqz test");

    ierr = rtseis_test_utils_design_groupDelay();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed groupDelay test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed groupDelay test"); 

    ierr = rtseis_test_utils_filters();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed filters test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed filters test");

    // Last line
    mkl_free_buffers();
    return EXIT_SUCCESS;
}
