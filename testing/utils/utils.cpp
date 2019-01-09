#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include "rtseis/config.h"
#include "rtseis/log.h"
#include "utils.hpp"

int main()
{
    int ierr;
    rtseis_utils_verbosity_setLoggingLevel(RTSEIS_SHOW_ALL);

    ierr = rtseis_test_utils_polynomial();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed polynomial tests");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed polynomial test");
 
    ierr = rtseis_test_utils_design_iir_ap();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed analog prototype design");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed analog prototype design");

    return EXIT_SUCCESS;
}
