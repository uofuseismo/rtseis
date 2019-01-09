#include <stdio.h>
#include <stdlib.h>
#include <string>
#define RTSEIS_LOGGING 1
#include "utils.hpp"
#include "rtseis/utils/design.h"
#include "rtseis/log.h"

int rtseis_test_utils_design_iir_ap(void)
{
    int ierr;
    double k;
    AnalogPrototype ap;
    ZPK zpk, zpkRef;
    std::vector<std::complex<double>> pref;
    std::vector<std::complex<double>> zref;
    // Test butterworth order 1
    ierr = ap.butter(1);
    zpk = ap.getTransferFunction();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "order 1 failed");
        return EXIT_FAILURE;
    }
    k = 1;
    pref.push_back(std::complex<double> (-1, 0));
    zpkRef = ZPK(zref, pref, k);
    if (!(zpkRef == zpk))
    {
        RTSEIS_ERRMSG("%s", "Failed order 1 butter design");
        zpkRef.print();
        zpk.print();
        return EXIT_FAILURE;
    }
    zpkRef.clear();
    // Test butterworth order 4
    ierr = ap.butter(5);
    zpk = ap.getTransferFunction();
    if (ierr != 0)
    {   
        RTSEIS_ERRMSG("%s", "order 4 failed");
        return EXIT_FAILURE;
    }
    pref.clear();
    zref.clear();
    pref.resize(5);
    k = 1;
    pref[0] = std::complex<double> (-0.3090169943749474,  0.95105651629515353);
    pref[1] = std::complex<double> (-0.80901699437494745, 0.58778525229247303);
    pref[2] = std::complex<double> (-1, 0);
    pref[3] = std::complex<double> (-0.80901699437494745,-0.58778525229247303);
    pref[4] = std::complex<double> (-0.30901699437494751,-0.95105651629515353);
    zpkRef = ZPK(zref, pref, k);
    if (!(zpkRef == zpk))
    {   
        RTSEIS_ERRMSG("%s", "Failed order 5 butter design");
        zpkRef.print();
        zpk.print();
        return EXIT_FAILURE;
    }
    zpkRef.clear();
    // Test order 1 cheby1
    ierr = ap.cheb1ap(1, 2.2);
    zpk = ap.getTransferFunction();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "order 1 failed");
        return EXIT_FAILURE;
    }
    zref.clear();
    pref.clear();
    pref.resize(1);
    k = 1.2313003041963828;
    pref[0] = std::complex<double> (-1.2313003041963828, 0);
    zpkRef = ZPK(zref, pref, k);
    if (!(zpkRef == zpk))
    {   
        RTSEIS_ERRMSG("%s", "Failed order 1 cheb1 design");
        zpkRef.print();
        zpk.print();
        return EXIT_FAILURE;
    }
    // Test order 6 cheby1
    ierr = ap.cheb1ap(6, 0.994);
    zpk = ap.getTransferFunction();
    if (ierr != 0)
    {   
        RTSEIS_ERRMSG("%s", "order 1 failed");
        return EXIT_FAILURE;
    }   
    zref.clear();
    pref.clear();
    pref.resize(6);
    k = 0.061620501119488615;
    pref[0] = std::complex<double> (-0.062314231644038744,+0.9935274525619241); 
    pref[1] = std::complex<double> (-0.17024564688613014,+0.72731257398980598);
    pref[2] = std::complex<double> (-0.23255987853016891,+0.26621487857211812);
    pref[3] = std::complex<double> (-0.23255987853016891,-0.26621487857211795);
    pref[4] = std::complex<double> (-0.17024564688613017,-0.72731257398980587);
    pref[5] = std::complex<double> (-0.062314231644038813,-0.99352745256192398);
    zpkRef = ZPK(zref, pref, k); 
    if (!(zpkRef == zpk))
    {   
        RTSEIS_ERRMSG("%s", "Failed order 1 cheb1 design");
        zpkRef.print();
        zpk.print();
        return EXIT_FAILURE;
    } 
    // Test order 2 cheby2 
    ierr = ap.cheb2ap(1, 1.1);
    zpk = ap.getTransferFunction();
    if (ierr != 0)
    {   
        RTSEIS_ERRMSG("%s", "order 1 failed");
        return EXIT_FAILURE;
    }
    zref.clear();
    pref.clear();
    pref.resize(1);
    k = 1.862583192806328;
    pref[0] = std::complex<double> (-1.862583192806328,0);
    zpkRef = ZPK(zref, pref, k);
    if (!(zpkRef == zpk))
    {
        RTSEIS_ERRMSG("%s", "Failed order 1 cheb2 design");
        zpkRef.print();
        zpk.print();
        return EXIT_FAILURE;
    }
    // Test order 6 cheby2 
    ierr = ap.cheb2ap(6, 1.2);
    zpk = ap.getTransferFunction();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "order 1 failed");
        return EXIT_FAILURE;
    }
    k = 0.8709635899560811;
    zref.resize(6);
    pref.resize(6);
    zref[0] = std::complex<double> (0, -1.035276180410083);
    zref[1] = std::complex<double> (0, -1.4142135623730951);
    zref[2] = std::complex<double> (0, -3.8637033051562737);
    zref[3] = std::complex<double> (0,  3.8637033051562737);
    zref[4] = std::complex<double> (0,  1.4142135623730951);
    zref[5] = std::complex<double> (0,  1.035276180410083); 
    pref[0] = std::complex<double> (-0.024686186266327684,-1.0305393933278832);
    pref[1] = std::complex<double> (-0.12492582633346083,-1.3973824335027194);
    pref[2] = std::complex<double> (-1.1553327165440157,-3.462761441343769);
    pref[3] = std::complex<double> (-1.1553327165440157,+3.462761441343769);
    pref[4] = std::complex<double> (-0.12492582633346083,+1.3973824335027194);
    pref[5] = std::complex<double> (-0.024686186266327684,+1.0305393933278832);
    zpkRef = ZPK(zref, pref, k);
    if (!(zpkRef == zpk))
    {
        RTSEIS_ERRMSG("%s", "Failed order 6 cheb2 design");
        zpkRef.print();
        zpk.print();
    }
    return EXIT_SUCCESS;
}
