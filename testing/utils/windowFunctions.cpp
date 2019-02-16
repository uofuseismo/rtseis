#include <cstdio>
#include <cstdlib>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "utils.hpp"
#include "rtseis/utilities/windowFunctions.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::WindowFunctions;

int windowFunctions_bartlett_test(void);
int windowFunctions_blackman_test(void);
int windowFunctions_hann_test(void);
int windowFunctions_kaiser_test(void);
int windowFunctions_hamming_test(void);

int rtseis_test_utils_windowFunctions(void)
{
    int ierr;

    ierr = windowFunctions_hamming_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed hamming test");
        return EXIT_FAILURE;
    }

    ierr = windowFunctions_blackman_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed blackman test");
        return EXIT_FAILURE;
    }

    ierr = windowFunctions_hann_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed hann test");
        return EXIT_FAILURE;
    }

    ierr = windowFunctions_bartlett_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed bartlett test");
        return EXIT_FAILURE;
    }

    ierr = windowFunctions_kaiser_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed kaiser test");
        return EXIT_FAILURE;
    }
    RTSEIS_INFOMSG("%s", "Passed window tests"); 
    return EXIT_SUCCESS;
}

int windowFunctions_bartlett_test(void)
{
    double win20[20], win19[19];
    const double w20[20] = {
                   0, 0.105263157894737, 0.210526315789474,
                   0.315789473684211,  0.421052631578947,   0.526315789473684,
                   0.631578947368421,  0.736842105263158,  0.842105263157895,
                   0.947368421052632, 0.947368421052632,  0.842105263157895,
                   0.736842105263158,  0.631578947368421,   0.526315789473684,
                   0.421052631578947, 0.315789473684211, 0.210526315789474,
                   0.105263157894737,
                   0};
    const double w19[19] ={0, 0.111111111111111, 0.222222222222222,
                 0.333333333333333, 0.444444444444444, 0.555555555555556,
                0.666666666666667,   0.777777777777778,   0.888888888888889,
                1.000000000000000, 0.888888888888889,   0.777777777777778,
                0.666666666666667, 0.555555555555556,   0.444444444444444,
                0.333333333333333,  0.222222222222222, 0.111111111111111,
                   0};
    int ierr = bartlett(20, win20);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute bartlett20");
        return EXIT_FAILURE;
    }
    for (int i=0; i<20; i++)
    {
        if (std::abs(w20[i] - win20[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Error in w20 %lf %lf", w20[i], win20[i]);
            return EXIT_FAILURE;
        }
    }
    ierr = bartlett(19, win19);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute bartlett19");
        return EXIT_FAILURE;
    }
    for (int i=0; i<19; i++)
    {
        if (std::abs(w19[i] - win19[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Error in w19 %lf %lf", w19[i], win19[i]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int windowFunctions_hann_test(void)
{
    double win20[20], win19[19];
    const double w20[20] = {0,
   0.027091379149683,
   0.105429745301803,
   0.226525920938787,
   0.377257256429600,
   0.541289672736166,
   0.700847712326485,
   0.838640785812870,
   0.939736875603244,
   0.993180651701361,
   0.993180651701361,
   0.939736875603244,
   0.838640785812870,
   0.700847712326485,
   0.541289672736166,
   0.377257256429600,
   0.226525920938787,
   0.105429745301803,
   0.027091379149683,
                   0};
    const double w19[19] ={0,
   0.030153689607046,
   0.116977778440511,
   0.250000000000000,
   0.413175911166535,
   0.586824088833465,
   0.750000000000000,
   0.883022221559489,
   0.969846310392954,
   1.000000000000000,
   0.969846310392954,
   0.883022221559489,
   0.750000000000000,
   0.586824088833465,
   0.413175911166535,
   0.250000000000000,
   0.116977778440511,
   0.030153689607046,
                   0
                   };
    int ierr = hann(20, win20);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hann20");
        return EXIT_FAILURE;
    }
    for (int i=0; i<20; i++)
    {
        if (std::abs(w20[i] - win20[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Error in w20 %lf %lf", w20[i], win20[i]);
            return EXIT_FAILURE;
        }
    }
    ierr = hann(19, win19);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hann19");
        return EXIT_FAILURE;
    }
    for (int i=0; i<19; i++)
    {
        if (std::abs(w19[i] - win19[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Error in w19 %lf %lf", w19[i], win19[i]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int windowFunctions_kaiser_test(void)
{
    double win20[20], win19[19];
    // beta = 2.5
    const double k20[20] = {0.303966229415369, 0.406333120777527,
                            0.511011861080447, 0.614155742009957,
                            0.711825746588635, 0.800182020135053,
                            0.875673636668776, 0.935216242888425,
                            0.976347754218060, 0.997353444300571,
                            0.997353444300571, 0.976347754218060,
                            0.935216242888425, 0.875673636668776,
                            0.800182020135053, 0.711825746588635,
                            0.614155742009957, 0.511011861080447,
                            0.406333120777527, 0.303966229415369};
    // beta = 5.5
    const double k19[19] = {0.023422141030647, 0.078225306346850,
                            0.166678217276960, 0.288508741610676,
                            0.436709600995736, 0.597737293640817,
                            0.753274215094738, 0.883279467813599,
                            0.969699049461490, 1.000000000000000,
                            0.969699049461490, 0.883279467813599,
                            0.753274215094738, 0.597737293640817,
                            0.436709600995736, 0.288508741610676,
                            0.166678217276960, 0.078225306346850,
                            0.023422141030647};
    int ierr;
#if __cplusplus > 201402L
    double tol = 1.e-14;
#else
    double tol = 1.e-6;
#endif
    ierr = kaiser(20, win20, 2.5);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute kaiser20");
        return EXIT_FAILURE;
    }
    for (int i=0; i<20; i++)
    {
        if (std::abs(k20[i] - win20[i]) > tol)
        {
            RTSEIS_ERRMSG("Error in k20 %lf %lf", k20[i], win20[i]);
            return EXIT_FAILURE;
        }
    }
    ierr = kaiser(19, win19, 5.5);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute kaiser19");
        return EXIT_FAILURE;
    }
    for (int i=0; i<19; i++)
    {
        if (std::abs(k19[i] - win19[i]) > tol)
        {
            RTSEIS_ERRMSG("Error in k19 %lf %lf", k19[i], win19[i]);
            return EXIT_FAILURE;
        }   
    }
    return EXIT_SUCCESS; 
}

int windowFunctions_blackman_test(void)
{
    double win20[20], win19[19];
    const double w20[20] = { 0,
   0.010222619901394,
   0.045068584273067,
   0.114390286966549,
   0.226899356333081,
   0.382380768463948,
   0.566665186596425,
   0.752034438175084,
   0.903492728253039,
   0.988846031037412,
   0.988846031037412,
   0.903492728253039,
   0.752034438175084,
   0.566665186596425,
   0.382380768463948,
   0.226899356333081,
   0.114390286966549,
   0.045068584273067,
   0.010222619901394,
                   0};
    const double w19[19] = {0,
   0.011437245056564,
   0.050869632653865,
   0.130000000000000,
   0.258000501503662,
   0.431648679170593,
   0.630000000000000,
   0.816914075772843,
   0.951129865842472,
   1.000000000000000,
   0.951129865842472,
   0.816914075772843,
   0.630000000000000,
   0.431648679170593,
   0.258000501503662,
   0.130000000000000,
   0.050869632653865,
   0.011437245056564,
                   0};
    int ierr = blackman(20, win20);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute blackman20");
        return EXIT_FAILURE;
    }
    for (int i=0; i<20; i++)
    {
        if (std::abs(w20[i] - win20[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Error in w20 %lf %lf", w20[i], win20[i]);
            return EXIT_FAILURE;
        }
    }
    ierr = blackman(19, win19);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute blackman19");
        return EXIT_FAILURE;
    }
    for (int i=0; i<19; i++)
    {
        if (std::abs(w19[i] - win19[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Error in w19 %lf %lf", w19[i], win19[i]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
};

int windowFunctions_hamming_test(void)
{
    double win20[20], win19[19];
    const double w20[20] = {0.080000000000000,
   0.104924068817708,
   0.176995365677659,
   0.288403847263684,
   0.427076675915232,
   0.577986498917273,
   0.724779895340366,
   0.851549522947841,
   0.944557925554985,
   0.993726199565252,
   0.993726199565252,
   0.944557925554985,
   0.851549522947841,
   0.724779895340366,
   0.577986498917273,
   0.427076675915232,
   0.288403847263684,
   0.176995365677659,
   0.104924068817708,
   0.080000000000000 };
    const double w19[19] = {0.080000000000000,
   0.107741394438482,
   0.187619556165270,
   0.310000000000000,
   0.460121838273212,
   0.619878161726788,
   0.770000000000000,
   0.892380443834730,
   0.972258605561518,
   1.000000000000000,
   0.972258605561518,
   0.892380443834730,
   0.770000000000000,
   0.619878161726788,
   0.460121838273212,
   0.310000000000000,
   0.187619556165270,
   0.107741394438482,
   0.080000000000000};
 
    int ierr = hamming(20, win20);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hamming20");
        return EXIT_FAILURE;
    }
    for (int i=0; i<20; i++)
    {
        if (std::abs(w20[i] - win20[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Error in w20 %lf %lf", w20[i], win20[i]);
            return EXIT_FAILURE;
        }
    }
    ierr = hamming(19, win19);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hamming19");
        return EXIT_FAILURE;
    }
    for (int i=0; i<19; i++)
    {
        if (std::abs(w19[i] - win19[i]) > 1.e-14)
        {
            RTSEIS_ERRMSG("Error in w19 %lf %lf", w19[i], win19[i]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
