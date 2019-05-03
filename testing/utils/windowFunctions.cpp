#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <exception>
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
int windowFunctions_sine_test(void);

int rtseis_test_utils_windowFunctions(void)
{
    int ierr;

    ierr = windowFunctions_sine_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed sine test");
        return EXIT_FAILURE;
    }

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
    std::vector<double> win20(20), win19(19); //double win20[20], win19[19];
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
    try
    {
        bartlett(20, win20.data());
    }
    catch (std::exception &e)
    {
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    try
    {
        bartlett(19, win19.data());
    }
    catch (std::exception &e)
    {   
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    std::vector<double> win20(20), win19(19); //double win20[20], win19[19];
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
    try
    {
        hann(20, win20.data());
    }
    catch (std::exception &e)
    {   
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    try
    {
        hann(19, win19.data());
    }
    catch (std::exception &e)
    {   
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    std::vector<double> win20(20), win19(19); //double win20[20], win19[19];
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
#ifdef __STDCPP_MATH_SPEC_FUNCS__ // __cplusplus > 201402L
    double tol = 1.e-14;
#else
    double tol = 1.e-6;
#endif
    try
    {
        kaiser(20, win20.data(), 2.5);
    }
    catch (std::exception &e)
    {
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    try
    {
        kaiser(19, win19.data(), 5.5);
    }
    catch (std::exception &e)
    {
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    std::vector<double> win20(20), win19(19); //double win20[20], win19[19];
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
    try
    {
        blackman(20, win20.data());
    }
    catch (std::exception &e)
    {
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    try
    {
        blackman(19, win19.data());
    }
    catch (std::exception &e)
    {
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    std::vector<double> win20(20), win19(19); //double win20[20], win19[19];
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
 
    try
    {
        hamming(20, win20.data());
    }
    catch (std::exception &e)
    {
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
    try
    {
        hamming(19, win19.data());
    }
    catch (std::exception &e)
    {
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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

int windowFunctions_sine_test(void)
{
    // N = 20
    // w20 = sin(linspace(0,N-1,N)*pi/(N-1))
    const double w20[20] = {0.000000000000000,
0.164594590280734,
0.324699469204683,
0.475947393037074,
0.614212712689668,
0.735723910673132,
0.837166478262529,
0.915773326655057,
0.969400265939330,
0.996584493006670,
0.996584493006670,
0.969400265939330,
0.915773326655057,
0.837166478262528,
0.735723910673132,
0.614212712689668,
0.475947393037074,
0.324699469204683,
0.164594590280734,
0.000000000000000};
    // N = 19 
    // w19 = sin(linspace(0,N-1,N)*pi/(N-1))
    const double w19[19] = {0.000000000000000,
0.173648177666930,
0.342020143325669,
0.500000000000000,
0.642787609686539,
0.766044443118978,
0.866025403784439,
0.939692620785908,
0.984807753012208,
1.000000000000000,
0.984807753012208,
0.939692620785908,
0.866025403784439,
0.766044443118978,
0.642787609686539,
0.500000000000000,
0.342020143325669,
0.173648177666930,
0.000000000000000};
/*
    const double w20[20] = {0.07845909572784494,
0.2334453638559054,
0.3826834323650898,
0.5224985647159488,
0.6494480483301837,
0.7604059656000309,
0.8526401643540922,
0.9238795325112867,
0.9723699203976766,
0.996917333733128,
0.996917333733128,
0.9723699203976767,
0.9238795325112867,
0.8526401643540923,
0.760405965600031,
0.6494480483301838,
0.5224985647159489,
0.3826834323650899,
0.23344536385590553,
0.07845909572784507};
    const double w19[19] = { 0.08257934547233232,
0.24548548714079912,
0.4016954246529694,
0.5469481581224268,
0.6772815716257411,
0.7891405093963936,
0.879473751206489,
0.9458172417006346,
0.9863613034027223,
1.0,
0.9863613034027224,
0.9458172417006346,
0.8794737512064891,
0.7891405093963936,
0.6772815716257411,
0.5469481581224268,
0.4016954246529698,
0.2454854871407995,
0.08257934547233267};
*/
    std::vector<double> win20(20), win19(19);
    try 
    {   
        sine(20, win20.data());
        sine(19, win19.data());
    }
    catch (std::exception &e) 
    {   
        RTSEIS_ERRMSG("Failed to compute bartlett20 %s", e.what());
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
