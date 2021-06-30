#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "private/throw.hpp"
#include "rtseis/utilities/windowFunctions.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities;

/*
std::vector<double> WindowFunctions::hamming(const int len)
{
    if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
    std::vector<double> window(len);
    hamming(len, window.data());
    return window;
}
*/

void WindowFunctions::hamming(const int len, double *winIn[])
{
    double *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("%s", "Length must be positive");}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    if (len == 2)
    {
        window[0] = 0.08;
        window[1] = 0.08;
        return;
    }
    ippsSet_64f(1, window, len);
    ippsWinHamming_64f_I(window, len);
    return;
}

/*
std::vector<float> WindowFunctions::hamming(const int len)
{
    if (len < 1)
    {
        throw std::invalid_argument("Length must be positive");
    }
    std::vector<float> window(len);
    hamming(len, window.data());
    return window;
}
*/

void WindowFunctions::hamming(const int len, float *winIn[])
{
    float *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("%s", "Length must be positive");}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    if (len == 2)
    {   
        window[0] = 0.08f;
        window[1] = 0.08f;
        return;
    }
    ippsSet_32f(1, window, len);
    ippsWinHamming_32f_I(window, len);
    return;
}
//============================================================================//

void WindowFunctions::hann(const int len, double *winIn[])
{
    double *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
    ippsSet_64f(1, window, len);
    ippsWinHann_64f_I(window, len);
    return;
}

void WindowFunctions::hann(const int len, float *winIn[])
{
    float *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
    ippsSet_32f(1, window, len);
    ippsWinHann_32f_I(window, len);
    return;
}

//============================================================================//

void WindowFunctions::bartlett(const int len, double *winIn[])
{
    double *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
    ippsSet_64f(1, window, len);
    ippsWinBartlett_64f_I(window, len);
    return;
}

void WindowFunctions::bartlett(const int len, float *winIn[])
{
    float *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
    ippsSet_32f(1, window, len);
    ippsWinBartlett_32f_I(window, len);
    return;
}

//============================================================================//

void WindowFunctions::sine(const int len, double *winIn[])
{
    double *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    // This will put zeros at end points of sine taper
    double pidmi = M_PI/static_cast<double> (len - 1);
    window[0] = 0;
    #pragma omp simd
    for (int i=1; i<len-1; i++)
    {
        window[i] = std::sin(pidmi*(static_cast<double> (i)));
    }
    window[len-1] = 0;
/*
    // I don't agree with scipy's implementation.
    double pidmi = M_PI/static_cast<double> (len);
    #pragma omp simd
    for (int i=0; i<len; i++)
    {
        window[i] = std::sin(pidmi*(static_cast<double> (i) + 0.5));
    }
*/
    return;
}

void WindowFunctions::sine(const int len, float *winIn[])
{
    float *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    // This will put zeros at end points of sine taper ala Wikipedia
    float pidmi = static_cast<float> (M_PI/static_cast<double> (len - 1));
    window[0] = 0;
    #pragma omp simd
    for (int i=1; i<len-1; i++)
    {   
        window[i] = std::sin(pidmi*(static_cast<float> (i)));
    }   
    window[len-1] = 0;
/*
    float pidmi = static_cast<float> (M_PI/static_cast<double> (len));
    #pragma omp simd
    for (int i=0; i<len; i++)
    {
        window[i] = std::sin(pidmi*(static_cast<float> (i) + 0.5f));
    }
*/
    return;
}

//============================================================================//

void WindowFunctions::kaiser(const int len, double *winIn[],
                             const double beta)
{
    double *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
#ifdef __STDCPP_MATH_SPEC_FUNCS__ //__cplusplus > 201402L
    double alpha = static_cast<double> (len - 1)/2;
    double i0betai = 1.0/std::cyl_bessel_i(0, beta);
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (int i=0; i<len; i++)
    {
        auto dn = static_cast<double> (i);
        double argSqrt = 1 - std::pow((dn - alpha)/alpha, 2);
        double arg = beta*std::sqrt(argSqrt);
        window[i] = std::abs(std::cyl_bessel_i(0, arg)*i0betai);
    }
#else
    double alpha = (2*beta)/static_cast<double> (len - 1);
    ippsSet_64f(1, window, len);
    IppStatus status = ippsWinKaiser_64f_I(window, len, alpha);
    if (status != ippStsNoErr)
    {
        RTSEIS_THROW_RTE("beta=%lf is too large for IPP", beta);
    }
#endif
    return;
}

void WindowFunctions::kaiser(const int len, float *winIn[],
                             const float beta)
{
    float *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    float alpha = (2*beta)/static_cast<float> (len - 1); 
    ippsSet_32f(1, window, len);
    IppStatus status = ippsWinKaiser_32f_I(window, len, alpha);
    if (status != ippStsNoErr)
    {
        RTSEIS_THROW_RTE("beta=%f is too large for IPP", beta);
    }
    return;
}

//============================================================================//

void WindowFunctions::blackman(const int len, double *winIn[])
{
    double *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
    ippsSet_64f(1, window, len);
    ippsWinBlackmanStd_64f_I(window, len);
    return;
}

void WindowFunctions::blackman(const int len, float *winIn[])
{
    float *window = *winIn;
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_THROW_IA("Length = %d must be positive", len);}
        RTSEIS_THROW_IA("%s", "window is NULL");
    }
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
    ippsSet_32f(1, window, len);
    ippsWinBlackmanStd_32f_I(window, len);
    return;
}
