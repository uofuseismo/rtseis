#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/windowFunctions.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities;

void WindowFunctions::hamming(const int len, std::vector<double> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
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
    ippsSet_64f(1, window.data(), len);
    ippsWinHamming_64f_I(window.data(), len);
    return;
}

void WindowFunctions::hamming(const int len, std::vector<float> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
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
    ippsSet_32f(1, window.data(), len); 
    ippsWinHamming_32f_I(window.data(), len);
    return;
}

//============================================================================//

void WindowFunctions::hann(const int len, std::vector<double> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
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
    ippsSet_64f(1, window.data(), len);
    ippsWinHann_64f_I(window.data(), len);
    return;
}

void WindowFunctions::hann(const int len, std::vector<float> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
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
    ippsSet_32f(1, window.data(), len);
    ippsWinHann_32f_I(window.data(), len);
    return;
}

//============================================================================//

void WindowFunctions::bartlett(const int len, std::vector<double> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
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
    ippsSet_64f(1, window.data(), len);
    ippsWinBartlett_64f_I(window.data(), len);
    return;
}

void WindowFunctions::bartlett(const int len, std::vector<float> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
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
    ippsSet_32f(1, window.data(), len);
    ippsWinBartlett_32f_I(window.data(), len);
    return;
}

//============================================================================//

void WindowFunctions::sine(const int len, std::vector<double> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    double pidmi = M_PI/static_cast<double> (len);
    #pragma omp simd
    for (int i=0; i<len; i++)
    {
        window[i] = std::sin(pidmi*(static_cast<double> (i) + 0.5));
    }
    return;
}

void WindowFunctions::sine(const int len, std::vector<float> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    float pidmi = static_cast<float> (M_PI/static_cast<double> (len));
    #pragma omp simd
    for (int i=0; i<len; i++)
    {
        window[i] = std::sin(pidmi*(static_cast<float> (i) + 0.5f));
    }
    return;
}

//============================================================================//

void WindowFunctions::kaiser(const int len, std::vector<double> &window,
                             const double beta)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
#if __cplusplus > 201402L
    double alpha = static_cast<double> (len - 1)/2;
    double i0betai = 1.0/std::cyl_bessel_i(0, beta);
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (int i=0; i<len; i++)
    {
        double dn = static_cast<double> (i);
        double argSqrt = 1 - std::pow((dn - alpha)/alpha, 2);
        double arg = beta*std::sqrt(argSqrt);
        window[i] = std::abs(std::cyl_bessel_i(0, arg)*i0betai);
    }
#else
    double alpha = (2*beta)/static_cast<double> (len - 1);
    ippsSet_64f(1, window.data(), len);
    IppStatus status = ippsWinKaiser_64f_I(window.data(), len, alpha);
    if (status != ippStsNoErr)
    {
        window.resize(0);
        RTSEIS_ERRMSG("beta=%lf is too large", beta);
        std::invalid_argument("beta is too large");
        return;
    }
#endif
    return;
}

void WindowFunctions::kaiser(const int len, std::vector<float> &window,
                            const float beta)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
    if (len == 1)
    {
        window[0] = 1;
        return;
    }
    float alpha = (2*beta)/static_cast<float> (len - 1); 
    ippsSet_32f(1, window.data(), len);
    IppStatus status = ippsWinKaiser_32f_I(window.data(), len, alpha);
    if (status != ippStsNoErr)
    {
        window.resize(0);
        RTSEIS_ERRMSG("beta=%f too large", beta);
        std::invalid_argument("beta is too large"); 
        return; 
    }
    return;
}

//============================================================================//

void WindowFunctions::blackman(const int len, std::vector<double> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
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
    ippsSet_64f(1, window.data(), len);
    ippsWinBlackmanStd_64f_I(window.data(), len);
    return;
}

void WindowFunctions::blackman(const int len, std::vector<float> &window)
{
    if (len < 1)
    {
        window.resize(0);
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        throw std::invalid_argument("Length must be positive");
    }
    window.resize(len);
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
    ippsSet_32f(1, window.data(), len);
    ippsWinBlackmanStd_32f_I(window.data(), len);
    return;
}
