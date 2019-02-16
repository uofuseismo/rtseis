#include <stdio.h>
#include <stdlib.h>
#include <ipps.h>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/windowFunctions.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities;

int WindowFunctions::hamming(const int len, double window[])
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    if (len == 2)
    {
        window[0] = 0.08;
        window[1] = 0.08;
        return 0;
    }
    ippsSet_64f(1, window, len);
    ippsWinHamming_64f_I(window, len);
    return 0;
}

int WindowFunctions::hamming(const int len, float window[])
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    if (len == 2)
    {
        window[0] = 0.08f;
        window[1] = 0.08f;
        return 0;
    }
    ippsSet_32f(1, window, len); 
    ippsWinHamming_32f_I(window, len);
    return 0;
}

//============================================================================//

int WindowFunctions::hann(const int len, double window[])
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return 0;
    }
    ippsSet_64f(1, window, len);
    ippsWinHann_64f_I(window, len);
    return 0;
}

int WindowFunctions::hann(const int len, float window[])
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return 0;
    }
    ippsSet_32f(1, window, len);
    ippsWinHann_32f_I(window, len);
    return 0;
}

//============================================================================//

int WindowFunctions::bartlett(const int len, double window[])
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return 0;
    }
    ippsSet_64f(1, window, len);
    ippsWinBartlett_64f_I(window, len);
    return 0;
}

int WindowFunctions::bartlett(const int len, float window[])
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return 0;
    }
    ippsSet_32f(1, window, len);
    ippsWinBartlett_32f_I(window, len);
    return 0;
}

//============================================================================//

int WindowFunctions::kaiser(const int len, double window[],
                            const double beta)
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
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
    ippsSet_64f(1, window, len);
    IppStatus status = ippsWinKaiser_64f_I(window, len, alpha);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Internal failure");
        return -1;
    }
#endif
    return 0;
}

int WindowFunctions::kaiser(const int len, float window[],
                            const float beta)
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1; 
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    float alpha = (2*beta)/static_cast<float> (len - 1); 
    ippsSet_32f(1, window, len);
    IppStatus status = ippsWinKaiser_32f_I(window, len, alpha);
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Internal failure");
        return -1; 
    }
    return 0;
}

//============================================================================//

int WindowFunctions::blackman(const int len, double window[])
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return 0;
    }
    ippsSet_64f(1, window, len);
    ippsWinBlackmanStd_64f_I(window, len);
    return 0;
}

int WindowFunctions::blackman(const int len, float window[])
{
    if (len < 1 || window == nullptr)
    {
        if (len < 1){RTSEIS_ERRMSG("Length = %d must be positive", len);}
        if (window == nullptr){RTSEIS_ERRMSG("%s", "window is NULL");}
        return -1;
    }
    if (len == 1)
    {
        window[0] = 1;
        return 0;
    }
    if (len == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return 0;
    }
    ippsSet_32f(1, window, len);
    ippsWinBlackmanStd_32f_I(window, len);
    return 0;
}
