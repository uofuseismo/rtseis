#ifndef WINDOW_FUNCTIONS_HPP
#define WINDOW_FUNCTIONS_HPP
#include <string>
#include <cmath>
#include <stdexcept>
#ifdef WITH_IPP
#include <ipp/ipps.h>
#endif

namespace
{

void hamming(const int length, double *windowIn[])
{
    double *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    if (length == 2)
    {
        window[0] = 0.08;
        window[1] = 0.08;
        return;
    }
#ifdef WITH_IPP
    ippsSet_64f(1, window, length);
    ippsWinHamming_64f_I(window, length);
#else
    constexpr double a0{0.54};
    constexpr double a1{0.46};
    auto twopidmi = static_cast<double> ((2*M_PI)/(length - 1));
    for (int i = 0; i < length; ++i)
    {
        window[i] = a0 - a1*std::cos(twopidmi*i);
    }
#endif
}

void hamming(const int length, float *windowIn[])
{
    float *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    if (length == 2)
    {   
        window[0] = 0.08f;
        window[1] = 0.08f;
        return;
    }
#ifdef WITH_IPP
    ippsSet_32f(1, window, length);
    ippsWinHamming_32f_I(window, length);
#else
    constexpr float a0{0.54f};
    constexpr float a1{0.46f};
    auto twopidmi = static_cast<float> ((2*M_PI)/(length - 1));
    for (int i = 0; i < length; ++i)
    {
        window[i] = a0 - a1*std::cos(twopidmi*i);
    }
#endif
}
//============================================================================//

void hann(const int length, double *windowIn[])
{
    double *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    if (length == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
#ifdef WITH_IPP
    ippsSet_64f(1, window, length);
    ippsWinHann_64f_I(window, length);
#else
    constexpr double half{0.5};
    auto twopidmi = static_cast<double> ((2*M_PI)/(length - 1));
    window[0] = 0;
    for (int i = 1; i < length - 1; ++i)
    {
        window[i] = half*(1 - std::cos(twopidmi*i));
    }
    window[length - 1] = 0;
#endif
}

void hann(const int length, float *windowIn[])
{
    float *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    if (length == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
#ifdef WITH_IPP
    ippsSet_32f(1, window, length);
    ippsWinHann_32f_I(window, length);
#else
    constexpr float half{0.5};
    auto twopidmi = static_cast<float> ((2*M_PI)/(length - 1));
    window[0] = 0;
    for (int i = 1; i < length - 1; ++i)
    {
        window[i] = half*(1 - std::cos(twopidmi*i));
    }
    window[length - 1] = 0;
#endif
}

//============================================================================//

void bartlett(const int length, double *windowIn[])
{
    double *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    if (length == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
#ifdef WITH_IPP
    ippsSet_64f(1, window, length);
    ippsWinBartlett_64f_I(window, length);
#else
    auto n2 = static_cast<double> ((length - 1)/2.0);
    auto l2inv = static_cast<double> (2.0/(length - 1));
    window[0] = 0;
    for (int i = 1; i < length - 1; i++)
    {   
        window[i] = 1 - std::abs(i - n2)*l2inv;
    }   
    window[length - 1] = 0;
#endif
}

void bartlett(const int length, float *windowIn[])
{
    float *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    if (length == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
#ifdef WITH_IPP
    ippsSet_32f(1, window, length);
    ippsWinBartlett_32f_I(window, length);
#else
    auto n2 = static_cast<float> ((length - 1)/2.0f);
    auto l2inv = static_cast<float> (2.0f/(length - 1));
    window[0] = 0;
    for (int i = 1; i < length - 1; i++)
    {
        window[i] = 1 - std::abs(i - n2)*l2inv;
    }
    window[length - 1] = 0;
#endif
}

//============================================================================//

template<typename T>
void sine(const int length, T *windowIn[])
{
    T *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    // This will put zeros at end points of sine taper
    auto pidmi = static_cast<T> (M_PI/(length - 1));
    window[0] = 0;
    for (int i = 1; i < length - 1; i++)
    {
        window[i] = std::sin(pidmi*i);
    }
    window[length - 1] = 0;
}

//============================================================================//

void kaiser(const int length, double *windowIn[], const double beta)
{
    double *__restrict__ window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
#ifdef WITH_IPP
    double alpha = (2*beta)/static_cast<double> (length - 1);
    ippsSet_64f(1, window, length);
    IppStatus status = ippsWinKaiser_64f_I(window, length, alpha);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("beta = " + std::to_string(beta)
                               + " too large for IPP");
    }
#else
    double alpha = static_cast<double> (length - 1)/2;
    double i0betai = 1.0/std::cyl_bessel_i(0, beta);
    for (int i = 0; i < length; i++)
    {
        double argSqrt = 1 - std::pow((i - alpha)/alpha, 2);
        double arg = beta*std::sqrt(argSqrt);
        window[i] = std::abs(std::cyl_bessel_i(0, arg)*i0betai);
    }
#endif
}

void kaiser(const int length, float *windowIn[], const float beta)
{
    float *__restrict__ window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
#ifdef WITH_IPP
    float alpha = (2*beta)/static_cast<float> (length - 1);
    ippsSet_32f(1, window, length);
    IppStatus status = ippsWinKaiser_32f_I(window, length, alpha);
    if (status != ippStsNoErr)
    {
        throw std::runtime_error("beta = " + std::to_string(beta)
                               + " too large for IPP");
    }
#else
    float alpha = static_cast<float> (length - 1)/2;
    float i0betai = 1.0f/std::cyl_bessel_i(0, beta);
    for (int i = 0; i < length; i++)
    {
        float argSqrt = 1 - std::pow((i - alpha)/alpha, 2); 
        float arg = beta*std::sqrt(argSqrt);
        window[i] = std::abs(std::cyl_bessel_i(0, arg)*i0betai);
    }
#endif
}

//============================================================================//

void blackman(const int length, double *windowIn[])
{
    double *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    if (length == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
#ifdef WITH_IPP
    ippsSet_64f(1, window, length);
    ippsWinBlackmanStd_64f_I(window, length);
#else
    // Coefficients for the "not very serious proposal" Blackman window
    constexpr double alpha{0.16};
    constexpr double a0{(1 - alpha)/2};
    constexpr double a1{0.5};
    constexpr double a2{alpha/2};
    auto twopidmi  = static_cast<double> ((2*M_PI)/(length - 1));
    auto fourpidmi = static_cast<double> ((4*M_PI)/(length - 1));
    window[0] = 0;
    for (int i = 1; i < length - 1; ++i)
    {
        window[i] = a0 - a1*std::cos(twopidmi*i) + a2*std::cos(fourpidmi*i);
    }
    window[length - 1] = 0;
#endif
}

void blackman(const int length, float *windowIn[])
{
    float *window = *windowIn;
    if (length < 1 || window == nullptr)
    {
        if (length < 1)
        {
            throw std::invalid_argument("Length = " + std::to_string(length)
                                      + " must be positive");
        }
        throw std::invalid_argument("window is NULL");
    }
    if (length == 1)
    {
        window[0] = 1;
        return;
    }
    if (length == 2)
    {
        window[0] = 0;
        window[1] = 0;
        return;
    }
#ifdef WITH_IPP
    ippsSet_32f(1, window, length);
    ippsWinBlackmanStd_32f_I(window, length);
#else
    constexpr float alpha{0.16f};
    constexpr float a0{(1 - alpha)/2};
    constexpr float a1{0.5f};
    constexpr float a2{alpha/2};
    auto twopidmi  = static_cast<float> ((2*M_PI)/(length - 1));
    auto fourpidmi = static_cast<float> ((4*M_PI)/(length - 1));
    window[0] = 0;
    for (int i = 1; i < length - 1; ++i)
    {
        window[i] = a0 - a1*std::cos(twopidmi*i) + a2*std::cos(fourpidmi*i);
    }   
    window[length - 1] = 0;
#endif
}

}
#endif
