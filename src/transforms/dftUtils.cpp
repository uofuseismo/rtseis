#include <string>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <valarray>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/transforms/utilities.hpp"
#include "rtseis/log.h"


using namespace RTSeis::Transforms;

namespace
{
#pragma omp declare simd
double rem(const double x, const double y)
{
    double rem = x - y*static_cast<double> (static_cast<int> (x/y));
    return rem;
}
#pragma omp declare simd
float rem(const float x, const float y)
{
    float rem = x - y*static_cast<float> (static_cast<int> (x/y));
    return rem;
}

double getMin(const int n, const double p[])
{
    double pmin;
    ippsMin_64f(p, n, &pmin);
    return pmin;
}
float getMin(const int n, const float p[])
{
    float pmin;
    ippsMin_32f(p, n, &pmin);
    return pmin;
}

}

/// Magnitude
template<typename T>
std::vector<T>
DFTUtilities::magnitude(const std::vector<std::complex<T>> &z) noexcept
{
    std::vector<T> mag;
    if (z.empty()){return mag;}
    mag.resize(z.size());
    auto magPtr = mag.data();
    magnitude(static_cast<int> (z.size()), z.data(), &magPtr);
    return mag;
}

template<>
void DFTUtilities::magnitude(const int n, const std::complex<double> *z,
                             double *magIn[])
{
    if (n < 1){return;}
    auto mag = *magIn;
    if (z == nullptr){throw std::invalid_argument("z is NULL");}
    if (mag == nullptr){throw std::invalid_argument("magnitude is NULL");}
    auto pSrc = reinterpret_cast<const Ipp64fc *> (z);
    ippsMagnitude_64fc(pSrc, mag, n);
}

template<>
void DFTUtilities::magnitude(const int n, const std::complex<float> *z,
                             float *magIn[])
{
    if (n < 1){return;}
    auto mag = *magIn;
    if (z == nullptr){throw std::invalid_argument("z is NULL");}
    if (mag == nullptr){throw std::invalid_argument("magnitude is NULL");}
    auto pSrc = reinterpret_cast<const Ipp32fc *> (z);
    ippsMagnitude_32fc(pSrc, mag, n);
}

/// Unwrap phase angle
template<typename T>
std::vector<T> DFTUtilities::unwrap(const std::vector<T> &p,
                                    const T tol)
{
    std::vector<T> q;
    if (p.empty()){return q;}
    q.resize(p.size());
    auto qPtr = q.data();
    unwrap(static_cast<int> (p.size()), p.data(), &qPtr, tol); 
    return q;
}

template<typename T>
void DFTUtilities::unwrap(const int n, const T p[], T *qIn[],
                          const T tol)
{
    if (n <= 0){return;}
    auto q = *qIn;
    if (p == nullptr || q == nullptr)
    {
        if (p == nullptr){throw std::invalid_argument("p is NULL");}
        throw std::invalid_argument("q is NULL");
    }
    if (tol < 0)
    {
         throw std::invalid_argument("Tolerance = " + std::to_string(tol)
                                   + " cannot be negative");
    }
    auto pmin = getMin(n, p);
    const T twopi = static_cast<T> (2*M_PI);
    #pragma omp simd
    for (int i=0; i<n; i++)
    {
        q[i] = rem(p[i] - pmin, twopi) + pmin;
    }
    // Differentiate phases
    std::valarray<T> b(n);
    b[0] = q[0];
    #pragma omp simd
    for (int i=1; i<n; i++)
    {
        b[i] = q[i] - q[i-1];
    }
    T cumsume = 0;
    for (int i=0; i<n; i++)
    {
        // Locate jumps
        int ic = 0;
        if (b[i] > tol){ic =-1;}
        int id = 0;
        if (b[i] <-tol){id = 1;}
        // 2*pi jumps
        T e = twopi*static_cast<T> (ic + id);
        // Integrate to get corrections
        cumsume = cumsume + e; // Integrate to get corrections
        q[i] = q[i] + cumsume;
    }
}

/// Computes the phase
template<typename T>
std::vector<T>
DFTUtilities::phase(const std::vector<std::complex<T>> &z, 
                    const bool lwantDeg) noexcept
{
    std::vector<T> phi;
    if (z.empty()){return phi;} 
    phi.resize(z.size());
    auto phiPtr = phi.data();
    phase(static_cast<int> (z.size()), z.data(), &phiPtr, lwantDeg); 
    return phi;
}

/// Compute the phase
template<>
void DFTUtilities::phase(const int n, const std::complex<double> z[],
                         double *phiIn[],
                         const bool lwantDeg)
{
    if (n <= 0){return;} 
    auto phi = *phiIn;
    if (z == nullptr || phi == nullptr)
    {
        if (z == nullptr){throw std::invalid_argument("z is NULL");}
        throw std::invalid_argument("phi is NULL");
    }
    auto pSrc = reinterpret_cast<const Ipp64fc *> (z);
#ifdef DEBUG
    IppStatus status = ippsPhase_64fc(pSrc, phi, n); 
    assert(status == ippStsNoErr);
#else
    ippsPhase_64fc(pSrc, phi, n); 
#endif
    if (lwantDeg){ippsMulC_64f_I(180.0/M_PI, phi, n);}
}

template<>
void DFTUtilities::phase(const int n, const std::complex<float> z[],
                         float *phiIn[],
                         const bool lwantDeg)
{
    if (n <= 0){return;}
    auto phi = *phiIn;
    if (z == nullptr || phi == nullptr)
    {
        if (z == nullptr){throw std::invalid_argument("z is NULL");}
        throw std::invalid_argument("phi is NULL");
    }
    auto pSrc = reinterpret_cast<const Ipp32fc *> (z);
#ifdef DEBUG
    IppStatus status = ippsPhase_32fc(pSrc, phi, n);
    assert(status == ippStsNoErr);
#else
    ippsPhase_32fc(pSrc, phi, n);
#endif
    if (lwantDeg){ippsMulC_32f_I(static_cast<float> (180.0/M_PI), phi, n);}
    return;
}

/// Compute real/complex frequencies
template<typename T>
std::vector<T>
DFTUtilities::realToComplexDFTFrequencies(const int nSamples,
                                          const T samplingPeriod)
{
    if (nSamples < 1 || samplingPeriod <= 0)
    {
        if (nSamples < 1)
        {
            throw std::invalid_argument("nSamples = " + std::to_string(nSamples)
                                      + " must be positive");
        }
        throw std::invalid_argument("samplingPeriod = "
                                  + std::to_string(samplingPeriod)
                                  + " must be positive");
    }
    int nbins = nSamples/2 + 1;
    std::vector<T> freqs(nbins);
    auto freqsPtr = freqs.data();
    realToComplexDFTFrequencies(nSamples, samplingPeriod, nbins, &freqsPtr);
    return freqs;
}

template<typename T>
void DFTUtilities::realToComplexDFTFrequencies(const int nSamples,
                                               const T samplingPeriod,
                                               const int lengthFreqs,
                                               T *freqsIn[])
{
    if (nSamples < 1 || samplingPeriod <= 0)
    {
        if (nSamples < 1)
        {
            throw std::invalid_argument("nSamples = " + std::to_string(nSamples)
                                      + " must be positive");
        }
        throw std::invalid_argument("samplingPeriod = "
                                  + std::to_string(samplingPeriod)
                                  + " must be positive");
    }
    auto freqs = *freqsIn;
    if (freqs == nullptr){throw std::invalid_argument("freqs is NULL");}
    int nbins = nSamples/2 + 1;
    if (lengthFreqs < nbins)
    {
        throw std::invalid_argument("lengthFreqs = "
                                  + std::to_string(lengthFreqs)
                                  + " must be at least "
                                  + std::to_string(nbins));
    }  
    // Edge case
    if (nbins == 1)
    {
        freqs[0] = 0;
        return;
    }
    auto df
        = static_cast<T> (1.0/(static_cast<double> (nSamples)*samplingPeriod));
    #pragma omp simd
    for (auto i=0; i<nbins; ++i)
    {
        freqs[i] = df*i;
    }
}

int DFTUtilities::nextPowerOfTwo(const int n)
{
    if (n < 0)
    {
        throw std::invalid_argument("n = " + std::to_string(n)
                                 + " must be positive");
    }
    // Simple base cases
    if (n == 0){return 1;}
    if (n == 1){return 1;}
    // General computation
    double dn = static_cast<double> (n);
    // Get log2 and round up
    int itemp = static_cast<int> (std::ceil(std::log2(dn)));
    // Compute next power of 2 greater than n
    uint64_t n2t = static_cast<uint64_t> (std::pow(2, itemp));
    int n2 = static_cast<int> (n2t);
    // Catch any overflow problems associated with int32_t and uint64_t
    if (n2t != static_cast<uint64_t> (n2))
    {
        throw std::runtime_error("Overflow error in nextPowerOf2");
    }
    return n2;
}

/// fftshift
template<typename T> std::vector<T> 
RTSeis::Transforms::DFTUtilities::fftShift(const std::vector<T> &x)
{
    std::vector<T> y;
    if (x.empty()){return y;}
    int npts = static_cast<int> (x.size());
    y.resize(x.size());
    T *yPtr = y.data();
    fftShift(npts, x.data(), &yPtr); 
    return y;
}

template<typename T>
void RTSeis::Transforms::DFTUtilities::fftShift(
    const int n, const T x[], T *yIn[])
{
    T *y = *yIn; 
    // Base cases
    if (n == 1)
    {
        y[0] = x[0];
        return;
    }
    if (n == 2)
    {
        y[0] = x[1];
        y[1] = x[0];
        return;
    }
    // Figure out the copy indices
    int jf, ncopy1, ncopy2;
    int i1 = n/2;
    if (n%2 == 1)
    {
        i1 = n/2 + 1;
        ncopy1 = n - i1; // Tail shift
        ncopy2 = i1;
        jf = i1 - 1; 
    }
    else
    {
        i1 = n/2;
        ncopy1 = n/2;
        ncopy2 = ncopy1;
        jf = i1;  
    }
    // Now copy second half of x to start of y
    std::copy(&x[i1], &x[i1]+ncopy1, y);
    // And copy first half of x to second half of y
    std::copy(x, x+ncopy2, &y[jf]); 
}

///--------------------------------------------------------------------------///
///                           Instantiate Templates                          ///
///--------------------------------------------------------------------------///
template
std::vector<double> RTSeis::Transforms::DFTUtilities::magnitude(
    const std::vector<std::complex<double>> &z);
template
std::vector<float> RTSeis::Transforms::DFTUtilities::magnitude(
    const std::vector<std::complex<float>> &z);

template
std::vector<double> RTSeis::Transforms::DFTUtilities::unwrap(
    const std::vector<double> &z, const double tol);
template
std::vector<float> RTSeis::Transforms::DFTUtilities::unwrap(
    const std::vector<float> &z, const float tol);

template
void RTSeis::Transforms::DFTUtilities::unwrap(
    const int n, const double p[], double *qIn[], const double tol);
template
void RTSeis::Transforms::DFTUtilities::unwrap(
    const int n, const float p[], float *qIn[], const float tol);

template
std::vector<double> RTSeis::Transforms::DFTUtilities::phase(
    const std::vector<std::complex<double>> &z, const bool lwantDeg);
template
std::vector<float> RTSeis::Transforms::DFTUtilities::phase(
    const std::vector<std::complex<float>> &z, const bool lwantDeg);

template
std::vector<double> RTSeis::Transforms::DFTUtilities::realToComplexDFTFrequencies(
    const int nSamples, const double samplingPeriod);
template
std::vector<float> RTSeis::Transforms::DFTUtilities::realToComplexDFTFrequencies(
    const int nSamples, const float samplingPeriod);

template
void RTSeis::Transforms::DFTUtilities::realToComplexDFTFrequencies(
    const int nSamples, const double samplingPeriod,
    const int lengthFreqs, double *freqsIn[]);
template
void RTSeis::Transforms::DFTUtilities::realToComplexDFTFrequencies(
    const int nSamples, const float samplingPeriod,
    const int lengthFreqs, float *freqsIn[]);

template
std::vector<double> RTSeis::Transforms::DFTUtilities::fftShift<double>
    (const std::vector<double> &x);
template
void RTSeis::Transforms::DFTUtilities::fftShift<double>
    (const int n, const double x[], double *yIn[]);

template
std::vector<float> RTSeis::Transforms::DFTUtilities::fftShift<float>
    (const std::vector<float> &x);
template
void RTSeis::Transforms::DFTUtilities::fftShift<float>
    (const int n, const float x[], float *yIn[]);

template std::vector<std::complex<double>> 
RTSeis::Transforms::DFTUtilities::fftShift<std::complex<double>>
    (const std::vector<std::complex<double>> &x);
template
void RTSeis::Transforms::DFTUtilities::fftShift<std::complex<double>>
    (const int n, const std::complex<double> x[], std::complex<double> *yIn[]);

template std::vector<std::complex<float>>
RTSeis::Transforms::DFTUtilities::fftShift<std::complex<float>>
    (const std::vector<std::complex<float>> &x);
template
void RTSeis::Transforms::DFTUtilities::fftShift<std::complex<float>>
    (const int n, const std::complex<float> x[], std::complex<float> *yIn[]);

template std::vector<int>
RTSeis::Transforms::DFTUtilities::fftShift<int>
    (const std::vector<int> &x);
template
void RTSeis::Transforms::DFTUtilities::fftShift<int>
    (const int n, const int x[], int *yIn[]);

