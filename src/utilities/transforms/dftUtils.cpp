#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <valarray>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/transforms/utilities.hpp"
#include "rtseis/log.h"


using namespace RTSeis::Utilities::Transforms;

#pragma omp declare simd
static double rem(const double x, double y);

std::vector<double> DFTUtilities::unwrap(const std::vector<double> &p,
                                         const double tol) 
{
    std::vector<double> q;
    if (p.empty()){return q;}
    q.resize(p.size());
    double *qPtr = q.data();
    unwrap(static_cast<int> (p.size()), p.data(), &qPtr, tol); 
    return q;
}

void DFTUtilities::unwrap(const int n, const double p[], double *qIn[],
                          const double tol)
{
    if (n <= 0){return;}
    double *q = *qIn;
    if (p == nullptr || q == nullptr)
    {
        if (p == nullptr){RTSEIS_THROW_IA("%s", "p is NULL");}
        RTSEIS_THROW_IA("%s", "q is NULL");
    }
    if (tol < 0){RTSEIS_THROW_IA("Tolerance = %lf cannot be negative", tol);}
    double pmin;
    ippsMin_64f(p, n, &pmin);
    const double twopi = 2*M_PI;
    #pragma omp simd
    for (int i=0; i<n; i++)
    {
        q[i] = rem(p[i] - pmin, twopi) + pmin;
    }
    // Differentiate phases
    std::valarray<double> b(n);
    b[0] = q[0];
    #pragma omp simd
    for (int i=1; i<n; i++)
    {
        b[i] = q[i] - q[i-1];
    }
    double cumsume = 0.0;
    for (int i=0; i<n; i++)
    {
        // Locate jumps
        int ic = 0;
        if (b[i] > tol){ic =-1;}
        int id = 0;
        if (b[i] <-tol){id = 1;}
        // 2*pi jumps
        double e = twopi*static_cast<double> (ic + id);
        // Integrate to get corrections
        cumsume = cumsume + e; // Integrate to get corrections
        q[i] = q[i] + cumsume;
    }
    return;
}

std::vector<double>
DFTUtilities::phase(const std::vector<std::complex<double>> &z, 
                    const bool lwantDeg)
{
    std::vector<double> phi;
    if (z.empty()){return phi;} 
    phi.resize(z.size());
    double *phiPtr = phi.data();
    phase(static_cast<int> (z.size()), z.data(), &phiPtr, lwantDeg); 
    return phi;
}

void DFTUtilities::phase(const int n, const std::complex<double> z[],
                         double *phiIn[],
                         const bool lwantDeg)
{
    if (n <= 0){return;} 
    double *phi = *phiIn;
    if (z == nullptr || phi == nullptr)
    {
        if (z == nullptr){RTSEIS_THROW_IA("%s", "z is NULL");}
        RTSEIS_THROW_IA("%s", "phi is NULL");
    }
    const Ipp64fc *pSrc = reinterpret_cast<const Ipp64fc *> (z);
#ifdef DEBUG
    IppStatus status = ippsPhase_64fc(pSrc, phi, n); 
    assert(status == ippStsNoErr);
#else
    ippsPhase_64fc(pSrc, phi, n); 
#endif
    if (lwantDeg){ippsMulC_64f_I(180.0/M_PI, phi, n);}
    return;
}

std::vector<double>
DFTUtilities::realToComplexDFTFrequencies(const int nSamples,
                                          const double samplingPeriod)
{
    std::vector<double> freqs;
    if (nSamples < 1 || samplingPeriod <= 0)
    {
        if (nSamples < 1)
        {
            RTSEIS_THROW_IA("nSamples = %d must be positive", nSamples);
        }
        RTSEIS_THROW_IA("samplingPeriod = %lf must be positive",
                        samplingPeriod);
    }
    int nbins = nSamples/2 + 1;
    freqs.resize(nbins);
    double *freqsPtr = freqs.data();
    realToComplexDFTFrequencies(nSamples, samplingPeriod, nbins, &freqsPtr);
    return freqs;
}

void DFTUtilities::realToComplexDFTFrequencies(const int nSamples,
                                               const double samplingPeriod,
                                               const int lengthFreqs,
                                               double *freqsIn[])
{
    if (nSamples < 1 || samplingPeriod <= 0)
    {
        if (nSamples < 1)
        {
            RTSEIS_THROW_IA("nSamples = %d must be positive", nSamples);
        }
        RTSEIS_THROW_IA("samplingPeriod = %lf must be positive",
                        samplingPeriod);
    }
    double *freqs = *freqsIn;
    if (freqs == nullptr)
    {
        RTSEIS_THROW_IA("%s", "freqs is NULL");
    } 
    int nbins = nSamples/2 + 1;
    if (lengthFreqs < nbins)
    {
        RTSEIS_THROW_IA("lengthFreqs = %d must be at least = %d",
                        lengthFreqs, nbins);
    }  
    // Edge case
    if (nbins == 1)
    {
        freqs[0] = 0;
        return;
    }
    auto df = 1.0/(static_cast<double> (nSamples)*samplingPeriod);
    #pragma omp simd
    for (auto i=0; i<nbins; ++i)
    {
        freqs[i] = df*i;
    }
}

int DFTUtilities::nextPow2(const int n)
{
    if (n < 0)
    {
        RTSEIS_THROW_IA("n=%d must be positive", n);
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
        RTSEIS_THROW_RTE("%s", "Overflow error");
    }
    return n2;
}


#pragma omp declare simd
static double rem(const double x, double y)
{
    double rem = x - y*static_cast<double> (static_cast<int> (x/y));
    return rem;
}
