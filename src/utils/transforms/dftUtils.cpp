#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <cmath>
#include <algorithm>
#define RTSEIS_LOGGING 1
#include "rtseis/utils/transforms.hpp"
#include "rtseis/log.h"
#include <ipps.h>


using namespace RTSeis::Utils::Transforms;

#pragma omp declare simd
static double rem(const double x, double y);

int DFTUtils::unwrap(const int n, const double p[], double q[],
                     const double tol)
{
    if (n <= 0){return 0;}
    if (p == nullptr || q == nullptr)
    {
        if (p == nullptr){RTSEIS_ERRMSG("%s", "p is NULL");}
        if (q == nullptr){RTSEIS_ERRMSG("%s", "q is NULL");}
        return -1;
    }
    if (tol < 0)
    {
        RTSEIS_ERRMSG("Tolerance = %lf cannot be negative", tol);
        return -1;
    }
    double pmin;
    ippsMin_64f(p, n, &pmin);
    const double twopi = 2*M_PI;
    #pragma omp simd
    for (int i=0; i<n; i++)
    {
        q[i] = rem(p[i] - pmin, twopi) + pmin;
    }
    // Differentiate phases
    double *b = new double[n];
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
    delete[] b;
    return 0;
}

int DFTUtils::phase(const int n, const std::complex<double> z[],
                    double phi[],
                    const bool lwantDeg)
{
    if (n <= 0){return 0;} 
    if (z == nullptr || phi == nullptr)
    {
        if (z == nullptr){RTSEIS_ERRMSG("%s", "z is NULL");}
        if (phi == nullptr){RTSEIS_ERRMSG("%s", "phi is NULL");}
        return -1;
    }
    const Ipp64fc *pSrc = reinterpret_cast<const Ipp64fc *> (z);
    IppStatus status = ippsPhase_64fc(pSrc, phi, n); 
    if (status != ippStsNoErr)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute phi");
        return -1;
    }
    if (lwantDeg)
    {
        ippsMulC_64f_I(180.0/M_PI, phi, n); 
    }
    return 0;
}

int DFTUtils::nextPow2(const int n)
{
    if (n < 0)
    {
        RTSEIS_ERRMSG("n=%d must be positive", n);
        return 0;
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
        RTSEIS_ERRMSG("%s", "Overflow error");
        return 0;
    }
    return n2;
}


#pragma omp declare simd
static double rem(const double x, double y)
{
    double rem = x - y*static_cast<double> (static_cast<int> (x/y));
    return rem;
}
