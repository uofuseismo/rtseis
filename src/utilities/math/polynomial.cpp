#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <cfloat>
#define RTSEIS_LOGGING 1
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/math/polynomial.hpp"
#include "rtseis/log.h"
#include <mkl_lapacke.h>

/*
 This source code is originally from ISTI's ISCL which is distributed under the
 Apache 2 license.  It has been heavily modified to conform to C++.
*/

using namespace::RTSeis::Utilities::Math;

std::vector<std::complex<double>> 
Polynomial::polyval(const std::vector<std::complex<double>> &p,
                    const std::vector<std::complex<double>> &x)
{
    if (p.empty())
    {
        RTSEIS_THROW_IA("%s", "No points in coefficients in p");
    }
    auto norder = p.size() - 1;
    int nx = static_cast<int> (x.size());
    std::vector<std::complex<double>> y(nx);
    if (nx < 1){return y;}
    // Expand the constant case 
    if (norder == 0)
    {
        std::fill(y.begin(), y.end(), p[0]);
    }
    // Expand the linear case
    else if (norder == 1)
    {
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i] + p[1];
        }
    }
    // Expand the quadratic case
    else if (norder == 2)
    {
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[2] + x[i]*(p[1] + x[i]*p[0]);
        }
    }
    // Expand the cubic case
    else if (norder == 3)
    {
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[3] + x[i]*(p[2] + x[i]*p[1] + x[i]*x[i]*p[0]);
        }   
    }
    // Be more general 
    else
    {
#ifdef __INTEL_COMPILER 
        #pragma ivdep
#endif
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i];
        }
        for (auto j=1; j<norder; j++)
        {
#ifdef __INTEL_COMPILER
            #pragma ivdep
#endif
            for (auto i=0; i<nx; i++)
            {
                y[i] = (p[j] + y[i])*x[i];
            }
        }
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[norder] + y[i];
        }
    }
    return y;
}

std::vector<double>
Polynomial::polyval(const std::vector<double> &p,
                    const std::vector<double> &x)
{
    if (p.empty())
    {
        RTSEIS_THROW_IA("%s", "No points in coefficients in p");
    }
    int norder = static_cast<int> (p.size()) - 1;
    int nx = static_cast<int> (x.size());
    std::vector<double> y(nx);
    if (nx < 1){return y;}
    // Expand the constant case 
    if (norder == 0)
    {
        std::fill(y.begin(), y.end(), p[0]);
    }
    // Expand the linear case
    else if (norder == 1)
    {
        #pragma omp simd
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i] + p[1];
        }
    }
    // Expand the quadratic case
    else if (norder == 2)
    {
        #pragma omp simd
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[2] + x[i]*(p[1] + x[i]*p[0]);
        }
    }
    // Expand the cubic case
    else if (norder == 3)
    {
        #pragma omp simd
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[3] + x[i]*(p[2] + x[i]*p[1] + x[i]*x[i]*p[0]);
        }   
    }
    // Be more general 
    else
    {   
        #pragma omp simd
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i];
        }
        for (auto j=1; j<norder; j++)
        {
            #pragma omp simd
            for (auto i=0; i<nx; i++)
            {
                y[i] = (p[j] + y[i])*x[i];
            }
        }
        #pragma omp simd
        for (auto i=0; i<nx; i++)
        {
            y[i] = p[norder] + y[i];
        }
    }
    return y;
}
//----------------------------------------------------------------------------//
std::vector<std::complex<double>> 
Polynomial::poly(const std::vector<std::complex<double>> &p) noexcept
{
    int nord = static_cast<int> (p.size());
    std::vector<std::complex<double>> y(nord+1);
    // Special case
    if (nord == 0)
    {
        y[0] = std::complex<double> (1, 0);
        return y;
    } 
    // Linear case
    if (nord == 1)
    {
        y[0] = std::complex<double> (1, 0);
        y[1] =-p[0];
    } 
    const std::complex<double> zero(0, 0);
    const std::complex<double> zone(1, 0);
    std::vector<std::complex<double>> temp1(nord+1, zero);
    std::vector<std::complex<double>> temp2(nord+1, zero);
    // Initialize
    y[0] =-p[0];
    y[1] = zone;  //y =-p_1 + x 
    for (auto i=2; i<=nord; i++)
    {
        // x*(a_0 + a_1 x + ... + a_n x^{n-1}) = a_0 x + a_1 x^2 + ... + a_n x^i
        // shift right
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (auto j=i; j>=1; j--)
        {
            temp1[j] = y[j-1];
        }
        temp1[0] = zero;
        // -p_i*(a_0 + .... + a_n x^{i-1}) =-p_i a_0 - ... p_i a_n x^{i-1}
        // multiply
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (auto j=1; j<=i; j++)
        {
            temp2[j-1] =-y[j-1]*p[i-1];
        }
        temp2[i] = zero;
        // -p_i a_0 + p_i a_1 x + ... + p_i a_n x^{i-1} 
        //          -     a_0 x - ...                   - a_n x^i
        // difference previous two loops
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (auto j=1; j<=i+1; j++)
        {
            y[j-1] = temp1[j-1] + temp2[j-1];
        }   
    }
    // Purge any numerical junk
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (auto i=0; i<nord+1; i++)
    {
        if (std::abs(std::imag(y[i])) < DBL_EPSILON)
        {
            y[i] = std::complex<double> (std::real(y[i]), 0);
        }
    }
    // Reverse y for consistency with MatLab
    std::reverse(y.begin(), y.end());
    // Free space
    temp1.clear();
    temp2.clear();
    return y;
}

std::vector<double>
Polynomial::poly(const std::vector<double> &p) noexcept
{
    auto nord = p.size();
    std::vector<double> y(nord+1);
    // Special case
    if (nord == 0)
    {
        y[0] = 1;
        return y;
    } 
    // Linear case
    if (nord == 1)
    {
        y[0] = 1;
        y[1] =-p[0];
    } 
    constexpr double zero = 0;
    constexpr double one  = 1;
    std::vector<double> temp1(nord+1, zero);
    std::vector<double> temp2(nord+1, zero);
    // Initialize
    y[0] =-p[0];
    y[1] = one;   //y =-p_1 + x 
    for (auto i=2; i<=nord; i++)
    {
        // x*(a_0 + a_1 x + ... + a_n x^{n-1}) = a_0 x + a_1 x^2 + ... + a_n x^i
        // shift right
        #pragma omp simd
        for (auto j=i; j>=1; j--)
        {
            temp1[j] = y[j-1];
        }
        temp1[0] = zero;
        // -p_i*(a_0 + .... + a_n x^{i-1}) =-p_i a_0 - ... p_i a_n x^{i-1}
        // multiply
        #pragma omp simd
        for (auto j=1; j<=i; j++)
        {
            temp2[j-1] =-y[j-1]*p[i-1];
        }
        temp2[i] = zero;
        // -p_i a_0 + p_i a_1 x + ... + p_i a_n x^{i-1} 
        //          -     a_0 x - ...                   - a_n x^i
        // difference previous two loops
        #pragma omp simd
        for (auto j=1; j<=i+1; j++)
        {
            y[j-1] = temp1[j-1] + temp2[j-1];
        }
    }
    // Reverse y for consistency with MatLab
    std::reverse(y.begin(), y.end());
    // Free space
    temp1.clear();
    temp2.clear();
    return y;
}

std::vector<std::complex<double>> 
    Polynomial::roots(const std::vector<double> &coeffs)
{
    auto nc = static_cast<int> (coeffs.size());
    if (nc < 1)
    {
        RTSEIS_THROW_IA("%s", "No coefficients");
    }
    int nord = nc - 1;
    if (coeffs[0] == 0)
    {
        RTSEIS_THROW_IA("%s", "Highest order coefficient is zero");
    }
    // Set space for companion matrix
    int n   = nord;
    int lda = std::max(8, nord);
    auto *a = new double[lda*lda];
    memset(a, 0, static_cast<size_t> (lda*lda)*sizeof(double));
    double ami = 1.0/coeffs[0]; //coefficient on highest order term
    // Fill out the non-zeros of the companion matrix
    for (int i=1; i<nord+1; i++)
    {
        int indx = (i - 1)*lda + 0;
        a[indx] =-coeffs[i]*ami;
        // One in subdiagonal 
        if (i < nord)
        {
            indx = (i - 1)*lda + i;
            a[indx] = 1.0;
        }
    }
    double *wr = new double[static_cast<size_t> (std::max(32, n))];
    double *wi = new double[static_cast<size_t> (std::max(32, n))];
    double vl[1] = {0};
    double vr[1] = {0}; 
    const int ldvl = 1;
    const int ldvr = 1;
#ifdef DEBUG
    int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a, lda,
                             wr, wi, vl, ldvl, vr, ldvr);
    cassert(info == 0);     
#else
    LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a, lda,
                  wr, wi, vl, ldvl, vr, ldvr);
#endif
    std::vector<std::complex<double>> roots(n);
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (auto i=0; i<n; i++)
    {
        roots[i] = std::complex<double> (wr[i], wi[i]);
    }
    delete[] wr;
    delete[] wi;
    delete[] a;
    return roots;
}
