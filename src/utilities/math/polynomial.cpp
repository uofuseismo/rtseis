#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <cfloat>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/math/polynomial.hpp"
#include "rtseis/log.h"
#include <mkl_lapacke.h>

/*
 This source code is originally from ISTI's ISCL which is distributed under the
 Apache 2 license.  It has been heavily modified to conform to C++.
*/

using namespace::RTSeis::Utilities::Math;

int Polynomial::polyval(const std::vector<std::complex<double>> &p,
                        const std::vector<std::complex<double>> &x,
                        std::vector<std::complex<double>> &y)
{
    if (p.size() < 1)
    {
        RTSEIS_ERRMSG("%s", "No points in coefficients in p");
        y.resize(0);
        return -1; 
    }
    size_t norder = p.size() - 1;
    size_t nx = x.size();
    y.resize(nx);
    if (nx < 1){return 0;}
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
        for (size_t i=0; i<nx; i++)
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
        for (size_t i=0; i<nx; i++)
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
        for (size_t i=0; i<nx; i++)
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
        for (size_t i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i];
        }
        for (size_t j=1; j<norder; j++)
        {
#ifdef __INTEL_COMPILER
            #pragma ivdep
#endif
            for (size_t i=0; i<nx; i++)
            {
                y[i] = (p[j] + y[i])*x[i];
            }
        }
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (size_t i=0; i<nx; i++)
        {
            y[i] = p[norder] + y[i];
        }
    }
    return 0;
}

int Polynomial::polyval(const std::vector<double> &p,
                        const std::vector<double> &x,
                        std::vector<double> &y)
{
    if (p.size() < 1)
    {
        RTSEIS_ERRMSG("%s", "No points in coefficients in p");
        y.resize(0);
        return -1;
    }
    size_t norder = p.size() - 1;
    size_t nx = x.size();
    y.resize(nx);
    if (nx < 1){return 0;}
    // Expand the constant case 
    if (norder == 0)
    {
        std::fill(y.begin(), y.end(), p[0]);
    }
    // Expand the linear case
    else if (norder == 1)
    {
        #pragma omp simd
        for (size_t i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i] + p[1];
        }
    }
    // Expand the quadratic case
    else if (norder == 2)
    {
        #pragma omp simd
        for (size_t i=0; i<nx; i++)
        {
            y[i] = p[2] + x[i]*(p[1] + x[i]*p[0]);
        }
    }
    // Expand the cubic case
    else if (norder == 3)
    {
        #pragma omp simd
        for (size_t i=0; i<nx; i++)
        {
            y[i] = p[3] + x[i]*(p[2] + x[i]*p[1] + x[i]*x[i]*p[0]);
        }   
    }
    // Be more general 
    else
    {   
        #pragma omp simd
        for (size_t i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i];
        }
        for (size_t j=1; j<norder; j++)
        {
            #pragma omp simd
            for (size_t i=0; i<nx; i++)
            {
                y[i] = (p[j] + y[i])*x[i];
            }
        }
        #pragma omp simd
        for (size_t i=0; i<nx; i++)
        {
            y[i] = p[norder] + y[i];
        }
    }
    return 0;
}
//----------------------------------------------------------------------------//
int Polynomial::poly(const std::vector<std::complex<double>> &p,
                     std::vector<std::complex<double>> &y)
{
    size_t nord = p.size();
    y.resize(nord+1);
    // Special case
    if (nord == 0)
    {
        y[0] = std::complex<double> (1, 0);
        return 0;
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
    for (size_t i=2; i<=nord; i++)
    {
        // x*(a_0 + a_1 x + ... + a_n x^{n-1}) = a_0 x + a_1 x^2 + ... + a_n x^i
        // shift right
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (size_t j=i; j>=1; j--)
        {
            temp1[j] = y[j-1];
        }
        temp1[0] = zero;
        // -p_i*(a_0 + .... + a_n x^{i-1}) =-p_i a_0 - ... p_i a_n x^{i-1}
        // multiply
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (size_t j=1; j<=i; j++)
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
        for (size_t j=1; j<=i+1; j++)
        {
            y[j-1] = temp1[j-1] + temp2[j-1];
        }   
    }
    // Purge any numerical junk
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<nord+1; i++)
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
    return 0;
}

int Polynomial::poly(const std::vector<double> &p,
                     std::vector<double> &y)
{
    size_t nord = p.size();
    y.resize(nord+1);
    // Special case
    if (nord == 0)
    {
        y[0] = 1;
        return 0;
    } 
    // Linear case
    if (nord == 1)
    {
        y[0] = 1;
        y[1] =-p[0];
    } 
    const double zero = 0;
    const double one  = 1;
    std::vector<double> temp1(nord+1, zero);
    std::vector<double> temp2(nord+1, zero);
    // Initialize
    y[0] =-p[0];
    y[1] = one;   //y =-p_1 + x 
    for (size_t i=2; i<=nord; i++)
    {
        // x*(a_0 + a_1 x + ... + a_n x^{n-1}) = a_0 x + a_1 x^2 + ... + a_n x^i
        // shift right
        #pragma omp simd
        for (size_t j=i; j>=1; j--)
        {
            temp1[j] = y[j-1];
        }
        temp1[0] = zero;
        // -p_i*(a_0 + .... + a_n x^{i-1}) =-p_i a_0 - ... p_i a_n x^{i-1}
        // multiply
        #pragma omp simd
        for (size_t j=1; j<=i; j++)
        {
            temp2[j-1] =-y[j-1]*p[i-1];
        }
        temp2[i] = zero;
        // -p_i a_0 + p_i a_1 x + ... + p_i a_n x^{i-1} 
        //          -     a_0 x - ...                   - a_n x^i
        // difference previous two loops
        #pragma omp simd
        for (size_t j=1; j<=i+1; j++)
        {
            y[j-1] = temp1[j-1] + temp2[j-1];
        }
    }
    // Reverse y for consistency with MatLab
    std::reverse(y.begin(), y.end());
    // Free space
    temp1.clear();
    temp2.clear();
    return 0;
}

int Polynomial::roots(const std::vector<double> &coeffs,
                      std::vector<std::complex<double>> &roots) 
{
    int nc = static_cast<int> (coeffs.size());
    roots.resize(0);
    if (nc < 1)
    {
        RTSEIS_ERRMSG("%s", "No coefficients");
        return -1;
    }
    int nord = nc - 1;
    if (coeffs[0] == 0)
    {
        RTSEIS_ERRMSG("%s", "Highest order coefficient is zero");
        return -1;
    }
    // Set space for companion matrix
    int n   = nord;
    int lda = std::max(8, nord);
    double *a = new double[lda*lda];
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
    int ierr = 0;
    int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a, lda,
                             wr, wi, vl, ldvl, vr, ldvr);
    roots.resize(static_cast<size_t> (n));
    if (info != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute eigenvalues");
        ierr = 1;
        std::fill(roots.begin(), roots.end(), std::complex<double> (0, 0));
    }
    else
    {
#ifdef __INTEL_COMPILER
        #pragma ivdep
#endif
        for (int i=0; i<n; i++)
        {
            roots[i] = std::complex<double> (wr[i], wi[i]);
        }
    }
    delete[] wr;
    delete[] wi;
    delete[] a;
    return ierr;
}
