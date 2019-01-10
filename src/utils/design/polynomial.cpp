#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <complex>
#include <algorithm>
#include <float.h>
#define RTSEIS_LOGGING 1
#include "rtseis/utils/polynomial.hpp"
#include "rtseis/log.h"
#include <mkl_lapacke.h>

using namespace::RTSeis::Utils::Math;

/*!
 * @defgroup rtseis_utils_polynomial Polynomial
 * @brief Utility functions for polynomial handling.
 *        This code is originally from ISTI's ISCL and has been
 *        modified to conform with C++.  Function names have also been
 *        changed to conform with rtseis's naming convention.
 * @copyright ISTI distributed under the Apache 2 license.
 * @ingroup rtseis_utils
 */

/*!
 * @brief Evaluates the polynomial
 *        \f[
 *            p(x) = p_{n_{order}}
 *                 + x p_{n_{order}-1}
 *                 + \cdots
 *                 + x^{n_{order}} p_0
 *        \f]
 *        at points \f$ x_j, j=1,2,...,n_x \f$.
 * @param[in] p    The polynomial coefficients ordered such that the 
 *                 highest order coefficient comes first.  This has
 *                 dimension [order+1].
 * @param[in] x    The points at which to evaluate the polynomial.  This
 *                 has dimension [x.size()].
 * @param[out] y   \f$ y = p(x) \f$ evaluated at each \f$ x_i \f$.  This
 *                 has dimension [x.size()].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_polynomial 
 */
int Polynomial::polyval(const std::vector<std::complex<double>> p,
                        const std::vector<std::complex<double>> x,
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
        #pragma ivdep
        for (int i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i] + p[1];
        }
    }
    // Expand the quadratic case
    else if (norder == 2)
    {
        #pragma ivdep
        for (int i=0; i<nx; i++)
        {
            y[i] = p[2] + x[i]*(p[1] + x[i]*p[0]);
        }
    }
    // Expand the cubic case
    else if (norder == 3)
    {
        #pragma ivdep
        for (int i=0; i<nx; i++)
        {
            y[i] = p[3] + x[i]*(p[2] + x[i]*p[1] + x[i]*x[i]*p[0]);
        }   
    }
    // Be more general 
    else
    {
        #pragma ivdep
        for (int i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i];
        }
        for (int j=1; j<norder; j++)
        {
            #pragma ivdep
            for (int i=0; i<nx; i++)
            {
                y[i] = (p[j] + y[i])*x[i];
            }
        }
        #pragma ivdep
        for (int i=0; i<nx; i++)
        {
            y[i] = p[norder] + y[i];
        }
    }
    return 0;
}
/*!
 * @brief Evaluates the polynomial
 *        \f[
 *            p(x) = p_{n_{order}}
 *                 + x p_{n_{order}-1}
 *                 + \cdots
 *                 + x^{n_{order}} p_0
 *        \f]
 *        at points \f$ x_j, j=1,2,...,n_x \f$.
 * @param[in] p    The polynomial coefficients ordered such that the 
 *                 highest order coefficient comes first.  This has
 *                 dimension [order+1].
 * @param[in] x    The points at which to evaluate the polynomial.  This
 *                 has dimension [x.size()].
 * @param[out] y   \f$ y = p(x) \f$ evaluated at each \f$ x_i \f$.  This
 *                 has dimension [x.size()].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_polynomial 
 */
int Polynomial::polyval(const std::vector<double> p,
                        const std::vector<double> x,
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
        for (int i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i] + p[1];
        }
    }
    // Expand the quadratic case
    else if (norder == 2)
    {
        #pragma omp simd
        for (int i=0; i<nx; i++)
        {
            y[i] = p[2] + x[i]*(p[1] + x[i]*p[0]);
        }
    }
    // Expand the cubic case
    else if (norder == 3)
    {
        #pragma omp simd
        for (int i=0; i<nx; i++)
        {
            y[i] = p[3] + x[i]*(p[2] + x[i]*p[1] + x[i]*x[i]*p[0]);
        }   
    }
    // Be more general 
    else
    {   
        #pragma omp simd
        for (int i=0; i<nx; i++)
        {
            y[i] = p[0]*x[i];
        }
        for (int j=1; j<norder; j++)
        {
            #pragma omp simd
            for (int i=0; i<nx; i++)
            {
                y[i] = (p[j] + y[i])*x[i];
            }
        }
        #pragma omp simd
        for (int i=0; i<nx; i++)
        {
            y[i] = p[norder] + y[i];
        }
    }
    return 0;
}
/*!
 * @brief Returns a polynomial whos eroots are given by p.
 * @param[in] p    The polynomial roots.  
 * @param[in] y    The polynomial coefficients corresponding to the
 *                 roots of the given polynomial.  This has dimension
 *                 [p.size()+1] and is ordered so that the last coefficient
 *                 is the constant term and the first coefficient scales
 *                 the highest order polynomial.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_polynomial
 */
int Polynomial::poly(const std::vector<std::complex<double>> p,
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
        #pragma ivdep
        for (size_t j=i; j>=1; j--)
        {
            temp1[j] = y[j-1];
        }
        temp1[0] = zero;
        // -p_i*(a_0 + .... + a_n x^{i-1}) =-p_i a_0 - ... p_i a_n x^{i-1}
        // multiply
        #pragma ivdep
        for (size_t j=1; j<=i; j++)
        {
            temp2[j-1] =-y[j-1]*p[i-1];
        }
        temp2[i] = zero;
        // -p_i a_0 + p_i a_1 x + ... + p_i a_n x^{i-1} 
        //          -     a_0 x - ...                   - a_n x^i
        // difference previous two loops
        #pragma ivdep
        for (size_t j=1; j<=i+1; j++)
        {
            y[j-1] = temp1[j-1] + temp2[j-1];
        }   
    }
    // Purge any numerical junk
    #pragma ivdep
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
/*!
 * @copydoc Polynomial::poly
 * @ingroup rtseis_utils_polynomial
 */
int Polynomial::poly(const std::vector<double> p,
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
        #pragma ivdep
        for (size_t j=i; j>=1; j--)
        {
            temp1[j] = y[j-1];
        }
        temp1[0] = zero;
        // -p_i*(a_0 + .... + a_n x^{i-1}) =-p_i a_0 - ... p_i a_n x^{i-1}
        // multiply
        #pragma ivdep
        for (size_t j=1; j<=i; j++)
        {
            temp2[j-1] =-y[j-1]*p[i-1];
        }
        temp2[i] = zero;
        // -p_i a_0 + p_i a_1 x + ... + p_i a_n x^{i-1} 
        //          -     a_0 x - ...                   - a_n x^i
        // difference previous two loops
        #pragma ivdep
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
/*!
 * @brief Computes the roots of a polynomial:
 *        \f$ q(x) = c_0 x^p + c_1 x^{p-1} + \cdots + c_{p+1} \f$.
 *        where \f$ p \f$ is the polynomial order and \f$ coeffs[0] = c_0 \f$.
 *
 * @param[in] coeffs   The coefficients of the polynomial whose order 
 *                     is defined above.  Note, the coeffs[0] cannot be
 *                     0 and coeffs must have length at least 2.
 * @param[out] roots   The roots of the polynomial.  This has dimension
 *                     [coeffs.size() - 1].
 * @result 0 indicates success.
 * @ingroup rtseis_utils_polynomial
 */
int Polynomial::roots(const std::vector<double> coeffs,
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
    int lda = nord;
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
    double *wr = new double[static_cast<size_t> (n)];
    double *wi = new double[static_cast<size_t> (n)];
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
        #pragma ivdep
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
