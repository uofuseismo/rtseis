#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cfloat>
#include <algorithm>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/design/response.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"
#include "rtseis/utilities/math/polynomial.hpp"
#include "rtseis/utilities/math/convolve.hpp"
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/log.h"
#include <ipps.h>

/*
 This source code is originally from ISTI's ISCL which is distributed under the
 Apache 2 license.  It has been heavily modified to conform to C++.
*/

using namespace RTSeis::Utilities::FilterDesign;
using namespace RTSeis::Utilities::FilterRepresentations;

int Response::freqs(const BA &ba, const std::vector<double> &w,
                    std::vector<std::complex<double>> &h)
{
    size_t nw = w.size();
    h.resize(0);
    if (nw == 0){return 0;}
    // Get a handle on the numerator and denominator coefficients
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients();
    bool lzero = true;
    for (size_t ia=0; ia<a.size(); ia++)
    {
        if (a[ia] != 0)
        {
            lzero = false;
            break;
        }
    } 
    if (lzero)
    {
        RTSEIS_ERRMSG("%s", "a is entirely 0; division by zero");
        return -1;
    }
    // Compute s = i \omega
    std::vector<std::complex<double>> s;
    s.resize(nw);
    for (size_t i=0; i<nw; i++)
    {
        s[i] = std::complex<double> (0, w[i]); // s = i omega
    }
    // Evaluate the numerator and denominator polynoimals
    std::vector<std::complex<double>> bz;
    bz.resize(b.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<b.size(); i++){bz[i] = std::complex<double> (b[i], 0);}
    std::vector<std::complex<double>> hsNum;
    int ierr = Math::Polynomial::polyval(bz, s, hsNum);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hsNum");
        return -1;
    }
    std::vector<std::complex<double>> hsDen;
    std::vector<std::complex<double>> az;
    az.resize(a.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<a.size(); i++){az[i] = std::complex<double> (a[i], 0);}
    ierr = Math::Polynomial::polyval(az, s, hsDen);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hsDen");
        return -1;
    }
    // Compute the transfer function
    ierr = Math::VectorMath::divide(hsDen, hsNum, h);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Division failed");
        return -1;
    }
    return 0;
}

int Response::freqz(const BA &ba, const std::vector<double> &w,
                    std::vector<std::complex<double>> &h) 
{
    size_t nw = w.size();
    h.resize(0);
    if (nw == 0){return 0;} 
    // Get a handle on the numerator and denominator coefficients
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients();
    bool lzero = true;
    for (size_t ia=0; ia<a.size(); ia++)
    {   
        if (a[ia] != 0)
        {
            lzero = false;
            break;
        }
    }   
    if (lzero)
    {   
        RTSEIS_ERRMSG("%s", "a is entirely 0; division by zero");
        return -1; 
    }
    // Copy numerators and denominators in reverse order for consistency with
    // polyval and then compute z.  The polyval convention requires us to
    // compute z=e^{-i\omega} to compensate.  Normally, this is would
    // be z^{i \omega}.
    size_t nb = b.size();
    std::vector<std::complex<double>> bz;
    bz.resize(nb);
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<nb; i++)
    {
        bz[i] = b[nb-1-i];
    }
    size_t na = a.size();
    std::vector<std::complex<double>> az;
    az.resize(a.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<na; i++)
    {
        az[i] = a[na-1-i];
    }
    std::vector<std::complex<double>> z;
    z.resize(nw);
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<nw; i++)
    {
        z[i] = std::exp(std::complex<double> (0, -w[i]));
    }
    // Evaluate the numerator and denominator polynoimals
    std::vector<std::complex<double>> hzNum;
    int ierr = Math::Polynomial::polyval(bz, z, hzNum);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hzNum");
        return -1; 
    }
    std::vector<std::complex<double>> hzDen;
    ierr = Math::Polynomial::polyval(az, z, hzDen);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute hzDen");
        return -1;
    }
    // Compute the transfer function: H = Num/Den
    ierr = Math::VectorMath::divide(hzDen, hzNum, h); 
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Division failed");
        return -1;
    }
    return 0;
}

/*
int Response::groupDelay(const BA &ba,
                         std::vector<double> &gd,
                         const int n,
                         const bool lwhole)
{
    gd.resize(0);
    if (n < 0)
    {
        RTSEIS_ERRMSG("%s", "n must be positive");
        return -1;
    }
    if (n == 0){return 0;}
    if (fs <= 0)
    {
        RTSEIS_ERRMSG("fs = %lf must be postiive", fs);
        return -1;
    }
    // Discretize
    std::vector<double> w(n);
    double di = M_PI/static_cast<double> (n);
    if (lwhole){di = (2*M_PI)/static_cast<double> (n);}
    if (n == 1)
    {
        w[0] = di;
    }
    else
    {
        #pragma omp simd
        for (int i=0; i<n; i++)
        {
            w[i] = di*static_cast<double> (i);
        }
    }
    int ierr = groupDelay(ba, w, gd);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute group delay");
        return -1;
    }
    return 0;
}
*/

int Response::groupDelay(const BA &ba,
                         const std::vector<double> &w,
                         std::vector<double> &gd)
{
    gd.resize(0);
    // Nothing to do
    if (w.size() == 0){return 0;}
    // Get handles on data
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients();
    if (b.size() < 1 || a.size() < 1)
    {
        if (b.size() < 1){RTSEIS_ERRMSG("%s", "No elements in b");}
        if (a.size() < 1){RTSEIS_ERRMSG("%s", "No elements in a");}
        return -1;
    }
    // B/A = b[0]/A + b[1]/A + ... = b[0]*conj(A) + ... = b*conj(a)
    int ierr;
    std::vector<double> c;
    // Correlation is convolution, however, a is reversed.
    try
    {
        Math::Convolve::correlate(b, a, c,
                                  Math::Convolve::Mode::FULL,
                                  Math::Convolve::Implementation::DIRECT);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("%s", ia.what());
        return -1;
    }
    // Differentiate b*conj(a).  Note, that the order is reversed to conform
    // to polyval.
    int nc = static_cast<int> (c.size());
    std::vector<std::complex<double>> zc(nc);
    std::vector<std::complex<double>> zcr(nc);
    for (int i=0; i<static_cast<int>(c.size()); i++)
    {
        zc[ nc-1-i] = std::complex<double> (c[i], 0);
        // Differentiate with power rule
        zcr[nc-1-i] = std::complex<double> (c[i]*static_cast<double> (i), 0);
    }
    std::vector<std::complex<double>> z(w.size());
    #pragma omp simd
    for (size_t i=0; i<w.size(); i++)
    {
        std::complex<double> arg(0,-w[i]);
        z[i] = std::exp(arg);
    }
    // Evaluate the numerator and denominator polynomials at angular frequencies
    std::vector<std::complex<double>> num;
    ierr = Math::Polynomial::polyval(zcr, z, num);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute numerator");
        return -1;
    }
    std::vector<std::complex<double>> den;
    ierr = Math::Polynomial::polyval(zc, z, den);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed to compute numerator");
        return -1; 
    }
    // Check for singular elements
    for (size_t i=0; i<w.size(); i++)
    {
        if (std::abs(den[i]) < 10*DBL_EPSILON)
        {
            RTSEIS_WARNMSG("Group delay is singular at %lf", w[i]);
            num[i] = std::complex<double> (0, 0);
            den[i] = std::complex<double> (1, 0);
        }
    }
    // Divide gd = num/den
    std::vector<std::complex<double>> zgd;
    Math::VectorMath::divide(den, num, zgd);
    // Take real part and shift
    gd.resize(zgd.size());
    double shift = static_cast<double> (a.size() - 1);
    #pragma omp simd
    for (int i=0; i<static_cast<int> (gd.size()); i++)
    {
        gd[i] = std::real(zgd[i]) - shift;
    }
    return 0;
}
