#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cfloat>
#include <algorithm>
#define RTSEIS_LOGGING 1
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/design/response.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"
#include "rtseis/utilities/math/polynomial.hpp"
#include "rtseis/utilities/math/convolve.hpp"
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/log.h"
#include <ipps.h>


using namespace RTSeis::Utilities::FilterDesign;
using namespace RTSeis::Utilities::FilterRepresentations;

std::vector<std::complex<double>>
Response::freqs(const BA &ba, const std::vector<double> &w)
{
    size_t nw = w.size();
    std::vector<std::complex<double>> h;
    if (nw == 0){return h;}
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
        RTSEIS_THROW_IA("%s", "a is entirely 0; division by zero");
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
    Math::Polynomial::polyval(bz, s, hsNum);

    std::vector<std::complex<double>> az;
    az.resize(a.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<a.size(); i++){az[i] = std::complex<double> (a[i], 0);}
    std::vector<std::complex<double>> hsDen;
    Math::Polynomial::polyval(az, s, hsDen);
    // Compute the transfer function
    Math::VectorMath::divide(hsDen, hsNum, h);
    return h;
}

std::vector<std::complex<double>>
Response::freqz(const BA &ba, const std::vector<double> &w)
{
    size_t nw = w.size();
    std::vector<std::complex<double>> h;
    if (nw == 0){return h;} 
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
        RTSEIS_THROW_IA("%s", "a is entirely 0; division by zero");
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
    Math::Polynomial::polyval(bz, z, hzNum);

    std::vector<std::complex<double>> hzDen;
    Math::Polynomial::polyval(az, z, hzDen);
    // Compute the transfer function: H = Num/Den
    Math::VectorMath::divide(hzDen, hzNum, h); 
    return h;
}

std::vector<double>
Response::groupDelay(const BA &ba,
                    const std::vector<double> &w)
{
    std::vector<double> gd;
    // Nothing to do
    if (w.size() == 0){return gd;}
    // Get handles on data
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients();
    if (b.size() < 1 || a.size() < 1)
    {
        if (b.size() < 1){RTSEIS_THROW_IA("%s", "No elements in b");}
        if (a.size() < 1){RTSEIS_THROW_IA("%s", "No elements in a");}
        RTSEIS_THROW_IA("%s", "Invalid arguments");
    }
    // B/A = b[0]/A + b[1]/A + ... = b[0]*conj(A) + ... = b*conj(a)
    std::vector<double> c;
    // Correlation is convolution, however, a is reversed.
    c = Math::Convolve::correlate(b, a,
                                  Math::Convolve::Mode::FULL,
                                  Math::Convolve::Implementation::DIRECT);
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
    Math::Polynomial::polyval(zcr, z, num);

    std::vector<std::complex<double>> den;
    Math::Polynomial::polyval(zc, z, den);
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
    return gd;
}
