#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/filterDesign/iir.hpp"
#include "rtseis/utilities/filterDesign/enums.hpp"
#include "rtseis/utilities/filterDesign/analogPrototype.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"
#include "rtseis/utilities/math/polynomial.hpp"
#include "private/throw.hpp"
#include "rtseis/log.h"

//#define DEBUG 1

/*
 This source code is originally from ISTI's ISCL and is distrubted under the
 Apache 2 license.  It has been heavily modified to conform to C++.
*/

using namespace RTSeis::FilterRepresentations;
using namespace RTSeis::Utilities::FilterDesign;
using namespace RTSeis::Utilities::Math;

static int getNextWorstRealPoleIndex(const std::vector<bool> &pmask,
                                     const std::vector<bool> &isreal_pole,
                                     const std::vector<std::complex<double>> &p);
static int getWorstPoleIndex(const std::vector<bool> &pmask,
                             const std::vector<std::complex<double>> &p);
static int isreal_p_sum(const std::vector<bool> &pmask,
                        const std::vector<bool> &isreal_pole);
static int nearest_real_complex_idx(const std::vector<bool> &zmask,
                                    const std::vector<bool> &isreal_zero,
                                    const std::vector<std::complex<double>> &z,
                                    const std::complex<double> p1,
                                    const bool lreal);
static int nearest_realOrComplex_idx(const std::vector<bool> &zmask,
                                     const std::vector<std::complex<double>> &z,
                                     const std::complex<double> p1);
static int cmplxreal(const std::vector<std::complex<double>> &z,
                     std::vector<std::complex<double>> &p,
                     std::vector<bool> &isreal,
                     const double *tolIn = nullptr);

SOS IIR::designSOSIIRFilter(const int n, const double *W,
                            const double rp, const double rs,
                            const Bandtype btype,
                            const IIRPrototype ftype,
                            const IIRFilterDomain ldigital,
                            const SOSPairing pairing)
{
     RTSeis::FilterRepresentations::ZPK zpk;
     zpk = IIR::designZPKIIRFilter(n, W, rp, rs, btype, ftype, ldigital);
     SOS sos = IIR::zpk2sos(zpk, pairing);
     return sos;
}


//RTSeis::FilterRepresentations::BA
BA IIR::designBAIIRFilter(const int n, const double *W, 
                         const double rp, const double rs, 
                         const Bandtype btype,
                       const IIRPrototype ftype,
                       const IIRFilterDomain ldigital)
{
     RTSeis::FilterRepresentations::ZPK zpk;
     zpk = IIR::designZPKIIRFilter(n, W, rp, rs, btype, ftype, ldigital);
     BA ba = IIR::zpk2tf(zpk);
     return ba;
}

ZPK IIR::designZPKIIRFilter(const int n, const double *W,
                            const double rp, const double rs,
                            const Bandtype btype,
                            const IIRPrototype ftype,
                            const IIRFilterDomain ldigital)
{
    ZPK zpk;
    if (n < 1){RTSEIS_THROW_IA("Order = %d must be positive", n);}
    if (W == nullptr){RTSEIS_THROW_IA("%s", "W is NULL");}
    // Check ripples
    if (ftype == IIRPrototype::CHEBYSHEV1 && rp <= 0)
    {
        RTSEIS_THROW_IA("rp = %lf must be positive", rp);
    }
    if (ftype == IIRPrototype::CHEBYSHEV2 && rs <= 0)
    {
        RTSEIS_THROW_IA("rs = %lf must be positive", rs);
    }
    if (ldigital == IIRFilterDomain::ANALOG)
    {
        RTSEIS_WARNMSG("%s", "Analog filter design is untested");
    }
    // Check the low-corner frequency
    double wn[2];
    wn[0] = W[0];
    wn[1] = 0;
    if (wn[0] < 0){RTSEIS_THROW_IA("W[0] = %lf cannot be negative", wn[0]);}
    if (ldigital == IIRFilterDomain::DIGITAL && wn[0] > 1)
    {
        RTSEIS_THROW_IA("W[0]=%lf cannot exceed 1", wn[0]);
    }
    // Check the high-corner frequency
    if (btype == Bandtype::BANDPASS || btype == Bandtype::BANDSTOP)
    {
        wn[1] = W[1];
        if (wn[1] < wn[0])
        {
            RTSEIS_THROW_IA("W[1] = %lf can't be less than W[0] = %lf",
                             wn[0], wn[1]);
        }
        if (wn[1] > 1)
        {
            RTSEIS_THROW_IA("W[1] = %lf cannot exced 1", wn[1]);
        } 
    }
    // Create the analog prototype
    RTSeis::FilterRepresentations::ZPK zpkAp;
    if (ftype == IIRPrototype::BUTTERWORTH)
    {
        zpkAp =  AnalogPrototype::butter(n); // throws invalid argument
    }
    else if (ftype == IIRPrototype::BESSEL)
    {
        zpkAp = AnalogPrototype::bessel(n); // throws invalid argument
    }
    else if (ftype == IIRPrototype::CHEBYSHEV1)
    {
        zpkAp = AnalogPrototype::cheb1ap(n, rp); // throws invalid argument
    }
    else
    {
        zpkAp = AnalogPrototype::cheb2ap(n, rs); // throws invalid argument
    }
    // Pre-warp the frequencies
    double warped[2];
    double fs;
    if (ldigital == IIRFilterDomain::DIGITAL)
    {
        fs = 2.0;
        warped[0] = 2.0*fs*tan(M_PI*wn[0]/fs);
        warped[1] = 2.0*fs*tan(M_PI*wn[1]/fs);
    }
    else
    {
        fs = 1.0;
        warped[0] = wn[0];
        warped[1] = wn[1];
        if (wn[1] > 1.0)
        {
            RTSEIS_WARNMSG("%s", "wn[1] > 0; result may be unstable");
        }
    }
    // Transform to lowpass, bandpass, highpass, or bandstop 
    RTSeis::FilterRepresentations::ZPK zpktf;
    if (btype == Bandtype::LOWPASS)
    {
        zpktf = IIR::zpklp2lp(zpkAp, warped[0]); // throws
    }
    else if (btype == Bandtype::HIGHPASS)
    {
        zpktf = IIR::zpklp2hp(zpkAp, warped[0]); // throws
    }
    else if (btype == Bandtype::BANDPASS)
    {
        double bw = warped[1] - warped[0];
        double w0 = std::sqrt(warped[0]*warped[1]);
        zpktf = IIR::zpklp2bp(zpkAp, w0, bw); // throws
    }
    else if (btype == Bandtype::BANDSTOP)
    {
        double bw = warped[1] - warped[0];
        double w0 = std::sqrt(warped[0]*warped[1]);
        zpktf = IIR::zpklp2bs(zpkAp, w0, bw); // throws
    }
    // Find the discrete equivalent
    if (ldigital == IIRFilterDomain::DIGITAL)
    {
        zpk = zpkbilinear(zpktf, fs);
    }
    else
    {
        zpk = zpktf;
    }
    return zpk;
}

SOS IIR::zpk2sos(const RTSeis::FilterRepresentations::ZPK &zpk,
                 const SOSPairing pairing)
{
    SOS sos;
    ZPK zpkWork = zpk;
    double k = zpkWork.getGain();
    int npoles = zpkWork.getNumberOfPoles();
    int nzeros = zpkWork.getNumberOfZeros();
    if (npoles <= 0 || nzeros <= 0)
    {
        RTSEIS_WARNMSG("sos will simply scale data by %lf", k); 
        std::vector<double> bs({k, 0, 0});
        std::vector<double> as({1, 0, 0});
        sos = RTSeis::FilterRepresentations::SOS(1, bs, as);
        return sos;
    }
    // Get sizes and handles to zeros and poles
    int np = npoles + std::max(nzeros - npoles, 0);
    int nz = nzeros + std::max(npoles - nzeros, 0);
    int nSections = std::max(np, nz + 1)/2;
#ifdef DEBUG
    assert(nSections > 0);
#endif
    std::vector<std::complex<double>> zin = zpkWork.getZeros();
    std::vector<std::complex<double>> pin = zpkWork.getPoles();
    // Base case
    if (nSections == 1)
    {
        // (z - z_1)*(z - 0) = z^2 - z*z_1 + 0
        std::vector<double> bs({1, 0, 0});
        if (nzeros == 1)
        {
            RTSEIS_WARNMSG("%s", "nzeros = 1 not yet tested");
            if (std::abs(std::imag(zin[0])) > 1.e-15)
            {
                RTSEIS_WARNMSG("%s", "Taking real part of zin");
            }
            bs[1] =-std::real(zin[0]);
        }
        else if (nzeros == 2)
        {
            bs[1] =-2.0*std::real(zin[0]);
            bs[2] = std::pow(std::abs(zin[0]), 2);
        }
        // Introduce gain
        bs[0] = bs[0]*k;
        bs[1] = bs[1]*k;
        bs[2] = bs[2]*k;
        // (p - p_1)*(p - 0) = p^2 - p*p_1 + 0
        std::vector<double> as({1, 0, 0});
        if (npoles == 1)
        {
            RTSEIS_WARNMSG("%s", "npoles = 1 not yet tested");
            if (std::abs(std::imag(pin[0])) > 1.e-15)
            {
                RTSEIS_WARNMSG("%s", "Taking real part of pin");
            }
            as[1] =-std::real(pin[0]);
        }
        // Conjugate pairs: (p - p_1)*(p - conj(p_1)) 
        else if (npoles == 2)
        {
            as[1] =-2.0*std::real(pin[0]);
            as[2] = std::pow(std::abs(pin[0]), 2);
        }
        sos = RTSeis::FilterRepresentations::SOS(1, bs, as);
        //RTSEIS_WARNMSG("%s", "Summary of design");
        //sos.print(stdout);
        return sos;
    }
    if (np%2 == 1 && pairing == SOSPairing::NEAREST)
    {
        np = np + 1;
        nz = nz + 1;
    }
    if (np != nz || np < 1 || nz < 1)
    {
        if (np != nz){RTSEIS_ERRMSG("Error size inconsistent %d,%d!", np, nz);}
        if (np < 1){RTSEIS_ERRMSG("Error no zeros in use %d", nz);}
        if (nz < 1){RTSEIS_ERRMSG("Error no poles in use %d", np);}
        RTSEIS_THROW_IA("%s", "Invalid inputs");
    }
    std::vector<std::complex<double>> z;
    std::vector<bool> isreal_zero;
#ifdef DEBUG
    assert(cmplxreal(zin, z, isreal_zero, nullptr) == 0);
#else
    cmplxreal(zin, z, isreal_zero, nullptr);
#endif
    std::vector<std::complex<double>> p;
    std::vector<bool> isreal_pole;
#ifdef DEBUG
    assert(cmplxreal(pin, p, isreal_pole, nullptr) == 0);
#else
    cmplxreal(pin, p, isreal_pole, nullptr); 
#endif
    std::vector<bool> zmask(z.size(), false);
    std::vector<bool> pmask(p.size(), false);
    std::vector<BA> basAll;
    for (int is=0; is<nSections; is++)
    {
        std::complex<double> z1 = std::complex<double> (0, 0);
        // Select the `worst' pole
        int p1_idx = getWorstPoleIndex(pmask, p);
#ifdef DEBUG
        assert(p1_idx !=-1);
#endif
        std::complex<double> p1 = p[p1_idx];
        pmask[p1_idx] = true;
        int psum = isreal_p_sum(pmask, isreal_pole);
        // Initialize complex conjugates for real case
        std::complex<double> p2 = std::complex<double> (0, 0);
        std::complex<double> z2 = std::complex<double> (0, 0);
        // Pair that pole with a zero
        if (isreal_pole[p1_idx] && psum == 0)
        {
            // Special case to set a first order section
            int z1_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                  z, p1, true);
#ifdef DEBUG
            assert(z1_idx !=-1);
#endif
            z1 = z[z1_idx];
            zmask[z1_idx] = true;
        }
        else
        {
            int zsum = isreal_p_sum(zmask, isreal_zero);
            int z1_idx =-1;
            if (!isreal_pole[p1_idx] && zsum == 1)
            {
                // Special case to ensure choose a complex zero to pair
                // with so later (setting up a first-order section)
                z1_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                  z, p1, false);
#ifdef DEBUG
                assert(z1_idx !=-1);
                assert(isreal_zero[z1_idx]);
#endif
            }
            else
            {
                // Pair that pole with the closest zero (real or complex)
                z1_idx = nearest_realOrComplex_idx(zmask, z, p1);
#ifdef DEBUG
                assert(z1_idx !=-1);
#endif
            }
            z1 = z[z1_idx];
            zmask[z1_idx] = true;
            // Now that we have p1 and z1, figure out what p2 and z2 need to be
            if (!isreal_pole[p1_idx])
            {
                // Complex pole, complex zero -> complex conjugates
                if (!isreal_zero[z1_idx])
                {
                    p2 = std::conj(p1);
                    z2 = std::conj(z1);
                }
                // Complex pole and real zero
                else
                {
                    p2 = std::conj(p1);
                    int z2_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                          z, p1, true);
#ifdef DEBUG
                    assert(z2_idx !=-1);
                    assert(isreal_zero[z2_idx]);
#endif
                    z2 = z[z2_idx];
                    zmask[z2_idx] = true;
                }
            }
            // Real pole
            else
            {
                // Real pole, complex zero
                if (!isreal_zero[z1_idx])
                {
                    z2 = conj(z1);
                    int p2_idx = nearest_real_complex_idx(pmask, isreal_pole,
                                                          p, z1, true);
#ifdef DEBUG
                    assert(p2_idx !=-1);
                    assert(isreal_pole[p2_idx]);
#endif
                    p2 = p[p2_idx];
                    pmask[p2_idx] = true;
                }
                // Real pole, real zero
                else
                {
                    // Pick the next `worst' pole to use
                    int p2_idx = getNextWorstRealPoleIndex(pmask,
                                                           isreal_pole, p);
#ifdef DEBUG
                    assert(p2_idx !=-1);
                    assert(isreal_pole[p2_idx]);
#endif
                    p2 = p[p2_idx];
                    // Find a real zero to match the added pole
                    int z2_idx = nearest_real_complex_idx(zmask, isreal_zero,
                                                          z, p2, true);
#ifdef DEBUG
                    assert(z2_idx !=-1);
                    assert(isreal_zero[z2_idx]);
#endif
                    z2 = z[z2_idx];
                    zmask[z2_idx] = true; 
                    pmask[p2_idx] = true;
                }
            }
        }
        std::vector<std::complex<double>> ptemp({p1, p2});
        std::vector<std::complex<double>> ztemp({z1, z2});
        double ktemp = 1;
        if (is == nSections - 1){ktemp = k;}
        RTSeis::FilterRepresentations::ZPK zpkTemp
             = RTSeis::FilterRepresentations::ZPK(ztemp, ptemp, ktemp);
        BA baTemp;
        baTemp = zpk2tf(zpkTemp);
        // Save the transfer function of the secon dorder section
        basAll.push_back(baTemp);
    } // Loop on sections
    // Reality check
    for (size_t ip=0; ip<pmask.size(); ip++)
    {
        if (!pmask[ip])
        {
            RTSEIS_ERRMSG("Failed to find pole %ld; nsections=%d",
                          ip, nSections);
#ifdef DEBUG
            assert(false);
#endif
            return sos;
        }
    }
    for (size_t iz=0; iz<zmask.size(); iz++)
    {
        if (!zmask[iz])
        {
            RTSEIS_ERRMSG("Failed to find zero %ld", iz);
#ifdef DEBUG
            assert(false);
#endif
            return sos;
        }
    }
    // Construct the cascading transfer function with the `worst' pole last
    int js = nSections;
    std::vector<double> bs;
    std::vector<double> as;
    for (int is=0; is<nSections; is++)
    {
        js = js - 1;
        std::vector<double> bsTemp = basAll[js].getNumeratorCoefficients(); 
        std::vector<double> asTemp = basAll[js].getDenominatorCoefficients();
        for (size_t k=0; k<3; k++)
        {
            bs.push_back(bsTemp[k]);
            as.push_back(asTemp[k]);
        }
    }
    // Finally pack the second order sections
    sos = RTSeis::FilterRepresentations::SOS(nSections, bs, as);
    return sos;
}

ZPK IIR::tf2zpk(const RTSeis::FilterRepresentations::BA &ba)
{
    std::vector<double> b = ba.getNumeratorCoefficients();
    std::vector<double> a = ba.getDenominatorCoefficients();
    if (b.size() < 1)
    {
        RTSEIS_THROW_IA("%s", "No numerator coefficients");
    }
    if (a.size() < 1)
    { 
        RTSEIS_THROW_IA("%s", "No denominator coefficients");
    }
    double a0 = a[0];
    if (a0 == 0){RTSEIS_THROW_IA("%s", "a[0] = 0");}
    double b0 = b[0];
    if (b0 == 0){RTSEIS_THROW_IA("%s", "b[0] = 0");}
    // Normalize
    #pragma omp simd
    for (size_t i=0; i<b.size(); i++){b[i] = b[i]/a0;} 
    #pragma omp simd
    for (size_t i=0; i<a.size(); i++){a[i] = a[i]/a0;}
    // Compute gain
    double k = b[0]; 
    #pragma omp simd
    for (size_t i=0; i<b.size(); i++){b[i] = b[i]/k;}
    //  Compute the roots of the numerator
    std::vector<std::complex<double>> z; //(std::max(b.size(), ltwo) - lone);
    if (b.size() > 1)
    {
        z = Polynomial::roots(b);
    }
    else
    {
        RTSEIS_WARNMSG("%s", "Warning no zeros");
        z.resize(1);
        z[0] = std::complex<double> (1, 0);
    }
    // Compute the roots of the denominator
    std::vector<std::complex<double>> p; //(std::max(a.size(), ltwo) - lone);
    if (a.size() > 1)
    {
        p = Polynomial::roots(a);
    }
    else
    {
        RTSEIS_WARNMSG("%s", "No poles");
        p.resize(1);
        p[0] = std::complex<double> (1, 0); 
    }
    // Make a transfer function
    ZPK zpk(z, p, k);
    return zpk;
}
                 
BA IIR::zpk2tf(const RTSeis::FilterRepresentations::ZPK &zpk) noexcept
{
    BA ba;
    double k = zpk.getGain();
    if (k == 0){RTSEIS_WARNMSG("%s", "System gain is zero");}
    // Compute polynomial representation of zeros by expanding:
    //  (z - z_1)*(z - z_2)*...*(z - z_nzeros)
    std::vector<std::complex<double>> z = zpk.getZeros();
    std::vector<std::complex<double>> bz;
    bz = Polynomial::poly(z) ;//, bz); 
    // Compute polynomial representation of zeros by expanding:
    //  (z - z_1)*(z - z_2)*...*(z - z_nzeros)
    std::vector<std::complex<double>> p = zpk.getPoles();
    std::vector<std::complex<double>> az;
    az = Polynomial::poly(p); //, az);
    // Introduce gain into the numerator zeros
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<bz.size(); i++){bz[i] = k*bz[i];}
    // Take the real
    std::vector<double> b;
    b.resize(bz.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<bz.size(); i++){b[i] = std::real(bz[i]);}
    std::vector<double> a;
    a.resize(az.size());
#ifdef __INTEL_COMPILER
    #pragma ivdep
#endif
    for (size_t i=0; i<az.size(); i++){a[i] = std::real(az[i]);}
    // Pack into a transfer function struct
    try
    {
        ba = RTSeis::FilterRepresentations::BA(b, a);
    }
    catch (const std::invalid_argument &ia)
    {
        RTSEIS_ERRMSG("The impossible has happened %s", ia.what());
    }
    return ba;
}

ZPK IIR::zpkbilinear(const RTSeis::FilterRepresentations::ZPK zpk,
                     const double fs)
{
    size_t nzeros = zpk.getNumberOfZeros();
    size_t npoles = zpk.getNumberOfPoles();
    if (nzeros > npoles)
    {
        RTSEIS_THROW_IA("%s", "Cannot have more zeros than poles");
    }
    std::vector<std::complex<double>> z = zpk.getZeros();
    std::vector<std::complex<double>> p = zpk.getPoles();
    // Bilinear transform
    std::complex<double> fs2 = std::complex<double> (2*fs, 0);
    std::vector<std::complex<double>> z_bl;
    z_bl.resize(npoles);
    for (size_t i=0; i<nzeros; i++)
    {
        std::complex<double> znum = fs2 + z[i];
        std::complex<double> zden = fs2 - z[i];
        z_bl[i] = znum/zden; //(fs2 + z)/(fs2 - z)
    }
    // Any zeros at infinity get moved to Nyquist frequency
    for (size_t i=nzeros; i<npoles; i++)
    {
        z_bl[i] = std::complex<double> (-1, 0);
    }
    std::vector<std::complex<double>> p_bl;
    p_bl.resize(npoles);
    for (size_t i=0; i<npoles; i++)
    {
        std::complex<double> znum = fs2 + p[i];
        std::complex<double> zden = fs2 - p[i];
        p_bl[i] = znum/zden; //(fs2 + p)/(fs2 - p) 
    }
    // Compensate for gain change
    std::complex<double> znum = std::complex<double> (1, 0);
    for (size_t i=0; i<nzeros; i++)
    {
        std::complex<double> zdif = fs2 - z[i];
        znum = znum*zdif;
    }
    std::complex<double> zden = std::complex<double> (1, 0);
    for (size_t i=0; i<npoles; i++)
    {
        std::complex<double> zdif = fs2 - p[i];
        zden = zden*zdif;
    }
    std::complex<double> zk = znum/zden;
    double k = zpk.getGain();
    double k_bl = k*std::real(zk);
    RTSeis::FilterRepresentations::ZPK zpkbl(z_bl, p_bl, k_bl);
    return zpkbl;
}

ZPK IIR::zpklp2bp(const RTSeis::FilterRepresentations::ZPK &zpkIn,
                  const double w0,
                  const double bw)
{
    size_t nzeros = zpkIn.getNumberOfZeros();
    size_t npoles = zpkIn.getNumberOfPoles();
    if (w0 < 0 || bw < 0)
    {
        if (w0 < 0){RTSEIS_THROW_IA("w0=%lf must be non-negative", w0);}
        if (bw < 0){RTSEIS_THROW_IA("bw=%lf must be non-negative", bw);}
        RTSEIS_THROW_IA("%s", "Invalid inputs");
    }
    size_t nzeros_bp = nzeros + npoles;
    size_t npoles_bp = 2*npoles;
    // Define some constants
    std::complex<double> zw02 = std::complex<double> (w0*w0 , 0);
    std::complex<double> zbw2 = std::complex<double> (bw/2.0, 0);
    // Duplicate and zeros and shift from baseband to +w0 and -w0
    std::vector<std::complex<double>> z_bp;
    z_bp.resize(nzeros_bp);
    std::vector<std::complex<double>> z = zpkIn.getZeros();
    size_t j = nzeros;
    for (size_t i=0; i<nzeros; i++)
    {
        std::complex<double> zlp = zbw2*z[i];        //Scale zeros to desired bw
        std::complex<double> zlp2 = zlp*zlp;         //zlp^2
        std::complex<double> zdif = zlp2 - zw02;     //zlp^2 - w0^2
        std::complex<double> zsqrt = std::sqrt(zdif);//sqrt(zlp^2 - w0^2)
        z_bp[i] = zlp + zsqrt; //zlp+sqrt(zlp^2-w0^2)
        z_bp[j] = zlp - zsqrt; //zlp-sqrt(zlp^2-w0^2)
        j = j + 1;
    }
    // Move degree zeros to origin leaving degree zeros at infinity for BPF 
    std::complex<double> zero(0, 0);
    for (size_t i=2*nzeros; i<nzeros_bp; i++)
    {
        z_bp[i] = zero;
    }
    // Duplicate poles and shift from baseband to +w0 and -w0
    std::vector<std::complex<double>> p_bp;
    p_bp.resize(npoles_bp);
    std::vector<std::complex<double>> p = zpkIn.getPoles();
    j = npoles;
    for (size_t i=0; i<npoles; i++)
    {   
        std::complex<double> plp = zbw2*p[i];        //Scale poles to desired bw
        std::complex<double> plp2 = plp*plp;         //plp^2
        std::complex<double> pdif = plp2 - zw02;     //plp^2 - w0^2
        std::complex<double> psqrt = std::sqrt(pdif);//sqrt(plp^2 - w0^2)
        p_bp[i] = plp + psqrt; //plp+sqrt(plp^2-w0^2)
        p_bp[j] = plp - psqrt; //plp-sqrt(plp^2-w0^2)
        j = j + 1;
    }
    // Cancel out gain change from frequency scaling
    double k = zpkIn.getGain();
    double k_bp = k*std::pow(bw,npoles-nzeros); //k*bw*bw^degree
    RTSeis::FilterRepresentations::ZPK zpkOut(z_bp, p_bp, k_bp);
    return zpkOut;
}

ZPK IIR::zpklp2bs(const RTSeis::FilterRepresentations::ZPK &zpkIn,
                  const double w0,
                  const double bw)
{
    size_t nzeros = zpkIn.getNumberOfZeros();
    size_t npoles = zpkIn.getNumberOfPoles();
    if (w0 < 0 || bw < 0)
    {
        if (w0 < 0){RTSEIS_THROW_IA("w0=%lf must be non-negative", w0);}
        if (bw < 0){RTSEIS_THROW_IA("bw=%lf must be non-negative", bw);}
        RTSEIS_THROW_IA("%s", "Invalid inputs");
    }
    size_t nzeros_bs = 2*npoles;
    size_t npoles_bs = 2*npoles;
    // Define some constants
    std::complex<double> zw02 = std::complex<double> (w0*w0 , 0);
    std::complex<double> zbw2 = std::complex<double> (bw/2.0, 0);
    // Duplicate the zeros and shift from baseband to +w0 and -w0
    std::vector<std::complex<double>> z = zpkIn.getZeros();
    std::vector<std::complex<double>> z_bs;
    z_bs.resize(nzeros_bs);
    std::complex<double> zprod = std::complex<double> (1, 0); // Init product
    size_t j = nzeros;
    for (size_t i=0; i<nzeros; i++)
    {
        // Invert to a highpass filter with desired bandwidth
        std::complex<double> zhp = zbw2/z[i];         // bw2/2/z
        std::complex<double> zhp2 = zhp*zhp;          // zhp^2
        std::complex<double> zdif = zhp2 - zw02;      // zhp^2 - w0^2
        std::complex<double> zsqrt = std::sqrt(zdif); // sqrt(zhp^2 - w0^2)
        z_bs[i] = zhp + zsqrt; // zhp+sqrt(zhp^2-w0^2)
        z_bs[j] = zhp - zsqrt; // zhp-sqrt(zhp^2-w0^2)
        // Update product 
        zprod = (-z[i])*zprod;
        j = j + 1;
    }
    // Move any zeros at infinity to center of stopband
    j = 2*nzeros;
    for (size_t i=0; i < npoles - nzeros; i++)
    {
        z_bs[j] = std::complex<double> (0, w0);
        j = j + 1;
    }
    for (size_t i=0; i < npoles - nzeros; i++)
    {
        z_bs[j] = std::complex<double> (0,-w0);
        j = j + 1;
    }
    // Duplicate the poles and shift from baseband to +w0 and -w0
    std::vector<std::complex<double>> p = zpkIn.getPoles();
    std::complex<double> pprod = std::complex<double> (1, 0);
    std::vector<std::complex<double>> p_bs;
    p_bs.resize(npoles_bs);
    j = npoles;
    for (size_t i=0; i<npoles; i++)
    {
        // Invert to a highpass filter with desired bandwidth
        std::complex<double> php = zbw2/p[i];         // bw2/2/p
        std::complex<double> php2 = php*php;          // php^2
        std::complex<double> pdif = php2 - zw02;      // php^2 - w0^2
        std::complex<double> psqrt = std::sqrt(pdif); // sqrt(php^2 - w0^2)
        p_bs[i] = php + psqrt; // php+sqrt(php^2-w0^2)
        p_bs[j] = php - psqrt; // php-sqrt(php^2-w0^2)
        // Update product
        pprod = (-p[i])*pprod;
        j = j + 1;
    }
    // Cancel out gain change caused by inversion
    double k = zpkIn.getGain();
    std::complex<double> zdiv =  zprod/pprod;
    double k_bs = k*std::real(zdiv);
    RTSeis::FilterRepresentations::ZPK zpkOut(z_bs, p_bs, k_bs);
    return zpkOut;
}

ZPK IIR::zpklp2lp(const RTSeis::FilterRepresentations::ZPK &zpkIn,
                  const double w0)
{
    size_t nzeros = zpkIn.getNumberOfZeros();
    size_t npoles = zpkIn.getNumberOfPoles();
    if (w0 < 0){RTSEIS_THROW_IA("w0=%lf must be non-negative", w0);}
    if (npoles < nzeros)
    {
        RTSEIS_THROW_IA("%s", "BUG! npoles < nzeros not implemented");
    }
    std::vector<std::complex<double>> z = zpkIn.getZeros();
    std::vector<std::complex<double>> p = zpkIn.getPoles();
    // Scale all points radially from origin to shift cutoff frequency
    std::vector<std::complex<double>> z_lp(nzeros);
    for (size_t i=0; i<z.size(); i++)
    {
        z_lp[i] = w0*z[i];
    }
    std::vector<std::complex<double>> p_lp(npoles);
    for (size_t i=0; i<p.size(); i++)
    {
        p_lp[i] = w0*p[i];
    }
    // Each shifted pole decreases gain by w0, each shifted zero increases
    // it.  Cancel out the net change to keep overall gain the same
    double k = zpkIn.getGain(); 
    int ndeg = npoles - nzeros;
    double k_lp = k*std::pow(w0, ndeg);
    RTSeis::FilterRepresentations::ZPK zpkOut(z_lp, p_lp, k_lp);
    return zpkOut;
}

ZPK IIR::zpklp2hp(const RTSeis::FilterRepresentations::ZPK &zpkIn,
                  const double w0)
{
    size_t nzeros = zpkIn.getNumberOfZeros();
    size_t npoles = zpkIn.getNumberOfPoles();
    if (w0 < 0){RTSEIS_THROW_IA("w0=%lf must be non-negative", w0);}
    if (npoles < nzeros)
    {
        RTSEIS_THROW_IA("%s", "BUG! npoles < nzeros not implemented");
    }
    std::vector<std::complex<double>> z = zpkIn.getZeros();
    std::vector<std::complex<double>> p = zpkIn.getPoles();
    // Invert positions radially about unit circle to convert LPF TO HPF
    // Scale all points radially from origin to shift cutoff frequency
    std::vector<std::complex<double>> z_hp(std::max(nzeros, npoles));
    std::complex<double> zprod = std::complex<double> (1, 0);
    for (size_t i=0; i<z.size(); i++)
    {
        z_hp[i] = w0/z[i];
        zprod = (-z[i])*zprod;
    }
    // If lowpass had zeros at infinity, inverting moves them to origin
    for (size_t i=nzeros; i<npoles; i++)
    {
        z_hp[i] = std::complex<double> (0, 0);
    }
    // Invert positions radially about unit circle to convert LPF TO HPF
    // Scale all points radially from origin to shift cutoff frequency
    std::vector<std::complex<double>> p_hp(npoles);
    std::complex<double> pprod = std::complex<double> (1, 0);
    for (size_t i=0; i<p.size(); i++)
    {
        p_hp[i] = w0/p[i];
        pprod = (-p[i])*pprod;
    }
    // Compute scale factors and cancel out gain caused by inversion
    double k = zpkIn.getGain();
    std::complex<double> zt = zprod/pprod;
    double k_hp = k*std::real(zt);
    RTSeis::FilterRepresentations::ZPK zpkOut(z_hp, p_hp, k_hp);
    return zpkOut;
}
//============================================================================//
//                            Utility functions for zpk2sos                   //
//============================================================================//
/*!
 * @brief Utility function for getting index of next closest pole to
 *        unit circle.
 */
static int getNextWorstRealPoleIndex(const std::vector<bool> &pmask,
                                     const std::vector<bool> &isreal_pole,
                                     const std::vector<std::complex<double>> &p)
{
    const std::complex<double> zone = std::complex<double> (1, 0);
    double difMin = DBL_MAX;
    int indx =-1;
    for (size_t i=0; i<p.size(); i++)
    {
        if (!pmask[i] && isreal_pole[i])
        {
            double dif = std::abs(p[i] - zone);
            if (dif < difMin)
            {
                dif = difMin;
                indx = i;
            }
        }
    }
    return indx;
}
/*!
 * @brief Utility function for getting next pole closest to unit circle.
 */
static int getWorstPoleIndex(const std::vector<bool> &pmask,
                             const std::vector<std::complex<double>> &p)
{
    double difMin = DBL_MAX;
    int indx =-1;
    for (size_t i=0; i<pmask.size(); i++)
    {
        if (!pmask[i])
        {
            double dif = std::abs(1.0 - std::abs(p[i]));
            if (dif < difMin)
            {
                difMin = dif;
                indx = i;
            }
        }
    }
    return indx;
}
/*!
 * @brief Utility function for computing the sum of the reals of the poles
 *        array
 */
static int isreal_p_sum(const std::vector<bool> &pmask,
                        const std::vector<bool> &isreal_pole)
{
    int xsum = 0;
    for (size_t i=0; i<isreal_pole.size(); i++)
    {
        if (!pmask[i] && isreal_pole[i]){xsum = xsum + 1;}
    }
    return xsum;
}
/*!
 * @brief Utility function for finding the closest real or complex zero
 *        to the pole p1
 */
static int nearest_real_complex_idx(const std::vector<bool> &zmask,
                                    const std::vector<bool> &isreal_zero,
                                    const std::vector<std::complex<double>> &z,
                                    const std::complex<double> p1,
                                    const bool lreal)
{
    double difMin = DBL_MAX;
    int indx =-1;
    // Look for the closest real zero to the pole p1
    if (lreal)
    {
        for (size_t i=0; i<z.size(); i++)
        {
            if (isreal_zero[i] && !zmask[i])
            {
                double dif = std::abs(z[i] - p1);
                if (dif < difMin)
                {
                    difMin = dif;
                    indx = i;
                }
            }
        }
    }
    // Look for the closest complex zero to the pole p1
    else
    {
        for (size_t i=0; i<z.size(); i++)
        {
            if (!isreal_zero[i] && !zmask[i])
            {
                double dif = std::abs(z[i] - p1);
                if (dif < difMin)
                {
                    difMin = dif;
                    indx = i;
                }
            }
        }
    }
    return indx;
}

/*!
 * @brief Find the index of the next nearest real or complex zero to the
 *        given pole.
 */
static int nearest_realOrComplex_idx(const std::vector<bool> &zmask,
                                     const std::vector<std::complex<double>> &z,
                                     const std::complex<double> p1)
{
    double difMin = DBL_MAX;
    int indx =-1;
    for (size_t i=0; i<z.size(); i++)
    {
        if (!zmask[i])
        {
            double dif = std::abs(p1 - z[i]);
            if (dif < difMin)
            {
                difMin = dif;
                indx = i;
            }
        }
    }
    return indx;
}
static int cmplxreal(const std::vector<std::complex<double>> &z,
                     std::vector<std::complex<double>> &p,
                     std::vector<bool> &isreal,
                     const double *tolIn)
{
    p.resize(0);
    isreal.resize(0);
    if (z.size() < 1)
    {
        RTSEIS_ERRMSG("%s", "No elements in z");
        return -1;
    }
    // Set the tolerance
    double tol = 100.0*DBL_EPSILON;
    if (tolIn != nullptr){tol = *tolIn;}
    // Get the real poles and complex poles 
    std::vector<double> reals;
    std::vector<std::complex<double>> cmplxs;
    for (size_t i=0; i<z.size(); i++)
    {
        if (std::imag(z[i]) == 0)
        {
            reals.push_back(std::real(z[i]));
        }
        else
        {
            cmplxs.push_back(z[i]);
        }
    }
    size_t nreal  = reals.size();
    size_t ncmplx = cmplxs.size();
    if (nreal + ncmplx != z.size())
    {
        RTSEIS_ERRMSG("%s", "Algorithmic error");
        return -1;
    }
    // Sort the reals into ascending order and save them 
    nreal = 0;
    if (reals.size() > 0)
    {
        std::sort(reals.begin(), reals.end());
        for (size_t i=0; i<reals.size(); i++)
        {
            p.push_back(std::complex<double> (reals[i], 0));
            isreal.push_back(true);
            nreal = nreal + 1;
        }
    }
    // Pair the complex conjugate poles up and save them
    size_t nconj = 0;
    if (cmplxs.size() > 0)
    {
        // Sort complex numbers based on magnitudes
        std::sort(cmplxs.begin(), cmplxs.end(),
                  [](std::complex<double> a,  
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
        // Verify everyone has a buddy (conjugate pair)
        std::vector<int> lskip(cmplxs.size(), 0);
        for (size_t i=0; i<cmplxs.size(); i++)
        {
            if (lskip[i] == 1){continue;}
            // Find it's complex conjugate
            for (size_t j=i+1; j<cmplxs.size(); j++)
            {
                if (lskip[j] == 1){continue;}
                if (std::abs(cmplxs[i] - std::conj(cmplxs[j])) < tol)
                {
                    // Balance poles for numerical accuracy
                    std::complex<double> zhalf
                         = (cmplxs[i] + std::conj(cmplxs[j]))/2.0;
                    // Retain the `positive' pole for consistency 
                    if (std::imag(zhalf) >= 0)
                    {
                        p.push_back(zhalf); 
                        //p.push_back(std::conj(zhalf));
                    }
                    else
                    {
                        p.push_back(std::conj(zhalf)); 
                        //p.push_back(zhalf);
                    }
                    nconj = nconj + 1;
                    isreal.push_back(false);
                    //isreal.push_back(false);
                    lskip[i] = 1;
                    lskip[j] = 1;
                    break;
                }
            }
        }
        // I'm not dealing with unpaired
        int npaired = std::accumulate(lskip.begin(), lskip.end(), 0);
        if (npaired != static_cast<int> (cmplxs.size()))
        {
            RTSEIS_ERRMSG("All poles need complex conjugate pairs; %d %ld",
                          npaired, cmplxs.size());
            fprintf(stderr, "Input\n");
            for (size_t m=0; m<z.size(); m++)
            {
                fprintf(stderr, "%e + %+ei\n", z[m].real(), z[m].imag());
            }
            fprintf(stderr, "Found complexes\n");
            for (size_t m=0; m<cmplxs.size(); m++) 
            {
                fprintf(stderr, "%e + %+ei\n", cmplxs[m].real(), cmplxs[m].imag());
            }
            p.resize(0);
getchar();
            return -1;
        }
    }
    if (2*nconj + nreal != z.size())
    {
        RTSEIS_ERRMSG("%s", "Algorithmic failure");
        return -1;
    }
/*
printf("Input\n");
for (size_t i=0; i<z.size(); i++)
{
 printf("%+lf%+lfi\n", std::real(z[i]), std::imag(z[i]));
}
printf("Output\n");
for (size_t i=0; i<p.size(); i++)
{
 printf("%+lf%+lfi\n", std::real(p[i]), std::imag(p[i]));
}
*/
    return 0;
}
