#define RTSEIS_LOGGING 1
#include <cmath>
#include <algorithm>
#include "rtseis/utilities/filterRepresentations/zpk.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

class ZPK::ZPKImpl
{
public:
    /// Zeros
    std::vector<std::complex<double>> z;
    /// Poles
    std::vector<std::complex<double>> p;
    /// Gain
    double k = 0;
    /// Tolerance in checking equality.
    double tol = DEFAULT_TOL;

};

ZPK::ZPK(void) :
    pImpl(std::make_unique<ZPKImpl>())
{
    return;
}

ZPK::ZPK(const std::vector<std::complex<double>> &zeros,
         const std::vector<std::complex<double>> &poles,
         const double k) :
    pImpl(std::make_unique<ZPKImpl>())
{
    setZeros(zeros);
    setPoles(poles); 
    setGain(k);
    return;
}

ZPK::ZPK(const ZPK &zpk)
{
    *this = zpk;
    return;
}

ZPK::ZPK(ZPK &&zpk)
{
    *this = std::move(zpk);
    return;
}

ZPK::~ZPK(void) = default;

ZPK& ZPK::operator=(ZPK &&zpk)
{
    if (&zpk == this){return *this;}
    pImpl = std::move(zpk.pImpl);
    return *this;
}

ZPK& ZPK::operator=(const ZPK &zpk)
{
    if (&zpk == this){return *this;}
    pImpl = std::make_unique<ZPKImpl> ();
    //pImpl = std::unique_ptr<ZPKImpl> (new ZPKImpl());
    pImpl->z = zpk.pImpl->z;
    pImpl->p = zpk.pImpl->p;
    pImpl->k = zpk.pImpl->k;
    pImpl->tol = zpk.pImpl->tol;
    return *this;
}

bool ZPK::operator==(const ZPK &zpk) const
{
    if (pImpl->p.size() != zpk.pImpl->p.size()){return false;}
    if (pImpl->z.size() != zpk.pImpl->z.size()){return false;}
    for (size_t i=0; i<pImpl->p.size(); i++)
    {
        if (std::abs(pImpl->p[i] - zpk.pImpl->p[i]) > pImpl->tol)
        {
            return false;
        }
    }
    for (size_t i=0; i<pImpl->z.size(); i++)
    {
        if (std::abs(pImpl->z[i] - zpk.pImpl->z[i]) > pImpl->tol)
        {
            return false;
        }
    }
    if (std::abs(pImpl->k - zpk.pImpl->k) > pImpl->tol){return false;}
    return true;
}

bool ZPK::operator!=(const ZPK &zpk) const
{
    return !(*this == zpk); 
}

void ZPK::print(FILE *fout) const noexcept
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Gain: %.16lf\n", pImpl->k);
    fprintf(f, "Zeros:\n");
    for (size_t i=0; i<pImpl->z.size(); i++)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n",
                std::real(pImpl->z[i]), std::imag(pImpl->z[i])); 
    }
    fprintf(f, "Poles:\n");
    for (size_t i=0; i<pImpl->p.size(); i++)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n",
                std::real(pImpl->p[i]), std::imag(pImpl->p[i]));
    }
    return;
}

void ZPK::sortPoles(const bool ascending)
{
    if (ascending)
    {
        std::sort(pImpl->p.begin(), pImpl->p.end(),
                  // Begin lambda
                  [](std::complex<double> a,
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
    }
    else
    {
        std::sort(pImpl->p.begin(), pImpl->p.end(),
                  // Begin lambda
                  [](std::complex<double> a,
                     std::complex<double> b)
                  {
                     return (std::abs(a) > std::abs(b));
                  });
    }
    return; 
}

void ZPK::sortZeros(const bool ascending)
{
    if (ascending)
    {
        std::sort(pImpl->z.begin(), pImpl->z.end(),
                  // Begin lambda
                  [](std::complex<double> a,
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
    }   
    else
    {
        std::sort(pImpl->z.begin(), pImpl->z.end(),
                  // Begin lambda
                  [](std::complex<double> a,
                     std::complex<double> b)
                  {
                     return (std::abs(a) > std::abs(b));
                  });
    }   
    return; 
}

void ZPK::clear(void)
{
    pImpl->p.clear();
    pImpl->z.clear();
    pImpl->k = 0;
    pImpl->tol = DEFAULT_TOL;
    return;
}

void ZPK::setGain(const double k)
{
    if (k == 0){RTSEIS_WARNMSG("%s", "Gain is zero");}
    pImpl->k = k;
    return;
}

double ZPK::getGain(void) const
{
    return pImpl->k;
}

int ZPK::getNumberOfPoles(void) const
{
    return static_cast<int> (pImpl->p.size());
}

int ZPK::getNumberOfZeros(void) const
{
    return static_cast<int> (pImpl->z.size());
}

void ZPK::setPoles(const size_t n, std::complex<double> poles[])
{
    if (n > 0 && poles == nullptr)
    {
        RTSEIS_ERRMSG("%s", "Poles is null");
        pImpl->p.resize(0);
        return;
    }
    pImpl->p.resize(n);
    for (size_t i=0; i<n; i++)
    {
        pImpl->p[i] = poles[i];
    }
    return;
}

void ZPK::setPoles(const std::vector<std::complex<double>> &poles)
{
    pImpl->p = poles;
    return;
}

void ZPK::setZeros(const size_t n, std::complex<double> zeros[])
{
    if (n > 0 && zeros == nullptr)
    {
        RTSEIS_ERRMSG("%s", "Zeros is null");
        pImpl->z.resize(0);
        return;
    }
    pImpl->z.resize(n);
    for (size_t i=0; i<n; i++)
    {
        pImpl->z[i] = zeros[i];
    }
    return;
}

void ZPK::setZeros(const std::vector<std::complex<double>> &zeros)
{
    pImpl->z = zeros;
    return;
}

std::vector<std::complex<double>> ZPK::getPoles(void) const
{
    return pImpl->p;
}

std::vector<std::complex<double>> ZPK::getZeros(void) const
{
    return pImpl->z;
}

void ZPK::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    pImpl->tol = tol;
    return;
}
