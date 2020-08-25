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

ZPK::ZPK() :
    pImpl(std::make_unique<ZPKImpl>())
{
}

ZPK::ZPK(const std::vector<std::complex<double>> &zeros,
         const std::vector<std::complex<double>> &poles,
         const double k) :
    pImpl(std::make_unique<ZPKImpl>())
{
    setZeros(zeros);
    setPoles(poles); 
    setGain(k);
}

ZPK::ZPK(const ZPK &zpk)
{
    *this = zpk;
}

ZPK::ZPK(ZPK &&zpk) noexcept
{
    *this = std::move(zpk);
}

ZPK::~ZPK() = default;

ZPK& ZPK::operator=(ZPK &&zpk) noexcept
{
    if (&zpk == this){return *this;}
    pImpl = std::move(zpk.pImpl);
    return *this;
}

ZPK& ZPK::operator=(const ZPK &zpk)
{
    if (&zpk == this){return *this;}
    pImpl = std::make_unique<ZPKImpl> (*zpk.pImpl);
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
    for (auto & z : pImpl->z)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n",
                std::real(z), std::imag(z));
    }
    fprintf(f, "Poles:\n");
    for (auto & p : pImpl->p)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n",
                std::real(p), std::imag(p));
    }
}

[[maybe_unused]]
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
}

[[maybe_unused]]
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
}

void ZPK::clear() noexcept
{
    pImpl->p.clear();
    pImpl->z.clear();
    pImpl->k = 0;
    pImpl->tol = DEFAULT_TOL;
}

void ZPK::setGain(const double k)
{
    if (k == 0){RTSEIS_WARNMSG("%s", "Gain is zero");}
    pImpl->k = k;
}

double ZPK::getGain() const
{
    return pImpl->k;
}

int ZPK::getNumberOfPoles() const
{
    return static_cast<int> (pImpl->p.size());
}

int ZPK::getNumberOfZeros() const
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
}

void ZPK::setPoles(const std::vector<std::complex<double>> &poles)
{
    pImpl->p = poles;
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
}

void ZPK::setZeros(const std::vector<std::complex<double>> &zeros)
{
    pImpl->z = zeros;
}

std::vector<std::complex<double>> ZPK::getPoles() const
{
    return pImpl->p;
}

std::vector<std::complex<double>> ZPK::getZeros() const
{
    return pImpl->z;
}

[[maybe_unused]]
void ZPK::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    pImpl->tol = tol;
}
