#define RTSEIS_LOGGING 1
#include <cmath>
#include <algorithm>
#include "rtseis/utilities/filterRepresentations/zpk.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

struct ZPK::ZPKImpl
{
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
    pImpl_(new ZPKImpl())
{
    return;
}

ZPK::ZPK(const std::vector<std::complex<double>> &zeros,
         const std::vector<std::complex<double>> &poles,
         const double k) :
    pImpl_(new ZPKImpl())
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

ZPK::~ZPK(void)
{
    clear();
    return;
}

ZPK& ZPK::operator=(const ZPK &zpk)
{
    if (&zpk == this){return *this;}
    pImpl_ = std::unique_ptr<ZPKImpl> (new ZPKImpl());
    pImpl_->z = zpk.pImpl_->z;
    pImpl_->p = zpk.pImpl_->p;
    pImpl_->k = zpk.pImpl_->k;
    pImpl_->tol = zpk.pImpl_->tol;
    return *this;
}

bool ZPK::operator==(const ZPK &zpk) const
{
    if (pImpl_->p.size() != zpk.pImpl_->p.size()){return false;}
    if (pImpl_->z.size() != zpk.pImpl_->z.size()){return false;}
    for (size_t i=0; i<pImpl_->p.size(); i++)
    {
        if (std::abs(pImpl_->p[i] - zpk.pImpl_->p[i]) > pImpl_->tol)
        {
            return false;
        }
    }
    for (size_t i=0; i<pImpl_->z.size(); i++)
    {
        if (std::abs(pImpl_->z[i] - zpk.pImpl_->z[i]) > pImpl_->tol)
        {
            return false;
        }
    }
    if (std::abs(pImpl_->k - zpk.pImpl_->k) > pImpl_->tol){return false;}
    return true;
}

bool ZPK::operator!=(const ZPK &zpk) const
{
    return !(*this == zpk); 
}

void ZPK::print(FILE *fout)
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Gain: %.16lf\n", pImpl_->k);
    fprintf(f, "Zeros:\n");
    for (size_t i=0; i<pImpl_->z.size(); i++)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n",
                std::real(pImpl_->z[i]), std::imag(pImpl_->z[i])); 
    }
    fprintf(f, "Poles:\n");
    for (size_t i=0; i<pImpl_->p.size(); i++)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n",
                std::real(pImpl_->p[i]), std::imag(pImpl_->p[i]));
    }
    return;
}

void ZPK::sortPoles(const bool ascending)
{
    if (ascending)
    {
        std::sort(pImpl_->p.begin(), pImpl_->p.end(),
                  // Begin lambda
                  [](std::complex<double> a,
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
    }
    else
    {
        std::sort(pImpl_->p.begin(), pImpl_->p.end(),
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
        std::sort(pImpl_->z.begin(), pImpl_->z.end(),
                  // Begin lambda
                  [](std::complex<double> a,
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
    }   
    else
    {
        std::sort(pImpl_->z.begin(), pImpl_->z.end(),
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
    pImpl_->p.clear();
    pImpl_->z.clear();
    pImpl_->k = 0;
    pImpl_->tol = DEFAULT_TOL;
    return;
}

void ZPK::setGain(const double k)
{
    if (k == 0){RTSEIS_WARNMSG("%s", "Gain is zero");}
    pImpl_->k = k;
    return;
}

double ZPK::getGain(void) const
{
    return pImpl_->k;
}

int ZPK::getNumberOfPoles(void) const
{
    return static_cast<int> (pImpl_->p.size());
}

int ZPK::getNumberOfZeros(void) const
{
    return static_cast<int> (pImpl_->z.size());
}

void ZPK::setPoles(const size_t n, std::complex<double> poles[])
{
    if (n > 0 && poles == nullptr)
    {
        RTSEIS_ERRMSG("%s", "Poles is null");
        pImpl_->p.resize(0);
        return;
    }
    pImpl_->p.resize(n);
    for (size_t i=0; i<n; i++)
    {
        pImpl_->p[i] = poles[i];
    }
    return;
}

void ZPK::setPoles(const std::vector<std::complex<double>> &poles)
{
    pImpl_->p = poles;
    return;
}

void ZPK::setZeros(const size_t n, std::complex<double> zeros[])
{
    if (n > 0 && zeros == nullptr)
    {
        RTSEIS_ERRMSG("%s", "Zeros is null");
        pImpl_->z.resize(0);
        return;
    }
    pImpl_->z.resize(n);
    for (size_t i=0; i<n; i++)
    {
        pImpl_->z[i] = zeros[i];
    }
    return;
}

void ZPK::setZeros(const std::vector<std::complex<double>> &zeros)
{
    pImpl_->z = zeros;
    return;
}

std::vector<std::complex<double>> ZPK::getPoles(void) const
{
    return pImpl_->p;
}

std::vector<std::complex<double>> ZPK::getZeros(void) const
{
    return pImpl_->z;
}

void ZPK::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    pImpl_->tol = tol;
    return;
}
