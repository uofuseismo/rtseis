#define RTSEIS_LOGGING 1
#include <cmath>
#include <algorithm>
#include "rtseis/utilities/filterRepresentations/zpk.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;

ZPK::ZPK(void)
{
    clear();
    return;
}

ZPK::ZPK(const std::vector<std::complex<double>> &zeros,
         const std::vector<std::complex<double>> &poles,
         const double k)
{
    clear();
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
    clear();
    z_ = zpk.z_;
    p_ = zpk.p_;
    k_ = zpk.k_;
    tol_ = zpk.tol_;
    return *this;
}

bool ZPK::operator==(const ZPK &zpk) const
{
    if (p_.size() != zpk.p_.size()){return false;}
    if (z_.size() != zpk.z_.size()){return false;}
    for (size_t i=0; i<p_.size(); i++)
    {
        if (std::abs(p_[i] - zpk.p_[i]) > tol_){return false;}
    }
    for (size_t i=0; i<z_.size(); i++)
    {
        if (std::abs(z_[i] - zpk.z_[i]) > tol_){return false;}
    }
    if (std::abs(k_ - zpk.k_) > tol_){return false;}
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
    fprintf(f, "Gain: %.16lf\n", k_);
    fprintf(f, "Zeros:\n");
    for (size_t i=0; i<z_.size(); i++)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n", std::real(z_[i]), std::imag(z_[i])); 
    }
    fprintf(f, "Poles:\n");
    for (size_t i=0; i<p_.size(); i++)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n", std::real(p_[i]), std::imag(p_[i]));
    }
    return;
}

void ZPK::sortPoles(const bool ascending)
{
    if (ascending)
    {
        std::sort(p_.begin(), p_.end(),
                  // Begin lambda
                  [](std::complex<double> a,
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
    }
    else
    {
        std::sort(p_.begin(), p_.end(),
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
        std::sort(z_.begin(), z_.end(),
                  // Begin lambda
                  [](std::complex<double> a,
                     std::complex<double> b)
                  {
                     return (std::abs(a) < std::abs(b));
                  });
    }   
    else
    {
        std::sort(z_.begin(), z_.end(),
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
    p_.clear();
    z_.clear();
    k_ = 0;
    tol_ = defaultTol_;
    return;
}

void ZPK::setGain(const double k)
{
    k_ = k;
    return;
}

double ZPK::getGain(void) const
{
    return k_;
}

int ZPK::getNumberOfPoles(void) const
{
    return static_cast<int> (p_.size());
}

int ZPK::getNumberOfZeros(void) const
{
    return static_cast<int> (z_.size());
}

void ZPK::setPoles(const size_t n, std::complex<double> poles[])
{
    if (n > 0 && poles == nullptr)
    {
        RTSEIS_ERRMSG("%s", "Poles is null");
        p_.resize(0);
        return;
    }
    p_.resize(n);
    for (size_t i=0; i<n; i++)
    {
        p_[i] = poles[i];
    }
    return;
}

void ZPK::setPoles(const std::vector<std::complex<double>> &poles)
{
    p_ = poles;
    return;
}

void ZPK::setZeros(const size_t n, std::complex<double> zeros[])
{
    if (n > 0 && zeros == nullptr)
    {
        RTSEIS_ERRMSG("%s", "Zeros is null");
        z_.resize(0);
        return;
    }
    z_.resize(n);
    for (size_t i=0; i<n; i++)
    {
        z_[i] = zeros[i];
    }
    return;
}

void ZPK::setZeros(const std::vector<std::complex<double>> &zeros)
{
    z_ = zeros;
    return;
}

std::vector<std::complex<double>> ZPK::getPoles(void) const
{
    return p_;
}

std::vector<std::complex<double>> ZPK::getZeros(void) const
{
    return z_;
}

void ZPK::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    tol_ = tol;
    return;
}
