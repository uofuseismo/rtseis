#include <cmath>
#include <algorithm>
#include "rtseis/utilities/filterRepresentations/sos.hpp"
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;


SOS::SOS(void)
{
    clear();
    return;
}

SOS::SOS(const int ns,
         const std::vector<double> &bs,
         const std::vector<double> &as)
{
    clear();
    setSecondOrderSections(ns, bs, as);
    return;
}

SOS& SOS::operator=(const SOS &sos)
{
    if (&sos == this){return *this;}
    clear();
    bs_ = sos.bs_;
    as_ = sos.as_;
    ns_ = sos.ns_;
    tol_ = sos.tol_;
    return *this;
}

bool SOS::operator==(const SOS &sos) const
{
    if (bs_.size() != sos.bs_.size()){return false;}
    if (as_.size() != sos.as_.size()){return false;}
    if (ns_ != sos.ns_){return false;}
    for (size_t i=0; i<bs_.size(); i++)
    {
        if (std::abs(bs_[i] - sos.bs_[i]) > tol_){return false;}
    }
    for (size_t i=0; i<as_.size(); i++)
    {
        if (std::abs(as_[i] - sos.as_[i]) > tol_){return false;}
    }
    return true;
}

bool SOS::operator!=(const SOS &sos) const
{
    return !(*this == sos);
}

SOS::SOS(const SOS &sos)
{
   *this = sos;
   return;
}

SOS::~SOS(void)
{
    clear();
    return;
}

void SOS::clear(void)
{
    bs_.clear();
    as_.clear();
    ns_ = 0;
    tol_ = defaultTol_;
    return;
}

void SOS::print(FILE *fout)
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Numerator sections\n");
    for (int i=0; i<ns_; i++)
    {
        fprintf(f, "%+.16lf, %+.16lf, %+.16lf\n",
                bs_[3*i], bs_[3*i+1], bs_[3*i+2]);
    }
    fprintf(f, "Denominator sections\n");
    for (int i=0; i<ns_; i++)
    {
        fprintf(f, "%+.16lf, %+.16lf, %+.16lf\n",
                as_[3*i], as_[3*i+1], as_[3*i+2]);
    }
    return;
}

int SOS::setSecondOrderSections(const int ns,
                                const std::vector<double> &bs,
                                const std::vector<double> &as)
{
    clear();
    if (ns < 1)
    {
        RTSEIS_ERRMSG("%s", "No sections in SOS filter");
        return -1;
    }
    size_t ns3 = static_cast<size_t> (ns)*3;
    if (ns3 != bs.size())
    {
        RTSEIS_ERRMSG("bs.size() = %ld must equal 3*ns=%ld", bs.size(), ns3);
        return -1;
    }
    if (ns3 != as.size()) 
    {
        RTSEIS_ERRMSG("as.size() = %ld must equal 3*ns=%ld", as.size(), ns3);
        return -1;
    }
    for (int i=0; i<ns; i++)
    {
        if (bs[3*i] == 0)
        {
            RTSEIS_ERRMSG("Leading bs coefficient of section %d is zero", i);
            return -1;
        }
    }
    for (int i=0; i<ns; i++)
    {
        if (as[3*i] == 0)
        {
            RTSEIS_ERRMSG("Leading bs coefficient of section %d is zero", i);
            return -1;
        }
    }
    // It all checks out
    ns_ = ns;
    bs_ = bs;
    as_ = as;
    return 0;
}

std::vector<double> SOS::getNumeratorCoefficients(void) const
{
    return bs_;
}

std::vector<double> SOS::getDenominatorCoefficients(void) const
{
    return as_;
}

int SOS::getNumberOfSections(void) const
{
    return ns_;
}

void SOS::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    tol_ = tol;
    return;
}
