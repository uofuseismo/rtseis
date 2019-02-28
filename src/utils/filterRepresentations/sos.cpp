#include <cmath>
#include <algorithm>
#include "rtseis/utilities/filterRepresentations/sos.hpp"
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

struct SOS::SOSImpl
{
    /// Numerator coefficients
    std::vector<double> bs;
    /// Denominator coefficients
    std::vector<double> as;
    /// Number of sections
    int ns = 0;
    /// Equality tolerance
    double tol = DEFAULT_TOL;
};

SOS::SOS(void) :
    pImpl_(new SOSImpl())
{
    return;
}

SOS::SOS(const int ns,
         const std::vector<double> &bs,
         const std::vector<double> &as) :
    pImpl_(new SOSImpl())
{
    setSecondOrderSections(ns, bs, as);
    return;
}

SOS& SOS::operator=(const SOS &sos)
{
    if (&sos == this){return *this;}
    pImpl_ = std::unique_ptr<SOSImpl> (new SOSImpl());
    pImpl_->bs  = sos.pImpl_->bs;
    pImpl_->as  = sos.pImpl_->as;
    pImpl_->ns  = sos.pImpl_->ns;
    pImpl_->tol = sos.pImpl_->tol;
    return *this;
}

bool SOS::operator==(const SOS &sos) const
{
    if (pImpl_->bs.size() != sos.pImpl_->bs.size()){return false;}
    if (pImpl_->as.size() != sos.pImpl_->as.size()){return false;}
    if (pImpl_->ns != sos.pImpl_->ns){return false;}
    for (size_t i=0; i<pImpl_->bs.size(); i++)
    {
        if (std::abs(pImpl_->bs[i] - sos.pImpl_->bs[i]) > pImpl_->tol)
        {
            return false;
        }
    }
    for (size_t i=0; i<pImpl_->as.size(); i++)
    {
        if (std::abs(pImpl_->as[i] - sos.pImpl_->as[i]) > pImpl_->tol)
        {
            return false;
        }
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
    pImpl_->bs.clear();
    pImpl_->as.clear();
    pImpl_->ns = 0;
    pImpl_->tol = DEFAULT_TOL;
    return;
}

void SOS::print(FILE *fout)
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Numerator sections\n");
    for (int i=0; i<pImpl_->ns; i++)
    {
        fprintf(f, "%+.16lf, %+.16lf, %+.16lf\n",
                pImpl_->bs[3*i], pImpl_->bs[3*i+1], pImpl_->bs[3*i+2]);
    }
    fprintf(f, "Denominator sections\n");
    for (int i=0; i<pImpl_->ns; i++)
    {
        fprintf(f, "%+.16lf, %+.16lf, %+.16lf\n",
                pImpl_->as[3*i], pImpl_->as[3*i+1], pImpl_->as[3*i+2]);
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
    pImpl_->ns = ns;
    pImpl_->bs = bs;
    pImpl_->as = as;
    return 0;
}

std::vector<double> SOS::getNumeratorCoefficients(void) const
{
    return pImpl_->bs;
}

std::vector<double> SOS::getDenominatorCoefficients(void) const
{
    return pImpl_->as;
}

int SOS::getNumberOfSections(void) const
{
    return pImpl_->ns;
}

void SOS::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    pImpl_->tol = tol;
    return;
}
