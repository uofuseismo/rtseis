#define RTSEIS_LOGGING 1
#include <cmath>
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

struct FIR::FIRImpl
{
    /// Filter taps
    std::vector<double> pTaps;
    /// Default tolerance
    double tol = DEFAULT_TOL;
};

FIR::FIR(void) :
    pImpl_(new FIRImpl())
{
    return;
}

FIR::FIR(const std::vector<double> &pTaps) :
    pImpl_(new FIRImpl())
{
    setFilterTaps(pTaps);
    return;
}

FIR& FIR::operator=(const FIR &fir)
{
    if (&fir == this){return *this;}
    pImpl_ = std::unique_ptr<FIRImpl> (new FIRImpl());
    pImpl_->pTaps = fir.pImpl_->pTaps;
    pImpl_->tol   = fir.pImpl_->tol;
    return *this;
}

FIR::FIR(const FIR &fir)
{
    *this = fir;
    return;
}

FIR::~FIR(void)
{
    clear();
    return;
}

bool FIR::operator==(const FIR &fir) const
{
    if (pImpl_->pTaps.size() != fir.pImpl_->pTaps.size()){return false;}
    for (size_t i=0; i<pImpl_->pTaps.size(); i++)
    {   
        if (std::abs(pImpl_->pTaps[i] - pImpl_->pTaps[i]) > pImpl_->tol)
        {
            return false;
        }
    }   
    return true;
}

bool FIR::operator!=(const FIR &fir) const
{
    return !(*this == fir);
}

void FIR::print(FILE *fout)
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Filter Coefficients:\n");
    for (size_t i=0; i<pImpl_->pTaps.size(); i++)
    {
        fprintf(f, "%+.16lf\n", pImpl_->pTaps[i]);
    }
    return;
}

void FIR::clear(void)
{
    pImpl_->pTaps.clear();
    pImpl_->tol = DEFAULT_TOL;
    return;
}

void FIR::setFilterTaps(const size_t n, double pTaps[])
{
    if (n > 0 && pTaps == nullptr)
    {
        RTSEIS_ERRMSG("%s", "pTaps is null");
        pImpl_->pTaps.resize(0);
        return;
    }
    pImpl_->pTaps.resize(n);
    std::copy(pTaps, pTaps+n, pImpl_->pTaps.begin());
    return;
}

void FIR::setFilterTaps(const std::vector<double> &pTaps)
{
    pImpl_->pTaps = pTaps;
    return;
}

std::vector<double> FIR::getFilterTaps(void) const
{
    return pImpl_->pTaps;
}

void FIR::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    pImpl_->tol = tol;
    return;
}
