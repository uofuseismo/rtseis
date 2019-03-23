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

void FIR::print(FILE *fout) const noexcept
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

void FIR::clear(void) noexcept
{
    pImpl_->pTaps.clear();
    pImpl_->tol = DEFAULT_TOL;
    return;
}

int FIR::getNumberOfFilterTaps(void) const noexcept
{
    int ntaps = static_cast<int> (pImpl_->pTaps.size());
    return ntaps;
}

void FIR::setFilterTaps(const size_t n, const double pTaps[])
{
    if (n < 1)
    {
        pImpl_->pTaps.resize(0);
        throw std::invalid_argument("No filter taps");
    }
    if (n > 0 && pTaps == nullptr)
    {
        pImpl_->pTaps.resize(0);
        RTSEIS_ERRMSG("%s", "pTaps is null");
        throw std::invalid_argument("pTaps is NULL");
        return;
    }
    pImpl_->pTaps.resize(n);
    std::copy(pTaps, pTaps+n, pImpl_->pTaps.begin());
    return;
}

void FIR::setFilterTaps(const std::vector<double> &pTaps)
{
    setFilterTaps(pTaps.size(), pTaps.data());
    //pImpl_->pTaps = pTaps;
    return;
}

std::vector<double> FIR::getFilterTaps(void) const noexcept
{
    return pImpl_->pTaps;
}

void FIR::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    pImpl_->tol = tol;
    return;
}
