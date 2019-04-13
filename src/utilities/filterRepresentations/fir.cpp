#define RTSEIS_LOGGING 1
#include <cmath>
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

class FIR::FIRImpl
{
public:
    /// Filter taps
    std::vector<double> pTaps;
    /// Default tolerance
    double tol = DEFAULT_TOL;
};

FIR::FIR(void) :
    pImpl(std::make_unique<FIRImpl>())
{
    return;
}

FIR::FIR(const std::vector<double> &pTaps) :
    pImpl(std::make_unique<FIRImpl>())
{
    setFilterTaps(pTaps);
    return;
}

FIR& FIR::operator=(const FIR &fir)
{
    if (&fir == this){return *this;}
    pImpl = std::make_unique<FIRImpl> ();
    //pImpl = std::unique_ptr<FIRImpl> (new FIRImpl());
    pImpl->pTaps = fir.pImpl->pTaps;
    pImpl->tol   = fir.pImpl->tol;
    return *this;
}

FIR& FIR::operator=(FIR &&fir)
{
    if (&fir == this){return *this;}
    pImpl = std::move(fir.pImpl);
    return *this;
}

FIR::FIR(const FIR &fir)
{
    *this = fir;
    return;
}

FIR::~FIR(void) = default;

bool FIR::operator==(const FIR &fir) const
{
    if (pImpl->pTaps.size() != fir.pImpl->pTaps.size()){return false;}
    for (size_t i=0; i<pImpl->pTaps.size(); i++)
    {   
        if (std::abs(pImpl->pTaps[i] - pImpl->pTaps[i]) > pImpl->tol)
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
    for (size_t i=0; i<pImpl->pTaps.size(); i++)
    {
        fprintf(f, "%+.16lf\n", pImpl->pTaps[i]);
    }
    return;
}

void FIR::clear(void) noexcept
{
    pImpl->pTaps.clear();
    pImpl->tol = DEFAULT_TOL;
    return;
}

int FIR::getNumberOfFilterTaps(void) const noexcept
{
    int ntaps = static_cast<int> (pImpl->pTaps.size());
    return ntaps;
}

void FIR::setFilterTaps(const size_t n, const double pTaps[])
{
    if (n < 1)
    {
        pImpl->pTaps.resize(0);
        throw std::invalid_argument("No filter taps");
    }
    if (n > 0 && pTaps == nullptr)
    {
        pImpl->pTaps.resize(0);
        RTSEIS_ERRMSG("%s", "pTaps is null");
        throw std::invalid_argument("pTaps is NULL");
        return;
    }
    pImpl->pTaps.resize(n);
    std::copy(pTaps, pTaps+n, pImpl->pTaps.begin());
    return;
}

void FIR::setFilterTaps(const std::vector<double> &pTaps)
{
    setFilterTaps(pTaps.size(), pTaps.data());
    //pImpl->pTaps = pTaps;
    return;
}

std::vector<double> FIR::getFilterTaps(void) const noexcept
{
    return pImpl->pTaps;
}

void FIR::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    pImpl->tol = tol;
    return;
}
