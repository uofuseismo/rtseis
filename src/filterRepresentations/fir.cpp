#include <iostream>
#include <cmath>
#include "rtseis/filterRepresentations/fir.hpp"

using namespace RTSeis::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

class FIR::FIRImpl
{
public:
    /// Filter taps
    std::vector<double> pTaps;
    /// Default tolerance
    double tol = DEFAULT_TOL;
};

/// C'tor
FIR::FIR() :
    pImpl(std::make_unique<FIRImpl>())
{
}

/// Construct from filter taps 
FIR::FIR(const std::vector<double> &pTaps) :
    pImpl(std::make_unique<FIRImpl>())
{
    setFilterTaps(pTaps);
}

/// Copy assignment
FIR& FIR::operator=(const FIR &fir)
{
    if (&fir == this){return *this;}
    pImpl = std::make_unique<FIRImpl> (*fir.pImpl);
    return *this;
}

/// Move assignmnet
FIR& FIR::operator=(FIR &&fir) noexcept
{
    if (&fir == this){return *this;}
    pImpl = std::move(fir.pImpl);
    return *this;
}

/// Copy c'tor
FIR::FIR(const FIR &fir)
{
    *this = fir;
}

/// Move c'tor
FIR::FIR(FIR &&fir) noexcept
{
    *this = std::move(fir);
}

/// Destructor
FIR::~FIR() = default;

/// Equality operator
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

/// Inequality operator
bool FIR::operator!=(const FIR &fir) const
{
    return !(*this == fir);
}

/// Write to file
void FIR::print(FILE *fout) const noexcept
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Filter Coefficients:\n");
    for (const double pTap : pImpl->pTaps)
    {
        fprintf(f, "%+.16lf\n", pTap);
    }
}

/// Clear
void FIR::clear() noexcept
{
    pImpl->pTaps.clear();
    pImpl->tol = DEFAULT_TOL;
}

int FIR::getNumberOfFilterTaps() const noexcept
{
    return static_cast<int> (pImpl->pTaps.size());;
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
        throw std::invalid_argument("pTaps is NULL");
    }
    pImpl->pTaps.resize(n);
    std::copy(pTaps, pTaps+n, pImpl->pTaps.begin());
}

void FIR::setFilterTaps(const std::vector<double> &pTaps)
{
    setFilterTaps(pTaps.size(), pTaps.data());
}

std::vector<double> FIR::getFilterTaps() const noexcept
{
    return pImpl->pTaps;
}

[[maybe_unused]]
void FIR::setEqualityTolerance(const double tol)
{
    if (tol < 0)
    {
        std::cerr << "Tolerance is negative" << std::endl;
    }
    pImpl->tol = tol;
}
