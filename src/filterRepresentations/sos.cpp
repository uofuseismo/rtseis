#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "rtseis/filterRepresentations/sos.hpp"

using namespace RTSeis::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

class SOS::SOSImpl
{
public:
    /// Numerator coefficients
    std::vector<double> bs;
    /// Denominator coefficients
    std::vector<double> as;
    /// Number of sections
    int ns = 0;
    /// Equality tolerance
    double tol = DEFAULT_TOL;
};

/// Constructor
SOS::SOS() :
    pImpl(std::make_unique<SOSImpl>())
{
}

/// Constructor
SOS::SOS(const int ns,
         const std::vector<double> &bs,
         const std::vector<double> &as) :
    pImpl(std::make_unique<SOSImpl>())
{
    setSecondOrderSections(ns, bs, as);
}

/// Copy assignment
SOS& SOS::operator=(const SOS &sos)
{
    if (&sos == this){return *this;}
    pImpl = std::make_unique<SOSImpl> (*sos.pImpl);
    return *this;
}

/// Move assignment
SOS& SOS::operator=(SOS &&sos) noexcept
{
    if (&sos == this){return *this;}
    pImpl = std::move(sos.pImpl);
    return *this;
}

/// Equality operator
bool SOS::operator==(const SOS &sos) const
{
    if (pImpl->bs.size() != sos.pImpl->bs.size()){return false;}
    if (pImpl->as.size() != sos.pImpl->as.size()){return false;}
    if (pImpl->ns != sos.pImpl->ns){return false;}
    for (size_t i=0; i<pImpl->bs.size(); i++)
    {
        if (std::abs(pImpl->bs[i] - sos.pImpl->bs[i]) > pImpl->tol)
        {
            return false;
        }
    }
    for (size_t i=0; i<pImpl->as.size(); i++)
    {
        if (std::abs(pImpl->as[i] - sos.pImpl->as[i]) > pImpl->tol)
        {
            return false;
        }
    }
    return true;
}

/// Inequality operator
bool SOS::operator!=(const SOS &sos) const
{
    return !(*this == sos);
}

/// Copy c'tor
SOS::SOS(const SOS &sos)
{
   *this = sos;
}

/// Move c'tor
SOS::SOS(SOS &&sos) noexcept
{
    *this = std::move(sos);
}

/// Destructor
SOS::~SOS() = default;

/// Reset class
void SOS::clear() noexcept
{
    pImpl->bs.clear();
    pImpl->as.clear();
    pImpl->ns = 0;
    pImpl->tol = DEFAULT_TOL;
}

void SOS::print(FILE *fout) const noexcept
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Numerator sections\n");
    for (int i = 0; i < pImpl->ns; i++)
    {
        fprintf(f, "%+.16lf, %+.16lf, %+.16lf\n",
                pImpl->bs[3*i], pImpl->bs[3*i+1], pImpl->bs[3*i+2]);
    }
    fprintf(f, "Denominator sections\n");
    for (int i = 0; i < pImpl->ns; i++)
    {
        fprintf(f, "%+.16lf, %+.16lf, %+.16lf\n",
                pImpl->as[3*i], pImpl->as[3*i+1], pImpl->as[3*i+2]);
    }
}

void SOS::setSecondOrderSections(const int ns,
                                 const std::vector<double> &bs,
                                 const std::vector<double> &as)
{
    clear();
    if (ns < 1)
    {
        throw std::invalid_argument("No sections in SOS filter");
    }
    size_t ns3 = static_cast<size_t> (ns)*3;
    if (ns3 != bs.size())
    {
        throw std::invalid_argument("ba.size() = " + std::to_string(bs.size())
                                 + " must equal 3*ns = " + std::to_string(ns3));
    }
    if (ns3 != as.size()) 
    {
        throw std::invalid_argument("as.size() = " + std::to_string(as.size())
                                 + " must equal 3*ns = " + std::to_string(ns3));
    }
    for (int i=0; i<ns; i++)
    {
        if (bs[3*i] == 0)
        {
            throw std::invalid_argument("Leading bs coefficient of section "
                                     + std::to_string(i+1) + "is zero");
        }
    }
    for (int i=0; i<ns; i++)
    {
        if (as[3*i] == 0)
        {
            throw std::invalid_argument("Leading as coefficient of section "
                                     + std::to_string(i+1) + "is zero");
        }
    }
    // It all checks out
    pImpl->ns = ns;
    pImpl->bs = bs;
    pImpl->as = as;
}

std::vector<double> SOS::getNumeratorCoefficients() const noexcept
{
    return pImpl->bs;
}

std::vector<double> SOS::getDenominatorCoefficients() const noexcept
{
    return pImpl->as;
}

int SOS::getNumberOfSections() const noexcept
{
    return pImpl->ns;
}

void SOS::setEqualityTolerance(const double tol)
{
    if (tol < 0)
    {
        std::cerr << "Tolerance is negative" << std::endl;
    }
    pImpl->tol = tol;
}

std::ostream& RTSeis::FilterRepresentations::operator<<(
    std::ostream &os, const SOS &sos)
{
    std::stringstream result;
    result << "Numerator sections" << std::endl;
    auto nSections = sos.getNumberOfSections();
    auto bs = sos.getNumeratorCoefficients(); 
    for (int i = 0; i < nSections; ++i)
    {
        result << std::setprecision(16)
               << bs.at(3*i)   << ","
               << bs.at(3*i+1) << ","
               << bs.at(3*i+2) << std::endl;
    } 
    auto as = sos.getDenominatorCoefficients();
    result << "Denominator sections" << std::endl;
    for (int i = 0; i < nSections; ++i)
    {
        result << std::setprecision(16)
               << as.at(3*i)   << "," 
               << as.at(3*i+1) << "," 
               << as.at(3*i+2) << std::endl;
    }
    return os << result.str();
}
