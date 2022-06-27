#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "rtseis/filterRepresentations/ba.hpp"

using namespace RTSeis::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

class BA::BAImpl
{
public:
    /// Numerator coefficients
    std::vector<double> b; 
    /// Denominator coefficients 
    std::vector<double> a; 
    /// Default tolerance
    double tol = DEFAULT_TOL;
};

/// C'tor
BA::BA() :
    pImpl(std::make_unique<BAImpl>())
{
}

/// Copy c'tor
BA::BA(const BA &ba)
{
    *this = ba;
}

/// Move c'tor
BA::BA(BA &&ba) noexcept
{
    *this = std::move(ba);
}

/// Construct from numerator/denominator coefficients
BA::BA(const std::vector<double> &b, const std::vector<double> &a) :
    pImpl(std::make_unique<BAImpl>())
{
    setNumeratorCoefficients(b);
    setDenominatorCoefficients(a);
}

/// Copy assignment
BA& BA::operator=(const BA &ba)
{
    if (&ba == this){return *this;}
    pImpl = std::make_unique<BAImpl> (*ba.pImpl);
    return *this;
}

BA& BA::operator=(BA &&ba) noexcept
{
    if (&ba == this){return *this;}
    pImpl = std::move(ba.pImpl);
    return *this;
}

/// Destructor
BA::~BA() = default;

/// Equality
bool BA::operator==(const BA &ba) const
{
    if (pImpl->b.size() != ba.pImpl->b.size()){return false;}
    if (pImpl->a.size() != ba.pImpl->a.size()){return false;}
    for (size_t i = 0; i < pImpl->b.size(); i++)
    {
        if (std::abs(pImpl->b[i] - pImpl->b[i]) > pImpl->tol)
        {
            return false;
        }
    }
    for (size_t i = 0; i < pImpl->a.size(); i++)
    {
        if (std::abs(pImpl->a[i] - pImpl->a[i]) > pImpl->tol)
        {
            return false;
        }
    }
    return true;
}

/// Inequality
bool BA::operator!=(const BA &ba) const
{
    return !(*this == ba);
}

void BA::print(FILE *fout) const noexcept
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Numerator Coefficients:\n");
    for (const double b : pImpl->b)
    {
        fprintf(f, "%+.16lf\n", b);
    }
    if (!pImpl->a.empty())
    {
        fprintf(f, "Denominator Coefficients:\n");
        for (const double a : pImpl->a)
        {
            fprintf(f, "%+.16lf\n", a);
        }
    }
}

/// Reset class
void BA::clear() noexcept
{
    if (pImpl)
    {
        pImpl->b.clear();
        pImpl->a.clear();
        pImpl->tol = DEFAULT_TOL;
    }
}

[[maybe_unused]]
int BA::getNumberOfNumeratorCoefficients() const noexcept
{
    return static_cast<int> (pImpl->b.size());
}

[[maybe_unused]]
int BA::getNumberOfDenominatorCoefficients() const noexcept
{
    return static_cast<int> (pImpl->a.size());
}

void BA::setNumeratorCoefficients(const size_t n, const double b[])
{
    if (n < 1)
    {
        pImpl->b.resize(0);
        throw std::invalid_argument("b is empty");
    }
    if (n > 0 && b == nullptr)
    {
        pImpl->b.resize(0);
        throw std::invalid_argument("b is NULL");
    }
    pImpl->b.resize(n);
    std::copy(b, b+n, pImpl->b.begin());
}

void BA::setNumeratorCoefficients(const std::vector<double> &b)
{
    setNumeratorCoefficients(b.size(), b.data());
}

void BA::setDenominatorCoefficients(const size_t n, const double a[])
{
    if (n < 1)
    {
        pImpl->b.resize(0);
        throw std::invalid_argument("a is empty");
    }
    if (n > 0 && a == nullptr)
    {
        pImpl->a.resize(0);
        throw std::invalid_argument("a is NULL or empty");
    }
    if (n > 0 && a[0] == 0)
    {
        throw std::invalid_argument("a[0] = 0");
    }
    pImpl->a.resize(n);
    std::copy(a, a+n, pImpl->a.begin());
}

void BA::setDenominatorCoefficients(const std::vector<double> &a)
{
    setDenominatorCoefficients(a.size(), a.data());
}

std::vector<double> BA::getNumeratorCoefficients() const noexcept
{
    return pImpl->b;
}

std::vector<double> BA::getDenominatorCoefficients() const noexcept
{
    return pImpl->a;
}

void BA::setEqualityTolerance(const double tol)
{
    if (tol < 0)
    {
        std::cerr << "Tolerance is negative" << std::endl;
    }
    pImpl->tol = tol;
}

/*
bool BA::isZeroDenominator(void) const
{
    for (size_t i=0; i<a_.size(); i++)
    {
        if (a_[i] != 0){return false;}
    }
    return true;
}

bool BA::isFIR(void) const
{
    return isFIR_;
}
*/

std::ostream& RTSeis::FilterRepresentations::operator<<(
    std::ostream &os, const BA &ba)
{
    std::stringstream result;
    auto b = ba.getNumeratorCoefficients();
    result << "Numerator Coefficients:" << std::endl;
    for (const auto &bi : b)
    {
        result << std::setprecision(16) << bi << std::endl;
    }   
    auto a = ba.getDenominatorCoefficients();
    if (!a.empty())
    {   
        result << "Denominator Coefficients:" << std::endl;
        for (const auto &ai : a)
        {
            result << std::setprecision(16) << ai << std::endl;
        }
    }
    return os << result.str();
}
