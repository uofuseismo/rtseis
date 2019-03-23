#define RTSEIS_LOGGING 1
#include <cmath>
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;

#define DEFAULT_TOL 1.e-12

struct BA::BAImpl
{
    /// Numerator coefficients
    std::vector<double> b; 
    /// Denominator coefficients 
    std::vector<double> a; 
    /// Default tolerance
    double tol = DEFAULT_TOL;
};

BA::BA(void) :
    pImpl_(new BAImpl())
{
    return;
}

BA::BA(const std::vector<double> &b, const std::vector<double> &a) :
    pImpl_(new BAImpl())
{
    setNumeratorCoefficients(b);
    setDenominatorCoefficients(a);
    return;
}

BA& BA::operator=(const BA &ba)
{
    if (&ba == this){return *this;}
    pImpl_ = std::unique_ptr<BAImpl> (new BAImpl());
    pImpl_->b   = ba.pImpl_->b;
    pImpl_->a   = ba.pImpl_->a;
    pImpl_->tol = ba.pImpl_->tol;
    return *this;
}

BA::BA(const BA &ba)
{
    *this = ba;
    return;
}

BA::~BA(void)
{
    clear();
    return;
}

bool BA::operator==(const BA &ba) const
{
    if (pImpl_->b.size() != ba.pImpl_->b.size()){return false;}
    if (pImpl_->a.size() != ba.pImpl_->a.size()){return false;}
    for (size_t i=0; i<pImpl_->b.size(); i++)
    {
        if (std::abs(pImpl_->b[i] - pImpl_->b[i]) > pImpl_->tol)
        {
            return false;
        }
    }
    for (size_t i=0; i<pImpl_->a.size(); i++)
    {
        if (std::abs(pImpl_->a[i] - pImpl_->a[i]) > pImpl_->tol)
        {
            return false;
        }
    }
    return true;
}

bool BA::operator!=(const BA &ba) const
{
    return !(*this == ba);
}

void BA::print(FILE *fout) const noexcept
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Numerator Coefficients:\n");
    for (size_t i=0; i<pImpl_->b.size(); i++)
    {
        fprintf(f, "%+.16lf\n", pImpl_->b[i]);
    }
    if (pImpl_->a.size() > 0)
    {
        fprintf(f, "Denominator Coefficients:\n");
        for (size_t i=0; i<pImpl_->a.size(); i++)
        {
            fprintf(f, "%+.16lf\n", pImpl_->a[i]);
        }
    }
    return;
}

void BA::clear(void) noexcept
{
    if (pImpl_)
    {
        pImpl_->b.clear();
        pImpl_->a.clear();
        pImpl_->tol = DEFAULT_TOL;
    }
    return;
}

int BA::getNumberOfNumeratorCoefficients(void) const noexcept
{
    return static_cast<int> (pImpl_->b.size());
}

int BA::getNumberOfDenominatorCoefficients(void) const noexcept
{
    return static_cast<int> (pImpl_->a.size());
}

void BA::setNumeratorCoefficients(const size_t n, const double b[])
{
    if (n < 1)
    {
        pImpl_->b.resize(0);
        throw std::invalid_argument("b is empty");
    }
    if (n > 0 && b == nullptr)
    {
        pImpl_->b.resize(0);
        RTSEIS_ERRMSG("%s", "b is null");
        throw std::invalid_argument("b is NULL");
        return;
    }
    pImpl_->b.resize(n);
    std::copy(b, b+n, pImpl_->b.begin());
    return;
}

void BA::setNumeratorCoefficients(const std::vector<double> &b)
{
    setNumeratorCoefficients(b.size(), b.data());
    //pImpl_->b = b;
    return;
}

void BA::setDenominatorCoefficients(const size_t n, const double a[])
{
    if (n < 1)
    {
        pImpl_->b.resize(0);
        throw std::invalid_argument("a is empty");
    }
    if (n > 0 && a == nullptr)
    {
        pImpl_->a.resize(0);
        RTSEIS_ERRMSG("%s", "a is null");
        throw std::invalid_argument("a is NULL or empty");
        return;
    }
    if (n > 0 && a[0] == 0)
    {
        RTSEIS_ERRMSG("%s", "a[0] = 0");
        throw std::invalid_argument("a[0] = 0");
    }
    pImpl_->a.resize(n);
    std::copy(a, a+n, pImpl_->a.begin());
    return;
}

void BA::setDenominatorCoefficients(const std::vector<double> &a)
{
    setDenominatorCoefficients(a.size(), a.data());
/*
    if (pImpl_->a.size() > 0)
    {
        if (a[0] == 0)
        {
            RTSEIS_WARNMSG("%s", "a[0] = 0");
        }
    }
    pImpl_->a = a;
*/
    return;
}

std::vector<double> BA::getNumeratorCoefficients(void) const noexcept
{
    return pImpl_->b;
}

std::vector<double> BA::getDenominatorCoefficients(void) const noexcept
{
    return pImpl_->a;
}

void BA::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    pImpl_->tol = tol;
    return;
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
