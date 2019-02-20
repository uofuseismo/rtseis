#define RTSEIS_LOGGING 1
#include <cmath>
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::FilterRepresentations;

BA::BA(void)
{
    clear();
    return;
}

BA::BA(const std::vector<double> &b, const std::vector<double> &a)
{
    clear();
    setNumeratorCoefficients(b);
    setDenominatorCoefficients(a);
    return;
}

BA::BA(const std::vector<double> &firTaps)
{
    clear();
    setNumeratorCoefficients(firTaps);
    std::vector<double> a({1});
    setDenominatorCoefficients(a); 
    return;
}

BA& BA:: operator=(const BA &ba)
{
    if (&ba == this){return *this;}
    clear();
    b_ = ba.b_;
    a_ = ba.a_;
    tol_ = ba.tol_;
    isFIR_ = ba.isFIR_;
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
    if (b_.size() != ba.b_.size()){return false;}
    if (a_.size() != ba.a_.size()){return false;}
    for (size_t i=0; i<b_.size(); i++)
    {
        if (std::abs(b_[i] - ba.b_[i]) > tol_){return false;}
    }
    for (size_t i=0; i<a_.size(); i++)
    {
        if (std::abs(a_[i] - ba.a_[i]) > tol_){return false;}
    }
    if (isFIR_ != ba.isFIR_){return false;}
    return true;
}

bool BA::operator!=(const BA &ba) const
{
    return !(*this == ba);
}

void BA::print(FILE *fout)
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    if (!isFIR())
    {
        fprintf(f, "Numerator Coefficients:\n");
    }
    else
    {
        fprintf(f, "FIR Coefficients:\n");
    }
    for (size_t i=0; i<b_.size(); i++)
    {
        fprintf(f, "%+.16lf\n", b_[i]);
    }
    if (!isFIR())
    {
        fprintf(f, "Denominator Coefficients:\n");
        for (size_t i=0; i<a_.size(); i++)
        {
            fprintf(f, "%+.16lf\n", a_[i]);
        }
    }
    return;
}

void BA::clear(void)
{
    b_.clear();
    a_.clear();
    tol_ = defaultTol_;
    isFIR_ = false;
    return;
}

int BA::getNumberOfNumeratorCoefficients(void) const
{
    return static_cast<int> (b_.size());
}

int BA::getNumberOfDenominatorCoefficients(void) const
{
    return static_cast<int> (a_.size());
}

void BA::setNumeratorCoefficients(const size_t n, double b[])
{
    if (n > 0 && b == nullptr)
    {
        RTSEIS_ERRMSG("%s", "b is null");
        b_.resize(0);
        return;
    }
    b_.resize(n);
    for (size_t i=0; i<n; i++)
    {
        b_[i] = b[i];
    }
    return;
}

void BA::setNumeratorCoefficients(const std::vector<double> &b)
{
    b_ = b;
    return;
}

void BA::setDenominatorCoefficients(const size_t n, double a[])
{
    isFIR_ = false;
    if (n > 0 && a == nullptr)
    {
        RTSEIS_ERRMSG("%s", "a is null");
        a_.resize(0);
        return;
    }
    if (n > 0 && a[0] == 0){RTSEIS_WARNMSG("%s", "a[0] = 0");}
    a_.resize(n);
    for (size_t i=0; i<n; i++)
    {
        a_[i] = a[i];
    }
    if (a_.size() == 1)
    {
        if (std::abs(a_[0] - 1) < 1.e-14){isFIR_ = true;}
    }
    return;
}

void BA::setDenominatorCoefficients(const std::vector<double> &a)
{
    isFIR_ = false;
    if (a.size() > 0)
    {
        if (a[0] == 0){RTSEIS_WARNMSG("%s", "a[0] = 0");}
    }
    a_ = a;
    if (a_.size() == 1)
    {
        if (std::abs(a_[0] - 1) < 1.e-14){isFIR_ = true;}
    }
    return;
}

std::vector<double> BA::getNumeratorCoefficients(void) const
{
    return b_;
}

std::vector<double> BA::getDenominatorCoefficients(void) const
{
    return a_;
}

void BA::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    tol_ = tol;
    return;
}

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
