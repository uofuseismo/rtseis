#ifndef RTSEIS_UTILS_BA_HPP
#define RTSEIS_UTILS_BA_HPP 1
#include "rtseis/config.h"
#include <stdio.h>
#include <complex>
#include <vector>

#ifdef __cplusplus
class BA
{
    public:
        BA(void);
        BA(const std::vector<double> b, const std::vector<double> a);
        BA &operator=(const BA &ba)
        {
            b_ = ba.b_;
            a_ = ba.a_;
            return *this;
        }
        bool operator==(const BA &ba) const
        {
            double tol = 1.e-12;
            if (b_.size() != ba.b_.size()){return false;}
            if (a_.size() != ba.a_.size()){return false;}
            for (size_t i=0; i<b_.size(); i++)
            {
                if (std::abs(b_[i] - ba.b_[i]) > tol){return false;}
            }
            for (size_t i=0; i<a_.size(); i++)
            {
                if (std::abs(a_[i] - ba.a_[i]) > tol){return false;}
            }
            return true;
        }
        ~BA(void);
        void print(FILE *fout = stdout);
        void clear(void);
        int getNumberOfNumeratorCoefficients(void) const;
        int getNumberOfDenominatorCoefficients(void) const;
        void setNumeratorCoefficients(const size_t n, double b[]);
        void setNumeratorCoefficients(const std::vector<double> b);
        void setDenominatorCoefficients(const size_t n, double a[]);
        void setDenominatorCoefficients(const std::vector<double> a);
        std::vector<double> getNumeratorCoefficients(void) const;
        std::vector<double> getDenominatorCoefficients(void) const;
    private:
        /*!< The numerator coefficients. */
        std::vector<double> b_;
        /*!< The denoninator coefficients. */
        std::vector<double> a_;
};
#endif

#endif
