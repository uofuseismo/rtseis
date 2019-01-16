#ifndef RTSEIS_UTILS_BA_HPP
#define RTSEIS_UTILS_BA_HPP 1
#include "rtseis/config.h"
#include <complex>
#include <vector>

class BA
{
    public:
        BA(void);
        BA(const std::vector<double> firTaps);
        BA(const std::vector<double> b, const std::vector<double> a);
        BA(const BA &ba);
        BA &operator=(const BA &ba);
        bool operator==(const BA &ba) const
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
            return true;
        }
        bool operator!=(const BA &ba) const
        {
            return !(*this == ba);
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
        void setEqualityTolerance(const double tol = 1.e-12);
        bool isZeroDenominator(void) const;
        bool isFIR(void) const;
    private:
        /*!< Default tolerance. */
        const double defaultTol_ = 1.e-12;
        /*!< The numerator coefficients. */
        std::vector<double> b_;
        /*!< The denoninator coefficients. */
        std::vector<double> a_;
        /*!< Tolerance in checking equality. */
        double tol_ = defaultTol_;
        /*!< Determines if the filter is an FIR filter. */
        bool isFIR_ = false;
};

#endif
