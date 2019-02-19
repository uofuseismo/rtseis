#ifndef RTSEIS_UTILS_SOS_HPP
#define RTSEIS_UTILS_SOS_HPP 1
#include <cstdio>
#include <vector>
#include "rtseis/config.h"

namespace RTSeis
{

class SOS
{
    public:
        SOS(void);
        SOS(const SOS &sos);
        SOS(const int ns,
            const std::vector<double> &bs,
            const std::vector<double> &as);
        SOS &operator=(const SOS &sos);
        bool operator==(const SOS &sos) const;
        bool operator!=(const SOS &sos) const;
        ~SOS(void);
        void print(FILE *fout = stdout);
        void clear(void);
        int setSecondOrderSections(const int ns,
                                   const std::vector<double> &bs,
                                   const std::vector<double> &as);
        std::vector<double> getNumeratorCoefficients(void) const;
        std::vector<double> getDenominatorCoefficients(void) const;
        int getNumberOfSections(void) const;
        void setEqualityTolerance(const double tol = 1.e-12);
    private:
        /*!< The numerator sections. */
        std::vector<double> bs_; 
        /*!< The denominator sections. */
        std::vector<double> as_; 
        /*!< The number of sections. */
        int ns_ = 0; 
        /*!< Default tolerance. */
        const double defaultTol_ = 1.e-12;
        /*!< Tolerance in checking equality. */
        double tol_ = defaultTol_;
};

};

#endif
