#ifndef RTSEIS_UTILS_ZPK_HPP
#define RTSEIS_UTILS_ZPK_HPP 1
#include "rtseis/config.h"
#include <complex>
#include <vector>

class ZPK
{
    public:
        ZPK(void);
        ZPK(const std::vector<std::complex<double>> zeros,
            const std::vector<std::complex<double>> poles,
            const double k);
        ZPK &operator=(const ZPK &zpk);
        bool operator==(const ZPK &zpk) const
        {
            if (p_.size() != zpk.p_.size()){return false;}
            if (z_.size() != zpk.z_.size()){return false;}
            for (size_t i=0; i<p_.size(); i++)
            {
                if (std::abs(p_[i] - zpk.p_[i]) > tol_){return false;}
            }
            for (size_t i=0; i<z_.size(); i++)
            {
                if (std::abs(z_[i] - zpk.z_[i]) > tol_){return false;}
            }
            if (std::abs(k_ - zpk.k_) > tol_){return false;}
            return true;
        }
        bool operator!=(const ZPK &zpk) const
        {
            return !(*this == zpk);
        }
        ~ZPK(void);
        void sortPoles(bool ascending=true);
        void sortZeros(bool ascending=true);
        void print(FILE *fout = stdout);
        void clear(void);
        void setGain(const double k);
        double getGain(void) const;
        int getNumberOfPoles(void) const;
        int getNumberOfZeros(void) const;
        void setPoles(const size_t n, std::complex<double> poles[]);
        void setPoles(const std::vector<std::complex<double>> poles);
        void setZeros(const size_t n, std::complex<double> zeros[]);
        void setZeros(const std::vector<std::complex<double>> zeros);
        std::vector<std::complex<double>> getPoles(void) const;
        std::vector<std::complex<double>> getZeros(void) const;
        void setEqualityTolerance(const double tol = 1.e-12);
    private:
        /*!< The zeros. */
        std::vector<std::complex<double>> z_;
        /*!< The poles. */
        std::vector<std::complex<double>> p_;
        /*!< The gain. */
        double k_ = 0;
        /*!< Default tolerance. */
        const double defaultTol_ = 1.e-12;
        /*!< Tolerance in checking equality. */
        double tol_ = defaultTol_;
};

#endif
