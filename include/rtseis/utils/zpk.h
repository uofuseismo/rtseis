#ifndef RTSEIS_UTILS_ZPK_HPP
#define RTSEIS_UTILS_ZPK_HPP 1
#include "rtseis/config.h"
#include <stdio.h>
#include <complex>
#include <vector>

#ifdef __cplusplus
class ZPK
{
    public:
        ZPK(void);
        ZPK(const std::vector<std::complex<double>> zeros,
            const std::vector<std::complex<double>> poles,
            const double k);
        ZPK &operator=(const ZPK &zpk)
        {
            z_ = zpk.z_;
            p_ = zpk.p_;
            k_ = zpk.k_;
            return *this;
        }
        bool operator==(const ZPK &zpk) const
        {
            double tol = 1.e-12;
            if (p_.size() != zpk.p_.size()){return false;}
            if (z_.size() != zpk.z_.size()){return false;}
            for (size_t i=0; i<p_.size(); i++)
            {
                if (std::abs(p_[i] - zpk.p_[i]) > tol){return false;}
            }
            for (size_t i=0; i<z_.size(); i++)
            {
                if (std::abs(z_[i] - zpk.z_[i]) > tol){return false;}
            }
            if (std::abs(k_ - zpk.k_) > tol){return false;}
            return true;
        }
        ~ZPK(void);
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
    private:
        /*!< The zeros. */
        std::vector<std::complex<double>> z_;
        /*!< The poles. */
        std::vector<std::complex<double>> p_;
        /*!< The gain. */
        double k_ = 0;
};
#endif

#endif
