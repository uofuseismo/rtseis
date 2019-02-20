#ifndef RTSEIS_UTILS_ZPK_HPP
#define RTSEIS_UTILS_ZPK_HPP 1
#include <complex>
#include <vector>
#include "rtseis/config.h"

namespace RTSeis
{
namespace Utilities
{
namespace FilterRepresentations
{
class ZPK
{
    public:
        ZPK(void);
        ZPK(const std::vector<std::complex<double>> &zeros,
            const std::vector<std::complex<double>> &poles,
            const double k);
        ZPK &operator=(const ZPK &zpk);
        bool operator==(const ZPK &zpk) const;
        bool operator!=(const ZPK &zpk) const;
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
        void setPoles(const std::vector<std::complex<double>> &poles);
        void setZeros(const size_t n, std::complex<double> zeros[]);
        void setZeros(const std::vector<std::complex<double>> &zeros);
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
}; // End ZPK
}; // End Representations
}; // End Utilities
}; // End RTSeis

#endif
