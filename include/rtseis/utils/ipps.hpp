#ifndef RTSEIS_UTILS_IPPS_HPP
#define RTSEIS_UTILS_IPPS_HPP 1
#include "rtseis/config.h"
#include "rtseis/enums.h"
#include <complex>
#include <vector>
#include <memory>

namespace RTSeis
{

    class IPPSFIRFilter
    {
        enum class Implementation
        {
            DIRECT = 0, /*!< Direct-form implementation. */
            FFT = 1     /*!< FFT overlap and add implementation. */ 
        };
        public:
            IPPSFIRFilter(void);
            ~IPPSFIRFilter(void);
            int initialize(const int nb, const double b[],
                           const bool lisRealTime = false,
                           const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE,
                           Implementation implementation = Implementation::DIRECT);
            int apply(const int n, const double x[], double y[]);
            int apply(const int n, const float x[], float y[]);
            void clear(void);
            bool isRealTime(void) const{return lrt_;}
            bool isInitialized(void) const{return linit_;}
        private:
            /*!< The filter state for double precision. */
            void *pFIRspec64_ = nullptr;
            /*!< The double precision taps. */
            void *pTaps64_ = nullptr;
            /*!< The input delay line. */
            void *dlysrc64_ = nullptr;
            /*!< The output delay line. */
            void *dlydst64_ = nullptr;
            /*!< The filter state for single precision. */
            void *pFIRspec32_ = nullptr;
            /*!< The double precision taps. */
            void *pTaps32_ = nullptr;
            /*!< The intpu delay line. */
            void *dlysrc32_ = nullptr;
            /*!< The output delay line. */
            void *dlydst32_ = nullptr;
            /*!< Workspace. */
            void *pBuf_ = nullptr;
            /*!< A copy of the input taps. */
            double *tapsRef_ = nullptr;
            /*!< A copy of the initial conditions.  This has dimension
                 [order]. */
            double *zi_ = nullptr;
            /*!< The number of taps. */
            int tapsLen_ = 0; 
            /*!< The length of the delay line which is max(128, order+1). */
            int nbDly_ = 0;
            /*!< Size of pBuf. */
            int bufferSize_ = 0;
            /*!< Size of state. */
            int specSize_ = 0;
            /*!< Filter order. */
            int order_ = 0;
            /*!< Implementation. */
            Implementation implementation_ = Implementation::DIRECT; 
            /*!< Precision of module. */
            enum rtseisPrecision_enum precision_ = RTSEIS_DOUBLE;
            /*!< Flag indicating the module is for real-time. */
            bool lrt_ = false;
            /*!< Flag indicating the module is initialized. */
            bool linit_ = false;
    };

namespace Utils
{
namespace IPPS
{
    /* Divides res = num/den. */
    int Div(std::vector<std::complex<double>> den,
            std::vector<std::complex<double>> num,
            std::vector<std::complex<double>> &res);
    /* Multiplies x = val*x. */
    int MulC(const std::complex<double> val,
             std::vector<std::complex<double>> x);

};
};
};

#endif
