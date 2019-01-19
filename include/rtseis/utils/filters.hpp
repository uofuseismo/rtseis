#ifndef RTSEIS_UTILS_FILTERS_HPP
#define RTSEIS_UTILS_FILTERS_HPP 1
#include <stdio.h>
#include "rtseis/config.h"
#include "rtseis/enums.h"
#include <complex>
#include <vector>
#include <memory>

namespace RTSeis
{
namespace Utils
{
namespace Filters
{
    class Precision
    {
        public:
            Precision(const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE) :
                precision_(precision)
            {
                return;
            }
            ~Precision(void){return;}
            Precision& operator=(const Precision &precision)
            {
                if (&precision == this){return *this;}
                precision_ = precision.precision_;
                return *this;
            }
            int setPrecision(const enum rtseisPrecision_enum precision)
            {
                precision_ = precision;
                return 0;
            }
            enum rtseisPrecision_enum getPrecision(void) const{return precision_;}
            bool isDoublePrecision(void) const
            {
                if (precision_ == RTSEIS_DOUBLE){return true;}
                return false;
            }
            bool isFloatPrecision(void) const
            {
                if (precision_ == RTSEIS_FLOAT){return true;}
                return false;
            }
            void print(FILE *fptr)
            {
                FILE *f = stdout;
                if (fptr != nullptr){f = fptr;}
                if (isDoublePrecision())
                {
                    fprintf(f, "Class is double precision\n");
                }
                else
                {
                    fprintf(f, "Class is float precision\n"); 
                }
            }
        private: 
            enum rtseisPrecision_enum precision_ = RTSEIS_DOUBLE;
    };

    class RealTime
    {
        public:
            RealTime(const bool lrt = false) :
                lrt_(lrt)
            {
                return;
            }
            ~RealTime(void){return;}
            RealTime &operator=(const RealTime realTime)
            {
                if (&realTime == this){return *this;}
                lrt_ = realTime.lrt_;
                return *this;
            }
            void toggleRealTime(const bool lrt)
            {
                lrt_ = lrt;
                return;
            }
            bool isRealTime(void) const {return lrt_;}
            void print(FILE *fptr)
            {
                FILE *f = stdout;
                if (fptr != nullptr){f = fptr;}
                if (isRealTime())
                {
                    fprintf(f, "Class is for real-time processing\n");
                }
                else
                {
                    fprintf(f, "Class is for post-processing\n"); 
                }
            }
        private:
            bool lrt_ = false;
    };

    class Downsample : protected Precision, RealTime
    {
        public:
            Downsample(void);
            ~Downsample(void);
            Downsample(const Downsample &downsample);
            Downsample& operator=(const Downsample &downsample);
            int initialize(const int downFactor,
                           const bool lisRealTime = false,
                           const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE);
            int estimateSpace(const int n) const;
            int getDownsampleFactor(void) const;
            int setInitialConditions(const int phase);
            int apply(const int n, const double x[],
                      const int ny, int *nyDown, double y[]);
            int apply(const int n, const float x[],
                      const int ny, int *nyDown, float y[]);
            int resetInitialConditions(void);
            void clear(void);
        private:
            /*!< Initial conditions for phase. */
            int phase0_ = 0;
            /*!< Downsampling factor. */
            int downFactor_ = 0;
            /*!< The phase. */
            int phase_ = 0;
            /*!< Flag indicating the module is initialized. */
            bool linit_ = false;
    };

    class MedianFilter : protected Precision, RealTime
    {
        public:
            MedianFilter(void);
            ~MedianFilter(void);
            MedianFilter(const MedianFilter &median);
            MedianFilter& operator=(const MedianFilter &median);
            int initialize(const int n,
                           const bool lisRealTime = false,
                           const enum rtseisPrecision_enum precision = RTSEIS_DOUBLE);
            int getInitialConditionLength(void) const;
            int getGroupDelay(void) const;
            int setInitialConditions(const int nz, const double zi[]);
            int apply(const int n, const double x[], double y[]);
            int apply(const int n, const float x[], float y[]);
            int resetInitialConditions(void);
            void clear(void);
        private:
            /*!< Delay line source vector.  This has dimension [nwork_]. */
            void *dlysrc_ = nullptr;
            /*!< Delay line destination.  This has dimension [nwork_]. */
            void *dlydst_ = nullptr;
            /*!< Workspace for median filter.  This has dimension [bufferSize_]. */
            void *pBuf_ = nullptr;
            /*!< A reference of the saved initial conditions.  This has 
                 dimension [nwork_] though only the first maskSize_  - 1
                 points are valid. */
            double *zi_ = nullptr;
            /*!< The median filter window length. */
            int maskSize_ = 0;
            /*!< The workspace for the temporary arrays. */
            int nwork_ = 0;
            /*!< The size of the workspace buffer. */
            int bufferSize_ = 0;
            /*!< Flag indicating the module is initialized. */
            bool linit_ = false;
    };

    class FIRFilter
    {
        enum class Implementation
        {
            DIRECT = 0, /*!< Direct-form implementation. */
            FFT = 1     /*!< FFT overlap and add implementation. */ 
        };
        public:
            FIRFilter(void);
            ~FIRFilter(void);
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
}; /* End Filters. */
}; /* End Utils. */
}; /* End RTSeis. */

#endif /* FILTERS */
