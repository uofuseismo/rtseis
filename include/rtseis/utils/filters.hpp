#ifndef RTSEIS_UTILS_FILTERS_HPP
#define RTSEIS_UTILS_FILTERS_HPP 1
#include <stdio.h>
#include "rtseis/config.h"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utils
{

/*!
 * @defgroup rtseis_utils_filters Filter Implementations
 * @brief These are the core filter implementations to be used by higher
 *        level modules.
 * @copyright Ben Baker distributed under the MIT license.
 */    
namespace Filters
{
    /*!
     * @brief Defines the precision of the filter implementation.  In general,
     *        both float and double precision filters are implemented.
     * @ingroup rtseis_utils_filters
     */
    class Precision
    {
        public:
            /*!
             * @brief Default constructor.
             * @param[in] precision  The implementation precision.  By default 
             *                       the filters are applied in double
             *                       precision.
             */
            Precision(
                const RTSeis::Precision precision = RTSeis::Precision::DOUBLE) :
                precision_(precision)
            {
                return;
            }
            /*!
             * @brief Copy constructor.
             * @param[in] precision  Precision class from which to initialize.
             */
            Precision(const Precision &precision)
            {
                *this = precision;
            }
            /*!
             * @brief Default destructor.
             */
            ~Precision(void){return;}
            /*!
             * @brief Copy operator.
             * @param[in] precision  Precision class to copy to this class.
             * @result A deep copy of the input class.
             */
            Precision& operator=(const Precision &precision)
            {
                if (&precision == this){return *this;}
                precision_ = precision.precision_;
                return *this;
            }
            /*!
             * @brief Sets the filter precision.
             * @param[in] precision  The precision which the filter will
             *                       be applied.
             * @result 0 indicates success.
             */
            int setPrecision(const RTSeis::Precision precision)
            {
                precision_ = precision;
                return 0;
            }
            /*!
             * @brief Gets the filter precision.
             * @result The precision of the filter to apply.
             */
            RTSeis::Precision getPrecision(void) const{return precision_;}
            /*! 
             * @brief Determines if the filter is for double precision.
             * @retval True indicates a double floating point arithmetic filter
             *         should be applied.
             * @retval False indicates a double floating point arithmetic filter
             *         should not be applied.
             */
            bool isDoublePrecision(void) const
            {
                if (precision_ == RTSeis::Precision::DOUBLE){return true;}
                return false;
            }
            /*!
             * @brief Determines if the filter is for single precision.
             * @retval True indicates a single flaoting point arithmetic filter
             *         should be applied.
             * @retval False indicates a single floating point arithmetic filter
             *         should not be applied.
             */
            bool isFloatPrecision(void) const
            {
                if (precision_ == RTSeis::Precision::FLOAT){return true;}
                return false;
            }
            /*!
             * @brief Prints the class contents to file.
             * @param[in] fptr  File handle to which the contents of the class
             *                  are printed.  By default this is stdout.
             */ 
            void print(FILE *fptr = stdout)
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
            /*! The precision of the module. */
            RTSeis::Precision precision_ = RTSeis::Precision::DOUBLE;
    };

    /*!
     * @brief Class to help a filter distinguish whether it should be applied
     *        in a real-time or post-processing fashion.
     * @ingroup rtseis_utils_filters
     */
    class RealTime
    {
        public:
            /*!
             * @brief Default constructor.
             * @param[in] lrt   True indicates that this filter is for real-time
             *                  application.
             * @param[in] lrt   False indicates that this filter is for
             *                  post-processing.  This is the default.
             */
            RealTime(const bool lrt = false) :
                lrt_(lrt)
            {
                return;
            }
            /*!
             * @brief Default destructor.
             */
            ~RealTime(void){return;}
            /*!
             * @brief Copy operator.
             * @param[in] realTime  Real-time class to copy.
             * @result A deep copy of the class.
             */
            RealTime &operator=(const RealTime realTime)
            {
                if (&realTime == this){return *this;}
                lrt_ = realTime.lrt_;
                return *this;
            }
            /*!
             * @brief Enables/disables the class as being for real-time
             *        processing.
             * @param[in] lrt   True indicates that this filter is for real-time
             *                  application.
             * @param[in] lrt   False indicates that this filter is for
             *                  post-processing.
             */
            void toggleRealTime(const bool lrt)
            {
                lrt_ = lrt;
                return;
            }
            /*!
             * @brief Returns if the class is for real-time processing.
             * @retval True indicates that the class is for real-time
             *         processing.
             * @retval False indicates that the class is for post-processing.
             */
            bool isRealTime(void) const {return lrt_;}
            /*!
             * @brief Prints the contents of the class to a file.
             * @param[in] fptr   File pointer to which this class will print its
             *                   contents.  The default is stdout.
             */
            void print(FILE *fptr = stdout)
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
            /*! Flag indicating whether or not this is for real-time. */
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
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
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
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
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
            /*!< The workspace for the delay lines. */
            int nwork_ = 0;
            /*!< The size of the workspace buffer. */
            int bufferSize_ = 0;
            /*!< Flag indicating the module is initialized. */
            bool linit_ = false;
    };

    class SOSFilter : protected Precision, RealTime
    {
        public:
            SOSFilter(void);
            ~SOSFilter(void);
            SOSFilter(const SOSFilter &sos);
            SOSFilter& operator=(const SOSFilter &sos);
            int initialize(const int ns,
                           const double bs[],
                           const double as[],
                           const bool lisRealTime = false,
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
            int getInitialConditionLength(void) const;
            int setInitialConditions(const int nz, const double zi[]);
            int apply(const int n, const double x[], double y[]);
            int apply(const int n, const float x[], float y[]);
            int resetInitialConditions(void);
            void clear(void);
            int getNumberOfSections(void) const;
        private:
            /*!< The IIR filtering state. */
            void *pState_ = nullptr;
            /*!< A workspace buffer.  This has dimension [bufferSize_]. */
            void *pBuf_ = nullptr;
            /*!< The filter coefficients.  This has dimension [tapsLen_]. */
            void *pTaps_ = nullptr;
            /*!< Delay line source vector.  This has dimension [nwork_]. */
            void *dlysrc_ = nullptr;
            /*!< Delay line destination vector.  This has dimension [nwork_]. */
            void *dlydst_ = nullptr;
            /*!< A copy of the input filter numerator coefficients.  This
                 has dimension [3 x nsections_]. */
            double *bsRef_ = nullptr;
            /*!< A copy of the input filter denominator coefficients.  This
                 has dimension [3 x nsections_]. */
            double *asRef_ = nullptr;
            /*!< A copy of the initial conditions.  This has dimension
                 [2*nsections_]. */
            double *zi_ = nullptr;
            /*!< The number of sections. */
            int nsections_ = 0;
            /*!< The number of filter taps. */
            int tapsLen_ = 0;
            /*!< The workspace for the delay lines. */
            int nwork_ = 0; 
            /*!< The size of the workspace buffer. */
            int bufferSize_ = 0;
            /*!< Flag indicating the module is initialized. */
            bool linit_ = false;
    };

    class FIRFilter : protected Precision, RealTime
    {
        public:
            /*!
             * @brief Defines the implementation of the FIR filter.
             */
            enum Implementation
            {
                DIRECT = 0, /*!< Direct-form implementation. */
                FFT = 1     /*!< FFT overlap and add implementation. */ 
            };
        public:
            FIRFilter(void);
            ~FIRFilter(void);
            FIRFilter(const FIRFilter &fir);
            FIRFilter& operator=(const FIRFilter &fir);
            int initialize(const int nb, const double b[],
                           const bool lisRealTime = false,
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE,
                           Implementation implementation = Implementation::DIRECT);
            int getInitialConditionLength(void) const;
            int getInitialConditions(const int nz, double zi[]) const;
            int setInitialConditions(const int nz, const double zi[]);
            int apply(const int n, const double x[], double y[]);
            int apply(const int n, const float x[], float y[]);
            int resetInitialConditions(void);
            void clear(void);
        private:
            /*!< The filter state. */
            void *pFIRSpec_ = nullptr;
            /*!< The filter taps.  This has dimension [tapsLen_]. */
            void *pTaps_ = nullptr;
            /*!< The input delay line.  This has dimension [nwork_]. */
            void *dlysrc_ = nullptr;
            /*!< The output delay line.  This has dimension [nwork_]. */
            void *dlydst_ = nullptr;
            /*!< Workspace.  This has dimension [bufferSize_]. */
            void *pBuf_ = nullptr;
            /*!< A copy of the input taps. */
            double *tapsRef_ = nullptr;
            /*!< A copy of the initial conditions.  This has dimension
                 [order]. */
            double *zi_ = nullptr;
            /*!< The number of taps. */
            int tapsLen_ = 0; 
            /*!< The length of the delay line which is max(128, order+1). */
            int nwork_ = 0;
            /*!< Size of pBuf. */
            int bufferSize_ = 0;
            /*!< Size of state. */
            int specSize_ = 0;
            /*!< Filter order. */
            int order_ = 0;
            /*!< Implementation. */
            Implementation implementation_ = Implementation::DIRECT; 
            /*!< Flag indicating the module is initialized. */
            bool linit_ = false;
    };

    /*!
     * @defgroup rtseis_utils_filters_iiriir  Zero-Phase IIR Filtering
     * @brief Implements a zero-phase IIR filter.  This is for
     *        post-processing only.
     * @ingroup rtseis_utils_filters
     * @copyright Ben Baker distributed under the MIT license.
     */
    class IIRIIRFilter : protected Precision, RealTime
    {
        public:
            /*!
             * @brief Default constructor.  Note, this class is not yet 
             *        initialized and cannot be used.
             */
            IIRIIRFilter(void);
            /*!
             * @brief Default destructor.
             */
            ~IIRIIRFilter(void);
            /*!
             * @brief A copy constructor.
             * @param[in] iiriir  The IIRIIRFilter class to initialize from.
             */
            IIRIIRFilter(const IIRIIRFilter &iiriir);
            /*!
             * @brief A copy operator.
             * @param[in] iiriir  The IIRIIRFilter class to copy.
             * @result A deep copy of the iiriir filter class.
             */
            IIRIIRFilter &operator=(const IIRIIRFilter &iiriir);
            /*!
             * @brief Initializes the zero-phase IIR filter.
             * @param[in] nb   The number of numerator coefficients.  This must
             *                 be positive.
             * @param[in] b    The numerator (feedforward) coefficients.  This
             *                 is an array of dimension [nb].
             * @param[in] na   The number of denominator coefficients.  This
             *                 must be positive.
             * @param[in] a    The denominator (feedback) coefficients.  This is
             *                 an array of dimension [na].
             * @param[in] precision  The precision of the filter.  By default
             *                       this is double precision.
             * @result 0 indicates success.
             */
            int initialize(const int nb, const double b[],
                           const int na, const double a[],
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
            /*!
             * @brief Gets the length of the initial conditions array.
             * @result On successful exit this is the length of the initial
             *         conditions array.  On failure this is negative.
             */
            int getInitialConditionLength(void) const;
            /*!
             * @brief Sets the initial conditions.
             * @param[in] nz   The length of the initial conditions.  This
             *                 should equal getInitialConditionLength().
             * @param[in] zi   The initial conditions to set.  This has
             *                 has dimension [nwork_]. 
             * @result 0 indicates success.
             */
            int setInitialConditions(const int nz, const double zi[]);
            /*!
             * @brief Applies the zero-phase IIR filter to the data.  Note,
             *        the class must be initialized prior to using this function.
             * @param[in] n   Number of points in signal.
             * @param[in] x   The signal to filter.  This has dimension [n].
             * @param[out] y  The zero-phase IIR filtered signal.  This has
             *                dimension [n].
             * @result 0 indicates success.
             */
            int apply(const int n, const double x[], double y[]);
            /*!
             * @brief Applies the zero-phase IIR filter to the data.  Note,
             *        the class must be initialized prior to using this function.
             * @param[in] n   Number of points in signal.
             * @param[in] x   The signal to filter.  This has dimension [n].
             * @param[out] y  The zero-phase IIR filtered signal.  This has
             *                dimension [n].
             * @result 0 indicates success.
             */
            int apply(const int n, const float x[], float y[]); 
            /*!
             * @brief Resets the initial conditions to those set in
             *        setInitialConditions().  Note, this will not do anything
             *        as the final conditions are never extracted from the
             *        IIR filter.
             * @result 0 indicates success.
             */
            int resetInitialConditions(void); 
            /*!
             * @brief Clears memory and resets the filter.  After applying
             *        this function the filter must be re-initialized prior
             *        to being applied to the data.
             */
            void clear(void);
            /*!
             * @brief Utilility function to get the filter order.
             * @result The filter order.
             */
            int getFilterOrder(void) const;
        private:
            /*!< A pointer to the IIRIIR filter state. */
            void *pState_ = nullptr;
            /*!< The IIR filter taps.  This has dimension [2*(order_+1)]. */
            void *pTaps_ = nullptr;
            /*!< Workspace for applying the IIR IIR filter.  This has 
                 dimension [bufferSize_]. */
            void *pBuf_ = nullptr;
            /*!< The initial conditions.  This has dimension [nwork_]. */
            void *dlysrc_ = nullptr;
            /*!< A saved copy of the numerator coefficients.  This has
                 dimension [nbRef_]. */
            double *bRef_ = nullptr;
            /*!< A saved copy of the denominator coefficients.  This has
                 dimension [naRef_]. */
            double *aRef_ = nullptr;
            /*!< A saved copy of the initial conditions.  This has
                 dimension [nwork_]. */
            double *zi_ = nullptr;
            /*!< The initial condtiions.  This has dimension [nwork_]. */
            int nwork_ = 0;
            /*!< The buffer size. */
            int bufferSize_ = 0;
            /*!< The filter order.  This equals max(nbRef_, naRef_) - 1. */
            int order_ = 0;
            /*!< The number of numerator coefficients. */
            int nbRef_ = 0;
            /*!< The number of denominator coefficients. */
            int naRef_ = 0;
            /*!< Flag indicating that the initial conditions have been set. */
            bool lhaveZI_ = false;
            /*!< Flag indicating that the filter is initialized. */
            bool isInitialized_ = false; 

    };
}; /* End Filters. */
}; /* End Utils. */
}; /* End RTSeis. */

#endif /* FILTERS */
