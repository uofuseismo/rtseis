#ifndef RTSEIS_UTILS_FILTERS_HPP
#define RTSEIS_UTILS_FILTERS_HPP 1
#include <memory>
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
     * @defgroup rtseis_utils_filters_downsample Downsample
     * @brief This is the core implementation for downsampling a signal.
     * @copyright Ben Baker distributed under the MIT license.
     * @ingroup rtseis_utils_filters
     */
    class Downsample
    {
        public:
            /*!
             * @brief Default constructor.
             */
            Downsample(void);
            /*!
             * @brief Copy constructor.
             * @param[in] downsample  Downsampling class from which
             *                        to initialize.
             */
            Downsample(const Downsample &downsample);
            /*!
             * @brief Copy operator
             * @param[in] downsample  Downsampling class to copy.
             * @result A deep copy of the downsampling class.
             */
            Downsample& operator=(const Downsample &downsample);
            /*!
             * @brief Default destructor.
             */
            ~Downsample(void);
            /*!
             * @brief Initializes the downsampler.
             * @param[in] downFactor  The downsampling factor.  This must be
             *                        positive.  This will retain every 
             *                        (downFactor-1)'th sample.
             * @param[in] mode  The processing mode.  By default this
             *                  is for post-processing.
             * @param[in] precision   The precision of the filter.  By default
             *                        this is double precision.
             * @result 0 indicates success.
             */
            int initialize(const int downFactor,
                           const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
            /*!
             * @brief Determines if the module is initialized.
             * @retval True indicates that the module is initialized.
             * @retval False indicates that the module is not initialized.
             */
            bool isInitialized(void) const;
            /*!
             * @brief Estimates the space required to hold the downsampled
             *        signal.
             * @param[in] n   The length of the signal to downsample.
             * @result The number of points required to store the output signal.
             *         If negative then there was a failure.
             */
            int estimateSpace(const int n) const;
            /*!
             * @brief Gets the downsampling factor.
             * @result The downsampling factor.
             */
            int getDownsampleFactor(void) const;
            /*!
             * @brief Sets the initial conditions of the downsampler which is
             *        the phase.
             * @param[in] phase  Phase of downsampler.  This must be in the 
             *                   range [0, getDownsampleFactor()].
             * @result 0 indicates success.
             */
            int setInitialConditions(const int phase);
            /*!
             * @brief Applies the downsampler to the data.
             * @param[in] nx       The number data points in x.
             * @param[in] x        The signal to downsample.
             * @param[in] ny       The maximum number of samples in y.  One can
             *                     estimate ny by using estimateSpace(). 
             * @param[out] nyDown  The number of defined downsampled points
             *                     in y.
             * @param[out] y       The downsampled signal.  This has dimension
             *                     [ny] however  only the first [nyDown] points
             *                     are defined.
             * @result 0 indicates success.
             */
            int apply(const int n, const double x[],
                      const int ny, int *nyDown, double y[]);
            int apply(const int n, const float x[],
                      const int ny, int *nyDown, float y[]);
            /*!
             * @brief Resets the initial conditions to the phase set in 
             *        setInitialConditions.  If setInitialConditions was not
             *        called then this will set the phase to 0.
             * @result 0 indicates success.
             */
            int resetInitialConditions(void);
            /*! 
             * @brief Clears the module and resets all parameters.
             */
            void clear(void);
        private:
            /* Forward of class for PIMPL. */
            class DownsampleImpl;
            /* Pointer to the the filter implementation and parameters. */
            std::unique_ptr<DownsampleImpl> pDownsample_; 
    };

    /*!
     * @defgroup rtseis_utils_filters_median Median Filter
     * @brief This is the core implementation for median filtering.
     * @copyright Ben Baker distributed under the MIT license.
     * @ingroup rtseis_utils_filters
     */
    class MedianFilter
    {
        public:
            /*!
             * @brief Copy constructor.
             */
            MedianFilter(void);
            /*!
             * @brief Destructor.
             */
            ~MedianFilter(void);
            /*!
             * @brief Copy constructor.
             * @param[in] median  Median class from which to initialize.
             */
            MedianFilter(const MedianFilter &median);
            /*!
             * @brief Copy operator.
             * @param[in] median  The median class to copy.
             * @result A deep copy of the median filter class.
             */
            MedianFilter& operator=(const MedianFilter &median);
            /*!
             * @brief Initializes the median filter.
             * @param[in] n   The window size of the median filter.  This must
             *                be a positive and odd number.  If n is not odd
             *                then it's length will be increased by 1.
             * @param[in] mode  The processing mode.  By default this
             *                  is for post-processing.
             * @param[in] precision   The precision of the filter.  By default
             *                        this is double precision.
             * @result 0 indicates success.
             */
            int initialize(const int n,
                           const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
            /*!
             * @brief Determines if the module is initialized.
             * @retval True indicates that the module is initialized.
             * @retval False indicates that the module is not initialized.
             */
            bool isInitialized(void) const;
            /*!
             * @brief Utility routine to determine the initial condition length.
             * @retval A non-negative number is the length of the initial
             *         condition array.
             * @retval -1 Indicates failure.
             */
            int getInitialConditionLength(void) const;
            /*!
             * @brief Returns the group delay of the filter.  Note, that this
             *        shift is required to get a correspondence to Matlab.
             * @result The group delay.
             */
            int getGroupDelay(void) const;
            /*!
             * @brief Sets the initial conditions for the filter.  This should
             *        be called prior to filter application as it will reset
             *        the filter.
             * @param[in] nz   The median filter initial conditions.  This
             *                 should be equal to getInitialConditionLength().
             * @param[in] zi   The initial conditions.  This has dimension [nz].
             * @result 0 indicates success.
             */
            int setInitialConditions(const int nz, const double zi[]);
            /*!
             * @brief Appplies the median filter to the array x.
             * @param[in] n   Number of points in x.
             * @param[in] x   The signal to filter.  This has dimension [n].
             * @param[out] y  The filtered signal.  This has dimension [n].
             * @result 0 indicates success.
             */
            int apply(const int n, const double x[], double y[]);
            int apply(const int n, const float x[], float y[]);
            /*!
             * @brief Resets the initial conditions on the source delay line to
             *        the default initial conditions or the initial conditions
             *        set when MedianFilter::setInitialConditions() was called.
             * @result 0 indicates success.
             */ 
            int resetInitialConditions(void);
            /*! 
             * @brief Clears the module and resets all parameters.
             */
            void clear(void);
        private:
            /* Forward declaration of class for PIMPL. */
            class MedianFilterImpl;
            /* Pointer to the the filter implementation and parameters. */
            std::unique_ptr<MedianFilterImpl> pMedian_;
    };

    /*!
     * @defgroup rtseis_utils_filters_sos Second Order Sections
     * @brief This is the core implementation for second order section (biquad)
     *        infinite impulse response filtering.
     * @copyright Ben Baker distributed under the MIT license.
     * @ingroup rtseis_utils_filters
     */
    class SOSFilter
    {
        public:
            /*!
             * @brief Default constructor.
             */
            SOSFilter(void);
            /*!
             * @brief A copy constructor.
             * @param[in] sos  The SOS class from which to initialize.
             */
            SOSFilter(const SOSFilter &sos);
            /*!
             * @brief Copy operator.
             * @param[in] sos  The class to copy.
             * @result A deep copy of the input SOS class.
             */
            SOSFilter& operator=(const SOSFilter &sos);
            /*!
             * @brief Default destructor.
             */
            ~SOSFilter(void);
            /*!
             * @brief Initializes the second order section filter.
             * @param[in] ns    The number of second order sections.
             * @param[in] bs    Numerator coefficients.  This is an array of
             *                  dimension [3 x ns] with leading dimension 3.
             *                  There is a further requirement that b[3*is]
             *                  for \f$ i_s=0,1,\cdots,n_s-1 \f$ not be zero.
             * @param[in] as    Denominator coefficients.  This is an array of
             *                  dimension [3 x ns] with leading dimension 3. 
             *                  There is a further requirement that a[3*is]
             *                  for \f$ i_s=0,1,\cdots,n_s-1 \f$ not be zero.
             * @param[in] mode  The processing mode.  By default this
             *                  is for post-processing.
             * @param[in] precision   The precision of the filter.  By default
             *                        this is double precision.
             * @result 0 indicates success.
             */
            int initialize(const int ns,
                           const double bs[],
                           const double as[],
                           const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
            /*!
             * @brief Determines if the module is initialized.
             * @retval True indicates that the module is initialized.
             * @retval False indicates that the module is not initialized.
             */
            bool isInitialized(void) const;
            /*!
             * @brief Returns the length of the initial conditions.
             * @result The length of the initial condtions array.
             */
            int getInitialConditionLength(void) const;
            /*!
             * @brief Sets the initial conditions for the filter.  This should
             *        be called prior to filter application as it will reset
             *        the filter.
             * @param[in] nz   The second order section filter initial
             *                 conditions.  This should be equal to
             *                 getInitialConditionLength().
             * @param[in] zi   The initial conditions.  This has dimension [nz].
             * @result 0 indicates success.
             */
            int setInitialConditions(const int nz, const double zi[]);
            /*!
             * @brief Applies the second order section filter to the data.
             * @param[in] n   Number of points in signals.
             * @param[in] x   The signal to filter.  This has dimension [n].
             * @param[out] y  The filtered signal.  This has dimension [n].
             */
            int apply(const int n, const double x[], double y[]);
            int apply(const int n, const float x[], float y[]);
            /*!
             * @brief Resets the initial conditions on the source delay line
             *        to the default initial conditions or the initial
             *        conditions set when SOSFilter::setInitialConditions()
             *        was called.
             */
            int resetInitialConditions(void);
            /*! 
             * @brief Clears the module and resets all parameters.
             */
            void clear(void);
            int getNumberOfSections(void) const;
        private:
            /* Forward declaration of class for PIMPL. */
            class SOSFilterImpl;
            /* Pointer to the the filter implementation and parameters. */
            std::unique_ptr<SOSFilterImpl> pSOS_;
    };

    /*!
     * @defgroup rtseis_utils_filters_fir FIR Filtering
     * @brief This is the core implementation for FIR filtering.
     * @copyright Ben Baker distributed under the MIT license.
     * @ingroup rtseis_utils_filters
     */
    class FIRFilter
    {
        public:
            /*!
             * @brief Defines the implementation of the FIR filter.
             */
            enum Implementation
            {
                /*!< Direct-form implementation. */
                DIRECT = 0,
                /*!< FFT overlap and add implementation. */
                FFT = 1,
                /*!< The implementation will decide between 
                     DIRECT or FFT. */
                AUTO = 2
            };
        public:
            /*!
             * @brief Default constructor.
             */
            FIRFilter(void);
            /*!
             * @brief Copy constructor.
             * @param[in] fir   FIR class from which to initialize.
             */
            FIRFilter(const FIRFilter &fir);
            /*!
             * @brief Copy operator.
             * @param[in] fir   FIR class to copy.
             * @result A deep copy of the FIR class.
             */
            FIRFilter& operator=(const FIRFilter &fir);
            /*!
             * @brief Default destructor.
             */
            ~FIRFilter(void);
            /*!
             * @brief Initializes the FIR filter.
             * @param[in] nb    Number of numerator coefficients.
             * @param[in] b     Numerator coefficients.  This is an array of
             *                  dimension [nb].
             * @param[in] mode  The processing mode.  By default this
             *                  is for post-processing.
             * @param[in] precision   The precision of the filter.  By default
             *                        this is double precision.
             * @param[in] implementation  Defines the implementation.  This can
             *                            specify Implementation::DIRECT form
             *                            for direct form implementation,
             *                            Implementation::FFT for an FFT-based
             *                            overlap-add implementation which 
             *                            can be advantageous when
             *                            \f$ \log_2 L < N \f$ where \f$ L \f$
             *                            is the signal length and \f N \f$
             *                            the number of filter taps, or, auto
             *                            to let the computer decide.
             *                            The default is to use the direct form.
             * @result 0 indicates success.
             */
            int initialize(const int nb, const double b[],
                           const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                           const RTSeis::Precision precision = RTSeis::Precision::DOUBLE,
                           Implementation implementation = Implementation::DIRECT);
            /*!
             * @brief Determines if the module is initialized.
             * @retval True indicates that the module is initialized.
             * @retval False indicates that the module is not initialized.
             */
            bool isInitialized(void) const;
            /*!
             * @brief Utility routine to determine the initial condition length.
             * @result A non-negative number is the length of the initial
             *         condition array.
             */
            int getInitialConditionLength(void) const;
            /*!
             * @brief Sets the initial conditions for the filter.  This should
             *        be called prior to filter application as it will reset
             *        the filter.
             * @param[in] nz   The FIR filter initial condition length.
             *                 This should be equal to
             *                 getInitialConditionLength().
             * @param[in] zi   The initial conditions.  This has dimension [nz].
             * @result 0 indicates success.
             */
            int setInitialConditions(const int nz, const double zi[]);
            /*!
             * @brief Gets a copy of the initial conditions.
             * @param[in] nz   The FIR filter initial condition length.
             *                 This should be equal to
             *                 getInitialConditionLength().
             * @param[out] zi  The initial conditions.  This has dimension [nz].
             * @result 0 indicate success.
             */
            int getInitialConditions(const int nz, double zi[]) const;
            /*!
             * @brief Applies the FIR filter to the data.
             * @param[in] n   Number of points in signals.
             * @param[in] x   Signal to filter.  This has dimension [n].
             * @param[out] y  The filtered signal.  This has dimension [n].
             */
            int apply(const int n, const double x[], double y[]);
            int apply(const int n, const float x[], float y[]);
            /*!
             * @brief Resets the initial conditions on the source delay line to
             *        the default initial conditions or the initial conditions
             *        set when FIRFilter::setInitialConditions() was called.
             * @result 0 indicates success.
             */
            int resetInitialConditions(void);
            /*!
             * @brief Clears the module and resets all parameters.
             */
            void clear(void);
        private:
            /* Forward of class for PIMPL. */
            class FIRImpl;
            /* Pointer to the the filter implementation and parameters. */
            std::unique_ptr<FIRImpl> pFIR_;
    };

    /*!
     * @defgroup rtseis_utils_filters_mrfir Multirate FIR Filtering
     * @brief Implements the multi-rate finite impulse response filters. 
     *        This allows for upsampling and/or downsampling while
     *        filtering.  Note, that this behaves slightly differently
     *        than Matlab.  When the upsample factor is greater than 1
     *        then the FIR filter will be implicitly multiplied 
     *        by the upsampling factor.  Compare this with Matlab where
     *        the user must gain the FIR filter prior to calling upfirdn.
     * @ingroup rtseis_utils_filters
     */
    class MultiRateFIRFilter
    {
        public:
            /*!
             * @brief Default constructor.  Note, this class is not yet
             *        inititalized and cannot be used.
             */
            MultiRateFIRFilter(void);
            /*!
             * @brief Copy constructor.
             * @param[in] firmr  Multi-rate filter class from which to 
             *                   initialize this class.
             */
            MultiRateFIRFilter(const MultiRateFIRFilter &firmr);
            /*!
             * @brief Copy operator.
             * @param[in] firmr  The multi-rate filter class to copy.
             * @result A deep copy of hte multi-rate filter.
             */
            MultiRateFIRFilter &operator=(const MultiRateFIRFilter &firmr);
            /*!
             * @brief Default destructor.
             */
            ~MultiRateFIRFilter(void);
            /*!
             * @brief Initializes the multi-rate filtering.
             * @param[in] upFactor    The upsampling factor.  This will insert
             *                        upFactor - 1 zeros to the signal prior to
             *                        filtering.  This must be positive.  There
             *                        will be no upsampling effect if upFactor
             *                        is 1.
             * @param[in] downFactor  The downsampling factor.  This will remove
             *                        remove downFactor - 1 samples beteween
             *                        each upsampled output point.  This must be
             *                        positive.   There will be no downsampling
             *                        effect if downFactor is 1.
             * @param[in] nb          The number of FIR filter coefficients.
             * @param[in] b           The FIR filter coefficients.  This is
             *                        an array of dimension [nb].
             * @param[in] chunkSize   This is an optional tuning parameter.
             *                        The internal workspace size will be
             *                        max(upFactor, downFactor)*chunkSize.
             * @param[in] mode        The processing mode.  By default this
             *                        is for post-processing.
             * @param[in] precision   The precision of the filter.  By default
             *                        this is double precision.
             * @result 0 indicates success.
             */
            int initialize(const int upFactor, const int downFactor,
                           const int nb, const double b[],
                           const RTSeis::ProcessingMode mode
                               = RTSeis::ProcessingMode::POST_PROCESSING,
                           const RTSeis::Precision precision
                               = RTSeis::Precision::DOUBLE);
            int initialize(const int upFactor, const int downFactor,
                           const int nb, const double b[],
                           const int chunkSize,
                           const RTSeis::ProcessingMode mode 
                               = RTSeis::ProcessingMode::POST_PROCESSING,
                           const RTSeis::Precision precision 
                               = RTSeis::Precision::DOUBLE);
            /*!
             * @brief Determines if the module is initialized.
             * @retval True indicates that the module is initialized.
             * @retval False indicates that the module is not initialized.
             */
            bool isInitialized(void) const;
            /*!
             * @brief Gets the length of the initial conditions array.
             * @result On successful exit this is the length of the initial
             *         conditions array.  On failure this is negative.
             */
            int getInitialConditionLength(void) const;
            /*!
             * @brief Applies the zero-phase IIR filter to the data.  Note,
             *        the class must be initialized prior to using this
             *        function.
             * @param[in] nz  The length of the initial conditions array.  This
             *                can be determined by first calling
             *                getInitialConditionLength().
             * @param[in] zi  The initial conditions.  This has dimension [nz].
             * @param[in] upPhase    The initial upsampling phase.  By default
             *                       this is 0.
             * @param[in] downPhase  The initial downsampler phase.  By default
             *                       this is 0.
             * @result 0 indicates success.
             */
            int setInitialConditions(const int nz,
                                     const double zi[],
                                     const int upPhase = 0,
                                     const int downPhase = 0);
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
             * @brief Estimates the space required to store the output signal.
             * @param[in] n  The length of the input signal.
             * @result The array length required to store the output signal.
             *         If negative than error has occurred.
             */
            int estimateSpace(const int n) const;
            /*!
             * @brief Applies the multi-rate filter to the signal.
             * @param[in] n       Number of points in the signal.
             * @param[in] x       The signal to filter.  This has dimension [n].
             * @param[in] nywork  The workspace reserved to y.  An estimate 
             *                    can be obtained from estimateSpace().
             * @param[out] ny     The length of the output signal.
             * @param[out] y      The filtered signal.  This has dimension
             *                    [nywork] however only the first [ny] samples
             *                    are set.
             * @result 0 indicates success.
             */
            int apply(const int n, const double x[],
                      const int nywork, int *ny, double y[]);
            /*!< @copydoc apply */
            int apply(const int n, const float x[],
                      const int nywork, int *ny, float y[]);
            /*!
             * @brief Resets the initial conditions to those set in
             *        setInitialConditions().
             * @result 0 indicates success.
             */
            int resetInitialConditions(void);
            /*!
             * @brief Clears memory and resets the filter.  After applying
             *        this function the filter must be re-initialized prior
             *        to being applied to the data.
             */
            void clear(void);
        private:
            /* Forward of class for PIMPL. */
            class MultiRateFIRImpl;
            /* Pointer to the the filter implementation and parameters. */
            std::unique_ptr<MultiRateFIRImpl> pFIR_;
             
    };

    /*!
     * @defgroup rtseis_utils_filters_iiriir  Zero-Phase IIR Filtering
     * @brief Implements a zero-phase IIR filter.  This is for
     *        post-processing only.
     * @ingroup rtseis_utils_filters
     * @copyright Ben Baker distributed under the MIT license.
     */
    class IIRIIRFilter
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
             * @brief Determines if the module is initialized.
             * @retval True indicates that the module is initialized.
             * @retval False indicates that the module is not initialized.
             */
            bool isInitialized(void) const;
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
            /* Forward of class for PIMPL. */
            class IIRIIRImpl;
            /* Pointer to the the filter implementation and parameters. */
            std::unique_ptr<IIRIIRImpl> pIIRIIR_;
    };
}; /* End Filters. */
}; /* End Utils. */
}; /* End RTSeis. */

#endif /* FILTERS */
