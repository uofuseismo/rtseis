#ifndef RTSEIS_MODULES_CLASSICSTALTA_HPP
#define RTSEIS_MODULES_CLASSICSTALTA_HPP 1
#include "rtseis/config.h"
#include "rtseis/enums.h"
#include "rtseis/utils/filters.hpp"

namespace RTSeis
{
namespace Modules
{

/*!
 * @defgroup rtseis_modules_cSTALTA_parameters Parameters
 * @brief Defines the parameters for classic STA/LTA module.
 * @ingroup rtseis_modules_cSTALTA
 * @copyright Ben Baker distributed under the MIT license.
 */
class ClassicSTALTAParameters
{
    public:
        /*!
         * @brief Default constructor.  This module will not yet be usable
         *        until the parameters are set.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */ 
        ClassicSTALTAParameters(void);
        /*!
         * @brief Copy constructor.
         * @param[in] parameters  Class from which to initialize.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        ClassicSTALTAParameters(const ClassicSTALTAParameters &parameters);
        /*!
         * @brief Copy operator.
         * @param[in] parameters  Class to copy.
         * @result A deep copy of the STA/LTA parameter class.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        ClassicSTALTAParameters& operator=(const ClassicSTALTAParameters &parameters);
        /*!
         * @brief Initializes the Classic STA/LTA parameters.
         * @param[in] nsta  Number of samples in the short-term average
         *                  window.  This must be positive.
         * @param[in] nlta  Number of samples in the long-term average
         *                  window.  This must be greater than nlta.
         * @param[in] mode  Indicates whether or not this is for real-time.
         *                  By default this is for post-processing.
         * @param[in] precision  Defines the precision.  By default this
         *                       is a double precision module.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        ClassicSTALTAParameters(
            const int nsta, const int nlta,
            const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
            const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        /*!
         * @brief Initializes the Classic STA/LTA parameters.
         * @param[in] staWin  The short-term average window duration in 
         *                    seconds.  This must be non-negative.
         * @param[in] ltaWin  The long-term average window duration in 
         *                    seconds.  This must be greater than the STA
         *                    window length plus the half sampling period
         *                    i.e., \f$ LTA \gt STA + \frac{\Delta T}{2} \f$.
         * @param[in] dt      The sampling period in seconds.  This must be
         *                    positive.
         * @param[in] mode    Indicates whether or not this is for real-time.
         *                    By default this is for post-processing.
         * @param[in] precision  Defines the precision.  By default this
         *                       is a double precision module.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        ClassicSTALTAParameters(
            const double staWin, const double ltaWin,
            const double dt,
            const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
            const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        /*!
         * @brief Initializes the Classic STA/LTA parameters.
         * @param[in] nsta  Number of samples in the short-term average
         *                  window.  This must be positive.
         * @param[in] nlta  Number of samples in the long-term average
         *                  window.  This must be greater than nlta.
         * @param[in] chunkSize  A tuning parameter that defines the temporary
         *                       storage of the workspace arrays.  This should
         *                       be a power of 2 and positive.
         * @param[in] mode  Indicates whether or not this is for real-time.
         *                  By default this is for post-processing.
         * @param[in] precision  Defines the precision.  By default this
         *                       is a double precision module.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        ClassicSTALTAParameters(
            const int nsta, const int nlta,
            const size_t chunkSize,
            const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
            const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        /*!
         * @brief Initializes the Classic STA/LTA parameters.
         * @param[in] staWin  The short-term average window duration in 
         *                    seconds.  This must be non-negative.
         * @param[in] ltaWin  The long-term average window duration in 
         *                    seconds.  This must be greater than the STA
         *                    window length plus half the sampling period
         *                    i.e., \f$ LTA \gt STA + \frac{\Delta T}{2} \f$.
         * @param[in] chunkSize  A tuning parameter that defines the temporary
         *                       storage of the workspace arrays.  This should
         *                       be a power of 2 and positive. 
         * @param[in] dt   The sampling period in seconds.  This must be
         *                 positive.
         * @param[in] mode  Indicates whether or not this is for real-time.
         *                  By default this is for post-processing.
         * @param[in] precision  Defines the precision.  By default this
         *                       is a double precision module.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        ClassicSTALTAParameters(
            const double staWin, const double ltaWin,
            const double dt, 
            const size_t chunkSize,
            const RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
            const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
        /*!
         * @brief Default destructor.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */ 
        ~ClassicSTALTAParameters(void);
        /*!
         * @brief Clears variables in class and restores defaults.
         *        This class will have to be re-initialized to use again.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        virtual void clear(void);
        /*!
         * @brief Determines if the class parameters are valid and can be
         *        used to initialize the STA/LTA processing.
         * @retval True indicates that the parameters are valid.
         * @retval False indicates that the parameters are invalid.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        bool isValid(void) const;
        /*!
         * @brief Sets the chunksize.  The real-time FIR filter will be applied
         *        by looping over chunks of data.  There are three temporary
         *        storage arrays.  This parameter can help minimize the memory
         *        footprint of the module.
         * @param[in] chunkSize  The length of the scratch space arrays.
         * @result 0 indicates success.
         * @ingroup rtseis_modules_cSTALTA_parameters 
         */
        int setChunkSize(const size_t chunkSize);
        /*!
         * @brief Gets the chunksize.
         * @result The chunk size for the temporary arrays to minimize 
         *         temporary space overhead.
         * @ingroup rtseis_modules_cSTALTA_parameters 
         */
        size_t getChunkSize(void) const;
        /*!
         * @brief Sets the short-term and long-term window size in samples.
         * @param[in] nsta  Number of samples in the short-term average
         *                  window.  This must be positive.
         * @param[in] nlta  Number of samples in the long-term average
         *                  window.  This must be greater than nlta.
         * @result 0 indicates success.   On failure the number of samples
         *         in the short-term and long-term windows will remain 
         *         unchaged.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        int setShortTermAndLongTermWindowSize(const int nsta, const int nlta);
        /*!
         * @brief Sets the short-term and long-term window durations.
         * @param[in] staWin  The short-term average window duration in 
         *                    seconds.  This must be non-negative.
         * @param[in] ltaWin  The long-term average window duration in 
         *                    seconds.  This must be greater than the STA
         *                    window length plus the sampling period
         *                    i.e., \f$ LTA \gt STA + \frac{\Delta T}{2} \f$.
         * @param[in] dt      The sampling period in seconds.  This must be
         *                    positive.
         * @result 0 indicates success.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        int setShortTermAndLongTermWindowSize(const double staWin,
                                              const double ltaWin,
                                              const double dt);
        /*!
         * @brief Gets the number of samples in the long-term window.
         * @result The number of samples in the long-term window.
         *         If this is negative then an error has occurred.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        int getLongTermWindowSize(void) const;
        /*!
         * @brief Gets the number of samples in teh short-term window.
         * @result The number of samples in the short-term window.
         *         If this is negative then an error has occurred.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        int getShortTermWindowSize(void) const;
        /*!
         * @brief Enables the class as being for real-time application or not.
         * @param[in] mode  Indicates whether the module is for post-processing
         *                  or real-time processing.
         * @ingroup rtseis_modules_onebit_parameters
         */
        void setProcessingMode(const RTSeis::ProcessingMode mode);
        /*!
         * @brief Determines if the class is for real-time application.
         * @retval True indicates that the class is for real-time
         *         application.
         * @retval False indicates that the class is not for real-time
         *         application.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        RTSeis::ProcessingMode getProcessingMode(void) const;
        /*!
         * @brief Determines the precision of the class.
         * @result The precision with which the underlying copysign
         *         operation will be performed.
         * @ingroup rtseis_modules_cSTALTA_parameters
         */
        RTSeis::Precision getPrecision(void) const;
    private:
        /*!< Routine to validate the parameters. */
        void validate_(void);
        /*!< Default precision. */
        const RTSeis::Precision defaultPrecision_ = RTSeis::Precision::DOUBLE;
        /*!< The number of samples in the short term average window. */
        int nsta_ = 0;
        /*!< The number of samples in the long-term average window. */
        int nlta_ = 0;
        /*!< A tuning parameter that controls the chunk size.  The 
             filter will require temporary space. */
        size_t chunkSize_ = 1024;
        /*!< The precision of the module. */
        RTSeis::Precision precision_ = defaultPrecision_;
        /*!< Flag indicating this module is for real-time or post-processing. */
        RTSeis::ProcessingMode processingMode_ = RTSeis::ProcessingMode::POST_PROCESSING;
        /*!< Flag indicating that this is a valid module for processing. */
        bool isValid_ = false;
};

/*!
 * @defgroup rtseis_modules_cSTALTA Classic STA/LTA
 * @brief Computes the `classic' Short Term to Long Term Average using
 *        an FIR filter.
 * @ingroup rtseis_modules
 * @copyright Ben Baker distributed under the MIT license.
 */
class ClassicSTALTA : ClassicSTALTAParameters
{
     public:
        /*!
         * @brief Default constructor.  This module will not yet be usable
         *        until the parameters are set.
         * @ingroup rtseis_modules_cSTALTA
         */
        ClassicSTALTA(void);
        /*!
         * @brief Initializes the classic STA/LTA from the parameters.
         * @param[in] parameters  Parameters from which to initialize
         *                        the classic STA/LTA.
         * @ingroup rtseis_modules_cSTALTA
         */
        ClassicSTALTA(const ClassicSTALTAParameters parameters);
        /*!
         * @brief Copy constructor.
         * @param[in] cstalta  A Classic STA/LTA class from which this
         *                     class is initialized.
         * @ingroup rtseis_modules_cSTALTA
         */ 
        ClassicSTALTA(const ClassicSTALTA &cstalta);
        /*!
         * @brief Copy operator.
         * @param[in] cstalta  A classic STA/LTA class to copy.
         * @result A deep copy of the input class.
         * @ingroup rtseis_modules_cSTALTA
         */
        ClassicSTALTA& operator=(const ClassicSTALTA &cstalta);
        /*!
         * @brief Default destructor.  
         * @ingroup rtseis_modules_cSTALTA
         */
        ~ClassicSTALTA(void);
        /*!
         * @brief Returns the number of coefficients in the numerator initial
         *        conditions array.
         * @result The number of elements in the numerator initial condition
         *         array.  If negative then an error has occured. 
         * @ingroup rtseis_modules_cSTALTA
         */
        int getNumeratorInitialConditionLength(void) const;
        /*!
         * @brief Returns the number of coefficients in the denominator initial
         *        conditions array.
         * @result The number of elements in the denominator initial condition
         *         array.  If negative then an error has occured. 
         * @ingroup rtseis_modules_cSTALTA
         */
        int getDenominatorInitialConditionLength(void) const;
        /*!
         * @brief Sets the initial conditions on the filter.
         * @param[in] nzNum  The length of the numerator initial conditions
         *                   array.  This
         * @param[in] zNum   The numerator initial condition coefficients.
         *                   This has dimension [nzNum].
         * @param[in] nzDen  The length of the denominator initial conditions
         *                   array.
         * @param[in] zDen   The denominator initial condition coefficients.
         *                   This has dimension [nzDen].
         */
        int setInitialConditions(const int nzNum, const double zNum[],
                                 const int nzDen, const double zDen[]);
        /*!
         * @brief Computes the STA/LTA of the input signal.  This will reset
         *        the initial conditions prior to setting the new initial
         *        conditions.
         * @param[in] nx   Number of points in signal.
         * @param[in] x    The signal of which to compute the STA/LTA.  This has
         *                 dimension [nx].
         * @param[in] y    The STA/LTA signal.  This has dimension [nx].
         * @retval  0 indicates success.
         * @ingroup rtseis_modules_cSTALTA
         */ 
        int apply(const int nx, const double x[], double y[]);
        /*! 
         * @brief Computes the STA/LTA of the input signal.
         * @param[in] nx   Number of points in signal.
         * @param[in] x    The signal of which to compute the STA/LTA.  This has
         *                 dimension [nx].
         * @param[in] y    The STA/LTA signal.  This has dimension [nx].
         * @result 0 indicates success.
         * @ingroup rtseis_modules_cSTALTA
         */ 
        int apply(const int nx, const float x[], float y[]);
        /*!
         * @brief Resets the filter to the initial conditions specified
         *        by setInitialConditions() or the default initial conditions. 
         * @result 0 indicates success.
         */
        int resetInitialConditions(void);
        /*!
         * @brief Clears variables in class and restores defaults.
         *        This class will have to be re-initialized to use again.
         * @ingroup rtseis_modules_cSTALTA
         */
        void clear() override;
        /*!
         * @brief Determines if the class is for real-time application.
         * @retval If true then the class is for real-time application.
         * @retval If false then the class is not for real-time application.
         */
        bool isInitialized(void) const;
    private:
        /*!< Numerator FIR signal to keep track of the short-term average. */
        RTSeis::Utilities::Filters::FIRFilter firNum_;
        /*!< Denominator FIR signal to keep track of the long-term average. */
        RTSeis::Utilities::Filters::FIRFilter firDen_;
        /*!< Workspace array for holding input signal squared. */
        void *x2_ = nullptr;
        /*!< Workspace array for holding the numerator. */
        void *ynum_ = nullptr;
        /*!< Workspace array for holding the denominator. */
        void *yden_ = nullptr; 
        /*!< The STA/LTA parameters. */
        ClassicSTALTAParameters parms_;
        /*!< Flag indicating the module is intialized. */
        bool isInitialized_ = false;
        
};

};
};

#endif
