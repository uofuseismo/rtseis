#ifndef RTSEIS_UTILITIES_PICKERS_CLASSICSTALTA_HPP
#define RTSEIS_UTILITIES_PICKERS_CLASSICSTALTA_HPP 1
#include <memory>
#include "rtseis/enums.hpp"

namespace RTSeis
{
namespace Modules
{

/*!
 * @class ClassicSTALTAParameters classicSTALTAParameters.hpp "include/rtseis/utilities/pickers/classicSTALTAParameters.hpp"
 * @brief Defines the parameters for the classic STA/LTA module.
 * @ingroup rtseis_utils_pickers
 * @copyright Ben Baker distributed under the MIT license.
 */
class ClassicSTALTAParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.  This module will not yet be usable
     *        until the parameters are set.
     */ 
   ClassicSTALTAParameters();
    /*!
     * @brief Copy constructor.
     * @param[in] parameters  Class from which to initialize.
     */
    ClassicSTALTAParameters(const ClassicSTALTAParameters &parameters);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy operator.
     * @param[in] parameters  Class to copy.
     * @result A deep copy of the STA/LTA parameter class.
     */
    ClassicSTALTAParameters& operator=(const ClassicSTALTAParameters &parameters);
    /*! @} */

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
     */
    ClassicSTALTAParameters(int nsta, int nlta,
                            RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                            RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Initializes the Classic STA/LTA parameters.
     * @param[in] staWin  The short-term average window duration in 
     *                    seconds.  This must be non-negative.
     * @param[in] ltaWin  The long-term average window duration in 
     *                    seconds.  This must be greater than the STA
     *                    window length plus the half sampling period
     *                    i.e., \f$ LTA > STA + \frac{\Delta T}{2} \f$.
     * @param[in] dt      The sampling period in seconds.  This must be
     *                    positive.
     * @param[in] mode    Indicates whether or not this is for real-time.
     *                    By default this is for post-processing.
     * @param[in] precision  Defines the precision.  By default this
     *                       is a double precision module.
     */
    ClassicSTALTAParameters(double staWin, double ltaWin,
                            double dt,
                            RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                            RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
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
     */
    ClassicSTALTAParameters(int nsta, int nlta,
                            size_t chunkSize,
                            RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                            RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Initializes the Classic STA/LTA parameters.
     * @param[in] staWin  The short-term average window duration in 
     *                    seconds.  This must be non-negative.
     * @param[in] ltaWin  The long-term average window duration in 
     *                    seconds.  This must be greater than the STA
     *                    window length plus half the sampling period
     *                    i.e., \f$ LTA > STA + \frac{\Delta T}{2} \f$.
     * @param[in] chunkSize  A tuning parameter that defines the temporary
     *                       storage of the workspace arrays.  This should
     *                       be a power of 2 and positive. 
     * @param[in] dt   The sampling period in seconds.  This must be
     *                 positive.
     * @param[in] mode  Indicates whether or not this is for real-time.
     *                  By default this is for post-processing.
     * @param[in] precision  Defines the precision.  By default this
     *                       is a double precision module.
     */
    ClassicSTALTAParameters(double staWin, double ltaWin, double dt, 
                            size_t chunkSize,
                            RTSeis::ProcessingMode mode = RTSeis::ProcessingMode::POST_PROCESSING,
                            RTSeis::Precision precision = RTSeis::Precision::DOUBLE);

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */ 
    ~ClassicSTALTAParameters();
    /*!
     * @brief Clears variables in class and restores defaults.
     *        This class will have to be re-initialized to use again.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Determines if the class parameters are valid and can be
     *        used to initialize the STA/LTA processing.
     * @retval True indicates that the parameters are valid.
     * @retval False indicates that the parameters are invalid.
     */
    bool isValid() const noexcept;
    /*!
     * @brief Sets the chunksize.  The real-time FIR filter will be applied
     *        by looping over chunks of data.  There are three temporary
     *        storage arrays.  This parameter can help minimize the memory
     *        footprint of the module.
     * @param[in] chunkSize  The length of the scratch space arrays.
     * @throws std::invalid_argument if chunkSize is not positive.
     * @result 0 indicates success.
     */
    int setChunkSize(size_t chunkSize);
    /*!
     * @brief Gets the chunksize.
     * @result The chunk size for the temporary arrays to minimize 
     *         temporary space overhead.
     */
    size_t getChunkSize() const;
    /*!
     * @brief Sets the short-term and long-term window size in samples.
     * @param[in] nsta  Number of samples in the short-term average
     *                  window.  This must be positive.
     * @param[in] nlta  Number of samples in the long-term average
     *                  window.  This must be greater than nlta.
     * @result 0 indicates success.   On failure the number of samples
     *         in the short-term and long-term windows will remain 
     *         unchaged.
     */
    int setShortTermAndLongTermWindowSize(int nsta, int nlta);
    /*!
     * @brief Sets the short-term and long-term window durations.
     * @param[in] staWin  The short-term average window duration in 
     *                    seconds.  This must be non-negative.
     * @param[in] ltaWin  The long-term average window duration in 
     *                    seconds.  This must be greater than the STA
     *                    window length plus the sampling period
     *                    i.e., \f$ LTA > STA + \frac{\Delta T}{2} \f$.
     * @param[in] dt      The sampling period in seconds.  This must be
     *                    positive.
     * @result 0 indicates success.
     */
    int setShortTermAndLongTermWindowSize(double staWin,
                                          double ltaWin,
                                          double dt);
    /*!
     * @brief Gets the number of samples in the long-term window.
     * @result The number of samples in the long-term window.
     *         If this is negative then an error has occurred.
     */
    int getLongTermWindowSize() const;
    /*!
     * @brief Gets the number of samples in teh short-term window.
     * @result The number of samples in the short-term window.
     *         If this is negative then an error has occurred.
     */
    int getShortTermWindowSize() const;
    /*!
     * @brief Enables the class as being for real-time application or not.
     * @param[in] mode  Indicates whether the module is for post-processing
     *                  or real-time processing.
     */
    void setProcessingMode(RTSeis::ProcessingMode mode) noexcept;
    /*!
     * @brief Determines if the class is for real-time application.
     * @retval The processing mode which can be post-processing or real-time.
     */
    RTSeis::ProcessingMode getProcessingMode() const;
    /*!
     * @brief Determines the precision of the class.
     * @result The precision with which the underlying copysign
     *         operation will be performed.
     */
    RTSeis::Precision getPrecision() const;

private:
    /*!< Routine to validate the parameters. */
    void validate_();
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
class ClassicSTALTA
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
        ClassicSTALTA(const ClassicSTALTAParameters &parameters);
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
        int apply(const int nx, const double x[], double *y[]);
        /*! 
         * @brief Computes the STA/LTA of the input signal.
         * @param[in] nx   Number of points in signal.
         * @param[in] x    The signal of which to compute the STA/LTA.  This has
         *                 dimension [nx].
         * @param[in] y    The STA/LTA signal.  This has dimension [nx].
         * @result 0 indicates success.
         * @ingroup rtseis_modules_cSTALTA
         */ 
        int apply(const int nx, const float x[], float *y[]);
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
        void clear(void);
        /*!
         * @brief Determines if the class is for real-time application.
         * @retval If true then the class is for real-time application.
         * @retval If false then the class is not for real-time application.
         */
        bool isInitialized(void) const;
    private:
        class ClassicSTALTAImpl;
        std::unique_ptr<ClassicSTALTAImpl> pSTALTA_;
};

};
};

#endif
