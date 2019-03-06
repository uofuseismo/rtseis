#ifndef RTSEIS_POSTPROCESSING_SC_WAVEFORM
#define RTSEIS_POSTPROCESSING_SC_WAVEFORM 1
#include <memory>
#include <string>
#include <exception>
#ifndef RTSEIS_UTILS_MATH_CONVOLVE_HPP
#include "rtseis/utilities/math/convolve.hpp"
#endif
#ifndef RTSEIS_POSTPROCESSING_SC_TAPER
#include "rtseis/postProcessing/singleChannel/taper.hpp"
#endif

namespace RTSeis
{
/*!
 * @defgroup rtseis_postprocessing Post-Processing
 * @brief These are higher level algorithms that expedite post-processing
 *        of seismic data. 
 */
namespace PostProcessing
{
/*!
 * @defgroup rtseis_postprocessing_sc Single-Channel Processing
 * @brief This section contains single-channel post-processing algorithms.
 * @ingroup rtseis_postprocessing
 */
namespace SingleChannel
{
/*!
 * @class Waveform Waveform "include/rtseis/processing/singleChannel/postProcessing.hpp"
 * @brief This class is used for single-channel post-processing applications.
 * @ingroup rtseis_postprocessing_sc
 */
class Waveform : public std::exception
{
    public:
        /*!
         * @brief Default constructor.
         */
        Waveform(void);
        /*!
         * @brief Constructs a waveform from time series data.  The sampling
         *        period will be unity.
         * @param[in] x  Signal from which to construct time series.
         * @throws std::invalid_argument if x is empty.
         */
        explicit Waveform(const std::vector<double> &x);
        /*!
         * @brief Constructs a waveform from time series data.
         * @param[in] dt   Sampling period in seconds.
         * @param[in] x    Signal from which to construct time series.
         * @throws std::invalid_argument if the sampling period is not positive
         *         or x is empty.
         */
        Waveform(const double dt, const std::vector<double> &x);
        /*!
         * @brief Default destructor.
         */ 
        ~Waveform(void);
        /*!
         * @brief Sets a signal on waveform class.
         * @throws std::invalid_argument If x is empty.
         */
        void setData(const std::vector<double> &x);
        /*!
         * @brief Sets a waveform on the module.
         * @param[in] n   The number of points in the signal.
         * @param[in] x   The signal to set.  This is an array
         *                of dimension [n].
         * @throws std::invalid_argument if x is null or n is less than 1.
         */
        void setData(const size_t n, const double x[]);
        /*!
         * @brief Gets the processed waveform data.
         * @param[out] y  The processed waveform data.
         */
        void getData(std::vector<double> &y);
        /*!
         * @brief Computes the convolution \f$ x \ast s \f$ where the
         *        convolution sum is defined by
         *        \f$ y[k] = \sum_n x[n] s[n-k] \f$.
         * @param[in] s     The signal to convolve with x.
         * @param[in] mode  Defines the convolution output.
         * @param[in] implementation  Defines the implementation type.
         * @throws std::invalid_argument if s is empty or there is no data.
         */
        void convolve(const std::vector<double> &s,
            const RTSeis::Utilities::Math::Convolve::Mode mode = RTSeis::Utilities::Math::Convolve::Mode::FULL,
            const RTSeis::Utilities::Math::Convolve::Implementation implementation = RTSeis::Utilities::Math::Convolve::Implementation::AUTO);
        /*!
         * @brief Removes the mean from the data.
         * @throws std::invalid_argument if there is no data.
         *
         * @snippet testing/postProcessing/singleChannel.cpp ppSCDemeanExample
         */
        void demean(void);
        /*!
         * @brief Removes a best fitting line \f$ \hat{y} = a x + b \f$
         *        from the data by first computing \f$ a \f$ and \f$ b \f$
         *        then computing \f$ y - (a x + b) \f$.
         * @throws std::invalid_argument if there are less than 2 data points 
         *         in x.
         *
         * @snippet testing/postProcessing/singleChannel.cpp ppSCDetrendExample
         */
        void detrend(void);
        /*!
         * @brief Tapers the ends of a signal.
         * @param[in] pct  The percentage of the signal to which the taper will
         *                 be applied.  For example, 5 percent indicates that
         *                 the first 2.5 and final 2.5 percent of the signal
         *                  will be tapered. 
         * @param[in] window  Defines the window function used to generate
         *                    the taper.
         * @note The SAC tapering convention is used where the first and final
         *       samples of the tapered signal are set to 0.
         * @throw std::invalid_argument if the parameters are invalid.
         */
        void taper(const double pct = 5,
                   const TaperParameters::Type window = TaperParameters::Type::HAMMING);
    private:
        class DataImpl;
        std::unique_ptr<DataImpl> pData_;
}; // end waveform
}; // End SingleChannel
}; // End PostProcessing
}; // End RTSeis

# endif
