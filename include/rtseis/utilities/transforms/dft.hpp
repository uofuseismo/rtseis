#ifndef RTSEIS_UTILITIES_TRANSFORMS_DFT_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_DFT_HPP 1
#include <memory>
#include <complex>
#include "rtseis/utilities/transforms/enums.hpp"
#include "rtseis/enums.h"

namespace RTSeis
{
namespace Utilities
{
namespace Transforms 
{

/*!
 * @class DFT dft.hpp "includertseis/utilities/transforms/dft.hpp" 
 * @brief Handles the discrete Fourier transform of general complex signals.
 * @ingroup rtseis_utils_transforms_dft
 */
class DFT
{
public:
    /*! @name Constructors
     * @{ 
     */
    /*!
     * @brief Default constructor.
     */
    DFT(void);
    /*!
     * @brief Copy constructor.
     * @param[in] dftr2c  Class from which to initialize from.
     */
    DFT(const DFT &dft);
    /*!
     * @brief Copy operator.
     * @param[in] dftr2c  DFTRealToComplex class to copy.
     * @result A deep copy of the input class.
     */
    DFT& operator=(const DFT &dft);
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~DFT(void);
    /*! @} */

    /*!
      * @brief Initializes the real-to-complex Fourier transform class.
      * @param[in] length          The length of the signal to transform.  
      *                            This must be greater than or equal to 2.
      *                            If this is a power of 2 then the FFT will
      *                            be used instead of the DFT.
      *                            Note, length can be an overestimate.  
      *                            For example, if you have a signal of 
      *                            100 and you select length to be 128 then
      *                            on applying the transform the trailing 
      *                            28 samples of the input signal will be
      *                            zero-padded.
      * @param[in] implementation  Defines the Fourier transform
      *                            implementation.
      * @param[in] precision       Controls the precision with which the DFT
      *                            will be applied.
      * @throws std::invalid_argument if the arguments are invalid.
      */
     void initialize(const int length,
                     const FourierTransformImplementation implementation = FourierTransformImplementation::DFT,
                     const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
     /*!
      * @brief Fourier transforms a real-valued time domain signal to the
      *        frequency domain.   The transform is defined as 
      *        \f$ y(\omega) = \int e^{-i \omega t} x(t) dt \f$.
      * @param[in] n     Number of points in input signal.  This cannot
      *                  be negative and cannot exceed
      *                  getMaximumSignalLength().  If n is less than
      *                  getMaximumSignalLength() then it will be
      *                  zero-padded prior to transforming.
      * @param[in] x     Time domain signal to transform.  This has
      *                  dimension [n].
      * @param[in] maxy  The maximum number of points allocated to y.
      *                  This should be at least getTransformLength().
      * @param[out] y    The Fourier transformed signal.  This has
      *                  dimension [maxy] though only the first 
      *                  getTransformLength() points are defined. 
      *                  Here, y[0] is the zero-frequency and 
      *                  y[getTransformLength()-1] is the Nyquist.
      * @result 0 indicates success.
      * @throws std::invalid_argument if any arguments are invalid.
      * @throws std::runtime_error if the class is not initialized.
      */
     void forwardTransform(const int n,
                           const std::complex<double> x[],
                           const int maxy,
                           std::complex<double> y[]);
     /*! @copydoc forwardTransform */
     void forwardTransform(const int n,
                           const std::complex<float> x[],
                           const int maxy,
                           std::complex<float> y[]);
     /*!
      * @brief Inverse transforms complex-valued frequency domian signal
      *        to a real-valued time domain signal.  The transform is
      *        defined as
      *        \f$ y(t) = \frac{1}{n}
      *                   \int e^{ i \omega t} x(\omega) d \omega \f$.
      *        where \f$ n \f$ is the length of the transform given
      *        in getTransformLength().
      * 
      * @param[in] lenft   Length of the Fourier transformed signal.
      *                    This cannot be negative and cannot exceed
      *                    getTransformLength().  If lenft is less 
      *                    than getTransformLength() then it will be
      *                    zero-padded prior to transforming.
      * @param[in] x       The complex signal to inverse transform.
      *                    This has dimension [lenft].  Here, x[0] is
      *                    the zero-frequency and x[lenft-1] is the
      *                    Nyquist frequency.
      * @param[in] maxy    The maximum number of points allocated to y.
      *                    This should be at least
      *                    getMaximumInputSignalLength().
      *                    If it exceeds this length, then elements
      *                    beyond getMaximumSignalLength() will have
      *                    unpredictable values.
      * @param[out] y      The inverse-transformed time-domain signal.
      *                    This has dimension [maxy] however only
      *                    the first getMaximumSignalLength() points
      *                    are defined.
      * @throws std::invalid_argument if any arguments are invalid.
      * @throws std::runtime_error if the class is not initialized.
      */
     void inverseTransform(const int lenft,
                           const std::complex<double> x[],
                           const int maxy, std::complex<double> y[]);
     /*! @copydoc inverseTransform */
     void inverseTransform(const int lenft,
                           const std::complex<float> x[],
                           const int maxy, std::complex<float> y[]);
     /*!
      * @brief Gets the inverse transform length.
      * @result The length of the inverse DFT or FFT.
      * @throws std::runtime_error if the class is not initialized.
      */
     int getInverseTransformLength(void) const;
     /*!
      * @brief Gets the length of the transform.
      * @result The length of the DFT or FFT.
      * @throws std::runtime_error if the class is not initialized.
      */
     int getTransformLength(void) const;
     /*!
      * @brief Gets the maximum length of the input signal.
      * @result The maximum length of the input signal.
      * @throws std::runtime_error if the class is not initialized.
      */
     int getMaximumInputSignalLength(void) const;
     /*!
      * @brief Returns whether or not the class is initialized.
      * @retval True indicates the class is initialized.
      * @retval False indicates the class is not-initialized.
      */
    bool isInitialized(void) const noexcept;
    /*!
     * @brief Releases the memory on the module and resets the
     *        defaults.  The class must be reinitialized before
     *        using it again.
     */
    void clear(void) noexcept;
private:
    class DFTImpl;
    std::unique_ptr<DFTImpl> pImpl;
}; // End DFTRealToComplex

}; /* End transforms */
}; /* End utils */
}; /* End RTSeis */

#endif
