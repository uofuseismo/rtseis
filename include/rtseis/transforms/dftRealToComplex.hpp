#ifndef RTSEIS_TRANSFORMS_DFTREALTOCOMPLEX_HPP
#define RTSEIS_TRANSFORMS_DFTREALTOCOPMLEX_HPP 1
#include <memory>
#include <complex>
#include "rtseis/transforms/enums.hpp"
#include "rtseis/enums.hpp"

namespace RTSeis::Transforms 
{
/// @class DFTRealToComplex dftRealToComplex.hpp "includertseis/utilities/transforms/dftRealToComplex.hpp" 
/// @brief Handles the discrete Fourier transform of real signals.  This
///        leverages an efficiency that requires that the value of a negative
///        frequency is a complex conjugate of the corresponding positive
///        frequency.
/// @ingroup rtseis_utils_transforms
template<class T = double>
class DFTRealToComplex
{
public:
    /// @name Constructors
    /// @{ 
    /// @brief Default constructor.
    DFTRealToComplex();
    /// @brief Copy constructor.
    /// @param[in] dftr2c  DFT class from which to initialize this class.
    DFTRealToComplex(const DFTRealToComplex &dftr2c); 
    /// @brief Move constructor.
    /// @param[in,out] dftr2c  DFT class from which to initialize this class. 
    ///                        On exit, dftr2c's behavior is undefined.
    DFTRealToComplex(DFTRealToComplex &&dftr2c) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy operator.
    /// @param[in] dftr2c  DFTRealToComplex class to copy to this.
    /// @result A deep copy of the input class.
    DFTRealToComplex& operator=(const DFTRealToComplex &dftr2c);
    /// @brief Move operator.
    /// @param[in,out] dftr2c  The DFTRealToComplex class whose memory will be
    ///                        to this.  On exit, dftr2c's behavior is undefined.
    /// @result The memory from dftr2c moved to this.
    DFTRealToComplex& operator=(DFTRealToComplex &&dftr2c) noexcept; 
    /// @}

    /// @brief Initializes the real-to-complex Fourier transform class.
    /// @param[in] length          The length of the signal to transform.  
    ///                            This must be greater than or equal to 2.
    ///                            If this is a power of 2 then the FFT will
    ///                            be used instead of the DFT.
    ///                            Note, length can be an overestimate.  
    ///                            For example, if you have a signal of 
    ///                            100 and you select length to be 128 then
    ///                            on applying the transform the trailing 
    ///                            28 samples of the input signal will be
    ///                            zero-padded.
    /// @param[in] implementation  Defines the Fourier transform
    ///                            implementation.
    /// @param[in] precision       Controls the precision with which the DFT
    ///                            will be applied.
    /// @result 0 indicates success.
    void initialize(int length,
                    FourierTransformImplementation implementation = FourierTransformImplementation::DFT);
    /// @brief Fourier transforms a real-valued time domain signal to the
    ///        frequency domain.   The transform is defined as 
    ///        \f$ y(\omega) = \int e^{-i \omega t} x(t) dt \f$.
    /// @param[in] n     Number of points in input signal.  This cannot
    ///                  be negative and cannot exceed
    ///                  getMaximumSignalLength().  If n is less than
    ///                  getMaximumSignalLength() then it will be
    ///                  zero-padded prior to transforming.
    /// @param[in] x     Time domain signal to transform.  This has
    ///                  dimension [n].
    /// @param[in] maxy  The maximum number of points allocated to y.
    ///                  This should be at least getTransformLength().
    /// @param[out] y    The Fourier transformed signal.  This has
    ///                  dimension [maxy] though only the first 
    ///                  getTransformLength() points are defined. 
    ///                  Here, y[0] is the zero-frequency.  When 
    ///                  \c getMaximumSignalLength() is even 
    ///                  y[getTransformLength()-1] is the Nyquist.
    ///                  Otherwise, y[getTransformLength()-1] is
    ///                  (\c getMaximumSignalLength() - 1 )
    ///                 /(samplingPeriod*\c getMaximumSignalLength()).
    /// @throws std::invalid_argument if any arguments are invalid.
    /// @throws std::runtime_error if the class is not initialized.
    void forwardTransform(int n, const T x[],
                          int maxy, std::complex<T> *y[]);
    /// @brief Inverse transforms complex-valued frequency domian signal
    ///        to a real-valued time domain signal.  The transform is
    ///        defined as
    ///        \f$ y(t) = \frac{1}{n}
    ///                   \int e^{ i \omega t} x(\omega) d \omega \f$.
    ///        where \f$ n \f$ is the length of the transform given
    ///        in getTransformLength().
    /// 
    /// @param[in] lenft   Length of the Fourier transformed signal.
    ///                    This cannot be negative and cannot exceed
    ///                    getTransformLength().  If lenft is less 
    ///                    than getTransformLength() then it will be
    ///                    zero-padded prior to transforming.
    /// @param[in] x       The complex signal to inverse transform.
    ///                    This has dimension [lenft].  Here, x[0] is
    ///                    the zero-frequency and x[lenft-1] is the
    ///                    Nyquist frequency.
    /// @param[in] maxy    The maximum number of points allocated to y.
    ///                    This should be at least
    ///                    getMaximumInputSignalLength().
    ///                    If it exceeds this length, then elements
    ///                    beyond getMaximumSignalLength() will have
    ///                    unpredictable values.
    /// @param[out] y      The inverse-transformed time-domain signal.
    ///                    This has dimension [maxy] however only
    ///                    the first getMaximumSignalLength() points
    ///                    are defined.
    /// @throws std::invalid_argument if any arguments are invalid.
    /// @throws std::runtime_error if the class is not initialized.
    void inverseTransform(int lenft, const std::complex<T> x[],
                          int maxy, T *y[]);
    /// @result The length of the inverse DFT or FFT.
    /// @throws std::invalid_argument if the class is not initialized.
    [[nodiscard]] int getInverseTransformLength() const;
    /// @result The length of the DFT or FFT.
    /// @throws std::runtime_error if the class is not initialized.
    [[nodiscard]] int getTransformLength() const;
    /// @result The maximum length of the input signal.
    /// @throws std::invalid_argument if the class is not initialized.
    [[nodiscard]] int getMaximumInputSignalLength() const;
    /// @result True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;

    /// @name Destructors
    /// @{
    /// @brief Releases the memory on the module and resets the
    ///        defaults.  The class must be reinitialized before
    ///        using it again.
    /// @ingroup rtseis_utils_transforms_dftr2c
    void clear() noexcept;
    /// @brief Default destructor.
    ~DFTRealToComplex();
    /// @}
private:
    class DFTImpl;
    std::unique_ptr<DFTImpl> pImpl;
}; // End DFTRealToComplex
}
#endif
