#ifndef RTSEIS_TRANSFORMS_HILBERT_HPP
#define RTSEIS_TRANSFORMS_HILBERT_HPP 1
#include <memory>
#include <complex>
#include "rtseis/transforms/enums.hpp"

namespace RTSeis::Transforms
{
/// @class Hilbert hilbert.hpp "include/rtseis/utilities/transforms/hilbert.hpp"
/// @brief Computes the analytic signal.
/// @note The Hilbert transform is obtained from the analytic signal, 
///       \f$ x_a = x_r + i x_h \f$ where \f$ x_r \f$ is the input signal and
///       \f$ x_h \f$ is the Hilbert transform of the signal.
/// @ingroup rtseis_utils_transforms
template<class T = double>
class Hilbert
{
public:
    /// @name Constructors
    /// @{ 
    /// @brief Default constructor.
    Hilbert();
    /// @brief Copy constructor.
    /// @param[in] hilbert  The class from which to initialize this class.
    Hilbert(const Hilbert &hilbert);
    /// @brief Move constructor.
    /// @param[in,out] hilbert  The class from which to initialize this class.
    ///                         On exit, hilbert's behavior is undefined.
    Hilbert(Hilbert &&hilbert) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] hilbert  The hilbert class to copy.
    /// @result A deep copy of the input class.
    Hilbert& operator=(const Hilbert &hilbert);
    /// @brief Move assignment operator.
    /// @param[in,out] hilbert  The hilbert class whose memory will be moved
    ///                         to this.  On exit, hilbert's behavior is
    ///                         undefined.
    /// @result The memory from hilbert moved to this.
    Hilbert& operator=(Hilbert &&hilbert) noexcept;
    /// @}

    /*! @name Destructors
     * @{
     */
    /*!
      * @brief Destructor.
      */
    ~Hilbert();
    /*!
     * @brief Resets the module and releases all memory.
     */
    void clear() noexcept;

    /*!
     * @brief Checks if the class is initialized.
     * @result True indicates that the class is initialized.
     */
    bool isInitialized() noexcept;
    /*!
     * @brief Gets the length of the Hilbert transform.
     * @result The length of the Hilbert transform.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getTransformLength();
    /*!
     * @brief Initializes the Hilbert transformer.
     * @param[in] n   The number of samples in the signal to transform.
     * @param[in] precision  Defines the underlying precision of the Hilbert
     *                       transform calculation.
     * @throws std::invalid_argument if n is not positive.
     */
    void initialize(int n);
    /*!
     * @brief Computes the analytic signal of the signal.
     * @param[in] n   The number of samples in the signal to transform.
     *                This must equal \c getTransformLength().
     * @param[in] x   The signal of which to compute the analytic signal.
     *                This is an array of dimension [n].
     * @param[out] h  The analytic signal.  This is an array of dimension [n].
     * @throws std::invalid_argument n is incorrect or if x or h are NULL.
     * @throws std::runtime_error if the class is not initialized.
     * @note The real part of \f$ h \f$ is the input signal and the imaginary
     *       part of \f$ h \f$ is the Hilbert transform.
     */
    void transform(int n, const T x[], std::complex<T> *h[]);
private:
    class HilbertImpl;
    std::unique_ptr<HilbertImpl> pImpl;
};
}
#endif
