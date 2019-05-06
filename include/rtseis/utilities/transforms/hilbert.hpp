#ifndef RTSEIS_UTILITIES_TRANSFORMS_HILBERT_HPP
#define RTSEIS_UTILITIES_TRANSFORMS_HILBERT_HPP 1
#include <memory>
#include <complex>
#include "rtseis/enums.h"
#include "rtseis/utilities/transforms/enums.hpp"

namespace RTSeis
{
namespace Utilities
{
namespace Transforms
{

/*!
 * @class Hilbert hilbert.hpp "include/rtseis/utilities/transforms/hilbert.hpp"
 * @brief Computes the Hilbert transform of a signal.
 * @ingroup rtseis_utils_transforms_dft
 */
class Hilbert
{
public:
    /*! @name Constructors
     * @{ 
     */
    /*!
     * @brief Default constructor.
     */
    Hilbert();
    /*!
     * @brief Copy constructor.
     * @param[in] dftr2c  Class from which to initialize from.
     */
    Hilbert(const Hilbert &hilbert);
    /*! @} */

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy operator.
     * @param[in] dftr2c  DFTRealToComplex class to copy.
     * @result A deep copy of the input class.
     */
    Hilbert& operator=(const Hilbert &hilbert);
    /*! @} */

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
    void initialize(const int n,
                    const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);
    /*!
     * @brief Applies the Hilbert transform to a signal.
     * @param[in] n   The number of samples in the signal to transform.
     *                This must equal \c getTransformLength().
     * @param[in] x   The signal to Hilbert transform.  This is an array
     *                of dimension [n].
     * @param[out] h  The Hilbert transformed signal.  This is an array
     *                of dimension [n].
     * @throws std::invalid_argument n is incorrect or if x or h are NULL.
     * @throws std::runtime_error if the class is not initialized.
     */
    void transform(const int n, const double x[], std::complex<double> h[]);
    //void transform(const int n, const float x[], std::complex<float> h[]); 
private:
    class HilbertImpl;
    std::unique_ptr<HilbertImpl> pImpl;
};

} // End Transforms
} // End Utilities
} // End RTSeis

#endif
