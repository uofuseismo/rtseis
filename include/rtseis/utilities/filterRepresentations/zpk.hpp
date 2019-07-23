#ifndef RTSEIS_UTILITIES_FILTERREPRESENTATIONS_ZPK_HPP
#define RTSEIS_UTILITIES_FILTERREPRESENTATIONS_ZPK_HPP 1
#include <complex>
#include <vector>
#include <memory>

namespace RTSeis::Utilities::FilterRepresentations
{
/*!
 * @class ZPK zpk.hpp "include/rtseis/utilities/filterRepresentations/zpk.hpp"
 * @brief Represents an infinite impulse response filter in terms of
 *        zeros, poles, and a gain.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_fr
 */
class ZPK
{
public:
    /*! @name Constructors.
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    ZPK(void);
    /*!
     * @brief Constructs a ZPK class from the given zeros, poles, and gain.
     * @param[in] zeros   Zeros to set.
     * @param[in] poles   Poles to set.
     * @param[in] k       Gain to set.
     */
    ZPK(const std::vector<std::complex<double>> &zeros,
        const std::vector<std::complex<double>> &poles,
        const double k);
    /*!
     * @brief Copy constructor.
     * @param[in] zpk  Class from which to initialize this class.
     */
    ZPK(const ZPK &zpk);
    /*!
     * @brief Move constructor.
     * @param[in,out] zpk  Class to move to this class.
     *                     On exit, zpk's behavior will be undefined.
     */
    ZPK(ZPK &&zpk);
    /*! @} */
    
    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] zpk   ZPK class to copy.
     * @result A deep copy of the input ZPK class.
     */ 
    ZPK &operator=(const ZPK &zpk);
    /*!
     * @brief Move operator.
     * @param[in,out] zpk  Class to move.  On exit zpk's behavior will be
     *                     undefined.
     * @result The moved version of zpk.
     */
    ZPK &operator=(ZPK &&zpk);
    /*!
     * @brief Equality operator.
     * @param[in] zpk  Class to compare to this class.
     * @result True indicates that zpk equals this class within a
     *         given tolerance.
     */
    bool operator==(const ZPK &zpk) const;
    /*!
     * @brief Inequality operator.
     * @param[in] zpk  Class to compare to this class.
     * @result True indicates that zpk does not equal this class within a
     *         given tolerance.
     */
    bool operator!=(const ZPK &zpk) const;
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~ZPK(void);
    /*! 
     * @brief Clears the class.
     */
    void clear(void);
    /*! @} */

    /*!
     * @brief Sorts the poles.
     * @param[in] ascending  If true then the poles will be sorted in
     *                       increasing order of magnitude.
     *                       If false then the poles will be sorted in
     *                       decreasing order of magnitude.
     */
    void sortPoles(bool ascending=true);
    /*!
     * @brief Sorts the zeros.
     * @param[in] ascending  If true then the zeros will be sorted in
     *                       increasing order of magnitude.
     *                       If false then the zeors will be sorted in
     *                       decreasing order of magnitude.
     */
    void sortZeros(bool ascending=true);
    /*!
     * @brief Prints the ZPK structure.
     * @param[in] fout   File handle to print to.  If fout is NULL then this
     *                   will print to stdout.
     */
    void print(FILE *fout = stdout) const noexcept;
    /*!
     * @brief Sets the gain.
     * @param[in] k  The gain to set.
     */
    void setGain(const double k);
    /*!
     * @brief Gets the gain.
     * @result Gain of transfer function.
     */
    double getGain(void) const;
    /*!
     * @brief Gets the number of poles.
     * @result The number of poles in the transfer function.
     */
    int getNumberOfPoles(void) const;
    /*!
     * @brief Gets the number of zeros.
     * @result The number of zeros in the transfer function.
     */
    int getNumberOfZeros(void) const;
    /*!
     * @{
     * @brief Sets the poles in the transfer function.
     * @param[in] n      The number of poles.
     * @param[in] poles  The poles to set.  This is an array of dimension [n].
     */
    void setPoles(const size_t n, std::complex<double> poles[]);
    /*!
     * @brief Sets the poles in the transfer function.
     * @param[in] poles  The poles to set on the transfer function.
     */
    void setPoles(const std::vector<std::complex<double>> &poles);
    /*! @} */
    /*!
     * @brief Sets the zeros in the transfer function.
     * @param[in] n      The number of zeros.
     * @param[in] zeros  The zeros to set.  This is an array of
     *                   dimension [n].
     */
    void setZeros(const size_t n, std::complex<double> zeros[]);
    /*!
     * @brief Sets the zeros in the transfer function.
     * @param[in] zeros  The zeros to set on the transfer function.
     */
    void setZeros(const std::vector<std::complex<double>> &zeros);
    /*! @} */
    /*!
     * @brief Gets the poles.
     * @result A vector containing the poles.
     */
    std::vector<std::complex<double>> getPoles(void) const;
    /*!
     * @brief Gets the zeros.
     * @result A vector containing the zeros.
     */
    std::vector<std::complex<double>> getZeros(void) const;
    /*!
     * @brief Sets the tolerance in the equality.
     * @param[in] tol   Tolerance.
     */
    void setEqualityTolerance(const double tol = 1.e-12);
 private:
    class ZPKImpl;
    std::unique_ptr<ZPKImpl> pImpl;
}; // End ZPK
}
#endif
