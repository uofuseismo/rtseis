#define RTSEIS_LOGGING 1
#include <cmath>
#include "rtseis/utils/zpk.h"
#include "rtseis/log.h"

/*!
 * @defgroup rtseis_utils_design_iir_zpk ZPK
 * @brief Utility functions for holding a transfer function in 
 *        zero, pole, gain format.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_design
 */

/*!
 * @brief Default constructor. 
 * @ingroup rtseis_utils_design_iir_zpk
 */
ZPK::ZPK(void)
{
    clear();
    return;
}
/*!
 * @brief Constructs a ZPK class from the given zeros, poles, and gain.
 * @param[in] zeros   Zeros to set.
 * @param[in] poles   Poles to set.
 * @param[in] k       Gain to set.
 * @ingroup rtseis_utils_design_iir_zpk 
 */
ZPK::ZPK(const std::vector<std::complex<double>> zeros,
         const std::vector<std::complex<double>> poles,
         const double k)
{
    clear();
    setZeros(zeros);
    setPoles(poles); 
    setGain(k);
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_utils_design_iir_zpk
 */
ZPK::~ZPK(void)
{
    clear();
    return;
}
/*!
 * @brief Prints the ZPK structure.
 * @param[in] fout   File handle to print to.  If fout is NULL then this
 *                   will print to stdout.
 * @ingroup rtseis_utils_design_iir_zpk
 */
void ZPK::print(FILE *fout)
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Gain: %.16lf\n", k_);
    fprintf(f, "Zeros:\n");
    for (size_t i=0; i<z_.size(); i++)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n", std::real(z_[i]), std::imag(z_[i])); 
    }
    fprintf(f, "Poles:\n");
    for (size_t i=0; i<p_.size(); i++)
    {
        fprintf(f, "%+.16lf + %+.16lfi\n", std::real(p_[i]), std::imag(p_[i]));
    }
    return;
}
/*!
 * @brief Clears the structure.
 * @ingroup rtseis_utils_design_iir_zpk
 */
void ZPK::clear(void)
{
    p_.clear();
    z_.clear();
    k_ = 0;
    return;
}
/*!
 * @brief Sets the gain.
 * @ingroup rtseis_utils_design_iir_zpk
 */
void ZPK::setGain(const double k)
{
    k_ = k;
    return;
}
/*!
 * @brief Gets the gain.
 * @retval Gain of transfer function.
 * @ingroup rtseis_utils_design_iir_zpk
 */
double ZPK::getGain(void) const
{
    return k_;
}
/*!
 * @brief Gets the number of poles.
 * @retval Number of poles in transfer function.
 * @ingroup rtseis_utils_design_iir_zpk
 */
int ZPK::getNumberOfPoles(void) const
{
    return static_cast<int> (p_.size());
}
/*!
 * @brief Gets the number of zeros.
 * @retval The number of zeros in the transfer function.
 * @ingroup rtseis_utils_design_iir_zpk
 */
int ZPK::getNumberOfZeros(void) const
{
    return static_cast<int> (p_.size());
}
/*!
 * @brief Sets the poles in the transfer function.
 * @param[in] n      The number of poles.
 * @param[in] poles  The poles to set.  This is an array of dimension [n].
 * @ingroup rtseis_utils_design_iir_zpk
 */
void ZPK::setPoles(const size_t n, std::complex<double> poles[])
{
    if (n > 0 && poles == nullptr)
    {
        RTSEIS_ERRMSG("%s", "Poles is null");
        p_.resize(0);
        return;
    }
    p_.resize(n);
    for (size_t i=0; i<n; i++)
    {
        p_[i] = poles[i];
    }
    return;
}
/*!
 * @brief Sets the poles in the transfer function.
 * @param[in] poles  The poles to set on the transfer function.
 * @ingroup rtseis_utils_design_iir_zpk
 */
void ZPK::setPoles(const std::vector<std::complex<double>> poles)
{
    p_ = poles;
    return;
}
/*!
 * @brief Sets the zeros in the transfer function.
 * @param[in] n      The number of zeros.
 * @param[in] zeros  The zeros to set.  This is an array of dimension [n].
 * @ingroup rtseis_utils_design_iir_zpk
 */
void ZPK::setZeros(const size_t n, std::complex<double> zeros[])
{
    if (n > 0 && zeros == nullptr)
    {
        RTSEIS_ERRMSG("%s", "Zeros is null");
        z_.resize(0);
        return;
    }
    z_.resize(n);
    for (size_t i=0; i<n; i++)
    {
        z_[i] = zeros[i];
    }
    return;
}
/*!
 * @brief Sets the zeros in the transfer function.
 * @param[in] zeros  The zeros to set on the transfer function.
 * @ingroup rtseis_utils_design_iir_zpk
 */
void ZPK::setZeros(const std::vector<std::complex<double>> zeros)
{
    z_ = zeros;
    return;
}
/*!
 * @brief Gets the poles.
 * @ingroup rtseis_utils_design_iir_zpk
 */
std::vector<std::complex<double>> ZPK::getPoles(void) const
{
    return p_;
}
/*!
 * @brief Gets the zeros.
 * @ingroup rtseis_utils_design_iir_zpk
 */
std::vector<std::complex<double>> ZPK::getZeros(void) const
{
    return z_;
}
