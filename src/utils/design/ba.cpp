#define RTSEIS_LOGGING 1
#include <cmath>
#include "rtseis/utils/ba.h"
#include "rtseis/log.h"

/*!
 * @defgroup rtseis_utils_design_iir_ba BA
 * @brief Utility functions for holding a transfer function in 
 *        terms of numerator and denominator coefficients.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_design
 */

/*!
 * @brief Default constructor. 
 * @ingroup rtseis_utils_design_iir_ba
 */
BA::BA(void)
{
    clear();
    return;
}
/*!
 * @brief Constructs a transfer function class from the given numerator and 
*         denminator coefficients.
 * @param[in] b     The numerator coefficients.
 * @param[in] a     The denomiantor coefficients.
 * @ingroup rtseis_utils_design_iir_ba
 */
BA::BA(const std::vector<double> b, const std::vector<double> a)
{
    clear();
    setNumeratorCoefficients(b);
    setDenominatorCoefficients(a);
    return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_utils_design_iir_ba
 */
BA::~BA(void)
{
    clear();
    return;
}
/*!
 * @brief Prints the BA structure.
 * @param[in] fout   File handle to print to.  If fout is NULL then this
 *                   will print to stdout.
 * @ingroup rtseis_utils_design_iir_ba
 */
void BA::print(FILE *fout)
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Numerator Coefficients:\n");
    for (size_t i=0; i<b_.size(); i++)
    {
        fprintf(f, "%+.16lf\n", b_[i]);
    }
    fprintf(f, "Denominator Coefficients:\n");
    for (size_t i=0; i<a_.size(); i++)
    {
        fprintf(f, "%+.16lf\n", a_[i]);
    }
    return;
}
/*!
 * @brief Clears the structure.
 * @ingroup rtseis_utils_design_iir_ba
 */
void BA::clear(void)
{
    b_.clear();
    a_.clear();
    return;
}
/*!
 * @brief Gets the number of numerator coefficients.
 * @retval Number of numerator coefficients in transfer function.
 * @ingroup rtseis_utils_design_iir_ba
 */
int BA::getNumberOfNumeratorCoefficients(void) const
{
    return static_cast<int> (b_.size());
}
/*!
 * @brief Gets the number of denominator coefficients.
 * @retval The number of denominator coefficients in the transfer function.
 * @ingroup rtseis_utils_design_iir_ba
 */
int BA::getNumberOfDenominatorCoefficients(void) const
{
    return static_cast<int> (a_.size());
}
/*!
 * @brief Sets the numerator coefficients of the transfer function.
 * @param[in] n    The number of numerator coefficients.
 * @param[in] b    The numerator coefficients to set.  This is an array of
 *                 dimension [n].
 * @ingroup rtseis_utils_design_iir_ba
 */
void BA::setNumeratorCoefficients(const size_t n, double b[])
{
    if (n > 0 && b == nullptr)
    {
        RTSEIS_ERRMSG("%s", "b is null");
        b_.resize(0);
        return;
    }
    b_.resize(n);
    for (size_t i=0; i<n; i++)
    {
        b_[i] = b[i];
    }
    return;
}
/*!
 * @brief Sets the numerator coefficients of the transfer function.
 * @param[in] b   Numerator coefficients to set on the transfer function.
 * @ingroup rtseis_utils_design_iir_ba
 */
void BA::setNumeratorCoefficients(const std::vector<double> b)
{
    b_ = b;
    return;
}
/*!
 * @brief Sets the denominator coefficients of the transfer function.
 * @param[in] n    The number of numerator coefficients.
 * @param[in] a    The denominator coefficients to set.  This is an array of
 *                 dimension [n].
 * @ingroup rtseis_utils_design_iir_ba
 */
void BA::setDenominatorCoefficients(const size_t n, double a[])
{
    if (n > 0 && a == nullptr)
    {
        RTSEIS_ERRMSG("%s", "a is null");
        a_.resize(0);
        return;
    }
    if (n > 0 && a[0] == 0){RTSEIS_WARNMSG("%s", "a[0] = 0");}
    a_.resize(n);
    for (size_t i=0; i<n; i++)
    {
        a_[i] = a[i];
    }
    return;
}
/*! 
 * @brief Sets the numerator coefficients of the transfer function.
 * @param[in] b   Numerator coefficients to set on the transfer function.
 * @ingroup rtseis_utils_design_iir_ba
 */
void BA::setDenominatorCoefficients(const std::vector<double> a)
{
    if (a.size() > 0)
    {
        if (a[0] == 0){RTSEIS_WARNMSG("%s", "a[0] = 0");}
    }
    a_ = a;
    return;
}
/*!
 * @brief Gets the numerator coefficients.
 * @ingroup rtseis_utils_design_iir_ba
 */
std::vector<double> BA::getNumeratorCoefficients(void) const
{
    return b_;
}
/*!
 * @brief Gets the denominator coefficients.
 * @ingroup rtseis_utils_design_iir_ba
 */
std::vector<double> BA::getDenominatorCoefficients(void) const
{
    return a_;
}
