#define RTSEIS_LOGGING 1
#include <cmath>
#include <algorithm>
#include "rtseis/utilities/sos.hpp"
#include "rtseis/log.h"

/*!
 * @defgroup rtseis_utils_design_iir_sos SOS
 * @brief Utility functions for holding a transfer function as 
 *        a cascaded series of second order sections.
 * @copyright Ben Baker distributed under the MIT license.
 * @ingroup rtseis_utils_design
 */

/*!
 * @brief Default constructor. 
 * @ingroup rtseis_utils_design_iir_sos
 */
SOS::SOS(void)
{
    clear();
    return;
}
/*!
 * @brief Constructs an SOS class from the given second order sections.
 * @param[in] ns      Number of sections.
 * @param[in] zeros   Zeros to set.
 * @param[in] poles   Poles to set.
 * @param[in] k       Gain to set.
 * @ingroup rtseis_utils_design_iir_sos
 */
SOS::SOS(const int ns,
         const std::vector<double> bs,
         const std::vector<double> as)
{
    clear();
    setSecondOrderSections(ns, bs, as);
    return;
}
/*!
 * @brief Copy operator.
 * @param[in] sos  SOS filter to copy.
 * @ingroup rtseis_utils_design_iir_sos
 */
SOS& SOS::operator=(const SOS &sos)
{
    if (&sos == this){return *this;}
    clear();
    bs_ = sos.bs_;
    as_ = sos.as_;
    ns_ = sos.ns_;
    tol_ = sos.tol_;
    return *this;
}
/*!
 * @brief Equality operator.
 * @param sos  Class to compare to this class.
 * @result True indicates that sos equals this class to within a given
 *         tolerance.
 * @ingroup rtseis_utils_design_iir_sos
 */
bool SOS::operator==(const SOS &sos) const
{
    if (bs_.size() != sos.bs_.size()){return false;}
    if (as_.size() != sos.as_.size()){return false;}
    if (ns_ != sos.ns_){return false;}
    for (size_t i=0; i<bs_.size(); i++)
    {
        if (std::abs(bs_[i] - sos.bs_[i]) > tol_){return false;}
    }
    for (size_t i=0; i<as_.size(); i++)
    {
        if (std::abs(as_[i] - sos.as_[i]) > tol_){return false;}
    }
    return true;
}
/*!
 * @brief Inequality operator.
 * @param[in] sos  Class to compare to this class.
 * @result True indicates that sos does not equal this class to within a
 *         given tolerance.
 * @ingroup rtseis_utils_design_iir_sos 
 */
bool SOS::operator!=(const SOS &sos) const
{
    return !(*this == sos);
}
/*!
 * @brief Copy constructor.
 * @param[in] sos  SOS filter from which to copy.
 * @ingroup rtseis_utils_design_iir_sos
 */
SOS::SOS(const SOS &sos)
{
   *this = sos;
   return;
}
/*!
 * @brief Default destructor.
 * @ingroup rtseis_utils_design_iir_sos
 */
SOS::~SOS(void)
{
    clear();
    return;
}
/*!
 * @brief Clears the class.
 * @ingroup rtseis_utils_design_iir_sos
 */
void SOS::clear(void)
{
    bs_.clear();
    as_.clear();
    ns_ = 0;
    tol_ = defaultTol_;
    return;
}
/*!
 * @brief Prints the SOS structure.
 * @param[in] fout   File handle to print to.  If fout is NULL then this
 *                   will print to stdout.
 * @ingroup rtseis_utils_design_iir_sos
 */
void SOS::print(FILE *fout)
{
    FILE *f = stdout;
    if (fout != nullptr){f = fout;}
    fprintf(f, "Numerator sections\n");
    for (int i=0; i<ns_; i++)
    {
        fprintf(f, "%+.16lf, %+.16lf, %+.16lf\n",
                bs_[3*i], bs_[3*i+1], bs_[3*i+2]);
    }
    fprintf(f, "Denominator sections\n");
    for (int i=0; i<ns_; i++)
    {
        fprintf(f, "%+.16lf, %+.16lf, %+.16lf\n",
                as_[3*i], as_[3*i+1], as_[3*i+2]);
    }
    return;
}
/*!
 * @brief Sets the second order sections on the class.
 * @param[in] ns   The number of second order sections.
 * @param[in] bs   The numerator coefficients.  This has dimension [3 x ns]
 *                 with leading order 3.  Futhermore, the \f$ b[3 i_s] \f$
 *                 cannot be 0 for \f$ i=1,2, \cdots n_s \f$.
 * @param[in] as   The denominator coefficients.  This has dimension [3 x ns]
 *                 with leading order 3.  Futhermore, the \f$ b[3 i_s] \f$
 *                 cannot be 0 for \f$ i=1,2, \cdots n_s \f$.
 * @result 0 indicates success.
 * @ingroup rtseis_utils_design_iir_sos
 */
int SOS::setSecondOrderSections(const int ns,
                                const std::vector<double> bs,
                                const std::vector<double> as)
{
    clear();
    if (ns < 1)
    {
        RTSEIS_ERRMSG("%s", "No sections in SOS filter");
        return -1;
    }
    size_t ns3 = static_cast<size_t> (ns)*3;
    if (ns3 != bs.size())
    {
        RTSEIS_ERRMSG("bs.size() = %ld must equal 3*ns=%ld", bs.size(), ns3);
        return -1;
    }
    if (ns3 != as.size()) 
    {
        RTSEIS_ERRMSG("as.size() = %ld must equal 3*ns=%ld", as.size(), ns3);
        return -1;
    }
    for (int i=0; i<ns; i++)
    {
        if (bs[3*i] == 0)
        {
            RTSEIS_ERRMSG("Leading bs coefficient of section %d is zero", i);
            return -1;
        }
    }
    for (int i=0; i<ns; i++)
    {
        if (as[3*i] == 0)
        {
            RTSEIS_ERRMSG("Leading bs coefficient of section %d is zero", i);
            return -1;
        }
    }
    // It all checks out
    ns_ = ns;
    bs_ = bs;
    as_ = as;
    return 0;
}
/*!
 * @brief Returns the numerator coefficients of the SOS filter.
 * @brief A vector holding the numerator coefficients.  This has dimension
 *        [3 x ns] with leading dimension 3 where ns is given by
 *        SOS::getNumberOfSections().
 * @ingroup rtseis_utils_design_iir_sos
 */
std::vector<double> SOS::getNumeratorCoefficients(void) const
{
    return bs_;
}
/*!
 * @brief Gets the denominator coefficients of the SOS filter.
 * @brief A vector holding the denominator coefficients.  This has dimension
 *        [3 x ns] wiht leading dimension 3 where ns is given by
 *        SOS::getNumberOfSections().
 * @ingroup rtseis_utils_design_iir_sos
 */
std::vector<double> SOS::getDenominatorCoefficients(void) const
{
    return as_;
}
/*!
 * @brief Returns the number of sections.
 * @result The number of second order sections.
 * @ingroup rtseis_utils_design_iir_sos
 */
int SOS::getNumberOfSections(void) const
{
    return ns_;
}
/*!
 * @brief Sets the tolerance in the equality.
 * @param[in] tol   Tolerance.
 * @ingroup rtseis_utils_design_iir_sos
 */
void SOS::setEqualityTolerance(const double tol)
{
    if (tol < 0){RTSEIS_WARNMSG("%s", "Tolerance is negative");}
    tol_ = tol;
    return;
}
