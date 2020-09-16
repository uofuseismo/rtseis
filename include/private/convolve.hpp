#ifndef RTSEIS_PRIVATE_CONVOLVE_HPP
#define RTSEIS_PRIVATE_CONVOLVE_HPP
#include <cmath>
namespace
{

/*!
 * @brief Computes the start and stop indices of the convolution or correlation.
 *        The resulting copy would go from conv(res.first:res.second) where the
 *        upper limit, res.second, is non-inclusive.
 * @param[in] mode   Defines the convolution or correlation mode.
 * @param[in] n1     Length of the first array.
 * @param[in] n2     Length of the second array.
 * @result The start and stop index from which to copy the full convolution.
 * @ingroup rtseis_utils_convolve
 */
std::pair<int,int> computeTrimIndices(
    const RTSeis::Utilities::Math::Convolve::Mode mode,
    const int n1,  const int n2)
{
    int lc;
    int nLeft = 0;
    int nRight = 0;
    // Full
    if (mode == RTSeis::Utilities::Math::Convolve::Mode::FULL)
    {
        lc = n1 + n2 - 1; // Length of full convolution
        nLeft = 0;
        nRight = lc;
    }
    // Valid
    else if (mode == RTSeis::Utilities::Math::Convolve::Mode::VALID)
    {
        lc = std::max(n1, n2) - std::min(n1, n2) + 1;
        if (n2 > n1)
        {
            nLeft = n1/2 + 1;
            if (n1%2 == 0 && n2%2 == 0 && n1 != n2)
            {
                nLeft = std::max(0, nLeft - 1);
            }
            nRight = nLeft + lc;
        }
        else
        {
            nLeft = n2/2 + 1;
            if (n1%2 == 0 && n2%2 == 0 && n1 != n2)
            {
                nLeft = std::max(0, nLeft - 1);
            }
            nRight = nLeft + lc;
        }
    }
    // Same
    else if (mode == RTSeis::Utilities::Math::Convolve::Mode::SAME)
    {
        lc = std::max(n1, n2);
        if (n1 < n2)
        {
            nLeft = n1/2;
            if (n1%2 == 0 && n2%2 == 0){nLeft = std::max(0, nLeft - 1);}
            nRight = nLeft + lc;
        }
        else
        {
            nLeft = n2/2;
            if (n1%2 == 0 && n2%2 == 0){nLeft = std::max(0, nLeft - 1);}
            nRight = nLeft + lc;
        }
    }
    else
    {
        throw std::invalid_argument("Invalid trim type");
    }
#ifndef NDEBUG
    assert(nRight - nLeft == lc);
#endif
    std::pair<int,int> result(nLeft, nRight);
    return result;
}

}
#endif
