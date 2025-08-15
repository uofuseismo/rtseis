#ifndef RTSEIS_UTILITIES_MATH_POLYNOMIAL_HPP
#define RTSEIS_UTILITIES_MATH_POLYNOMIAL_HPP
#include <cmath>
#include <complex>
#include <algorithm>
#include "rtseis/vector.hpp"
namespace RTSeis::Utilitities::Math::Polynomial
{

/// @brief Evalutes the polynomial:
///           p(x) = p_0 x^n + p_1 x^{n-1} + \cdots + p_n  
///        at the given evaluation points.
/// @param[in] polynomialCoefficients  The polynomial coefficients
///                                    ordered as shown above.
/// @param[in] evaluationPoints The x_i's at which to evaluate the polynomial.
/// @result The polynomial evaluated at each evaluation point.
template<typename U>
RTSeis::Vector<U> 
evaluate(
    const RTSeis::Vector<U> &polynomialCoefficients,
    const RTSeis::Vector<U> &evaluationPoints) 
{
    if (polynomialCoefficients.empty())
    {
        throw std::invalid_argument(
           "No coefficients in polynomialCoeffiecients");
    }
    auto order = static_cast<int> (polynomialCoefficients.size()) - 1;
    auto nEvaluationPoints = static_cast<int> (evaluationPoints.size());
    RTSeis::Vector<U> y;
    if (nEvaluationPoints < 1){return y;}
    y.resize(nEvaluationPoints, 0);
    // Expand the constant case 
    if (order == 0) // Constant
    {
        std::fill(y.begin(), y.end(), polynomialCoefficients[0]);
    }
    else if (order == 1) // Linear
    {
        const auto xPtr = evaluationPoints.data();
        auto yPtr = y.data();
        const auto p0 = polynomialCoefficients[0];
        const auto p1 = polynomialCoefficients[1];
        for (int i = 0; i < nEvaluationPoints; i++)
        {
            yPtr[i] = p0*xPtr[i] + p1;
        }
    }
    else if (order == 2) // Quadratic
    {
        const auto xPtr = evaluationPoints.data();
        auto yPtr = y.data();
        const auto p0 = polynomialCoefficients[0];
        const auto p1 = polynomialCoefficients[1];
        const auto p2 = polynomialCoefficients[2];
        for (int i = 0; i < nEvaluationPoints; i++)
        {
            auto xi = xPtr[i];
            yPtr[i] = p2 + xi*(p1 + xi*p0);
        }
    }
    else // General case (Horner's rule)
    {
        const auto xPtr = evaluationPoints.data();
        auto yPtr = y.data();
        const auto p0 = polynomialCoefficients.at(0);
        for (auto i = 0; i < nEvaluationPoints; i++)
        {
            yPtr[i] = p0*xPtr[i];
        }
        for (auto j = 1; j < order; j++)
        {
            const auto pj = polynomialCoefficients[j];
            for (auto i = 0; i < nEvaluationPoints; i++)
            {
                yPtr[i] = (pj + yPtr[i])*xPtr[i];
            }
        }
        const auto pn = polynomialCoefficients[order];
        for (auto i = 0; i < nEvaluationPoints; i++)
        {
            yPtr[i] = pn + yPtr[i];
        }
    }
    return y;
}

/*
/// @brief Evalutes the polynomial:
///           p(x) = p_0 x^n + p_1 x^{n-1} + \cdots + p_n  
///        at the given evaluation points.
/// @param[in] polynomialCoefficients  The polynomial coefficients
///                                    ordered as shown above.
/// @param[in] evaluationPoints The x_i's at which to evaluate the polynomial.
/// @result The polynomial evaluated at each evaluation point.
template<typename U>
RTSeis::Vector<std::complex<U>> 
evaluate(
    const RTSeis::Vector<std::complex<U>> &polynomialCoefficients,
    const RTSeis::Vector<std::complex<U>> &evaluationPoints) 
{
    if (polynomialCoefficients.empty())
    {
        throw std::invalid_argument(
           "No coefficients in polynomialCoeffiecients");
    }
    auto order = static_cast<int> (polynomialCoefficients.size()) - 1;
    auto nEvaluationPoints = static_cast<int> (evaluationPoints.size());
    RTSeis::Vector<std::complex<U>> y;
    if (nEvaluationPoints < 1){return y;}
    y.resize(nEvaluationPoints, 0);
    // Expand the constant case 
    if (order == 0) // Constant
    {
        std::fill(y.begin(), y.end(), polynomialCoefficients[0]);
    }
    else if (order == 1) // Linear
    {
        const auto xPtr = evaluationPoints.data();
        auto yPtr = y.data();
        const auto p0 = polynomialCoefficients[0];
        const auto p1 = polynomialCoefficients[1];
        for (int i = 0; i < nEvaluationPoints; i++)
        {
            yPtr[i] = p0*xPtr[i] + p1;
        }
    }
    else if (order == 2) // Quadratic
    {
        const auto xPtr = evaluationPoints.data();
        auto yPtr = y.data();
        const auto p0 = polynomialCoefficients[0];
        const auto p1 = polynomialCoefficients[1];
        const auto p2 = polynomialCoefficients[2];
        for (int i = 0; i < nEvaluationPoints; i++)
        {
            auto xi = xPtr[i];
            yPtr[i] = p2 + xi*(p1 + xi*p0);
        }
    }
    else // General case (Horner's rule)
    {
        const auto xPtr = evaluationPoints.data();
        auto yPtr = y.data();
        const auto p0 = polynomialCoefficients.at(0);
        for (auto i = 0; i < nEvaluationPoints; i++)
        {
            yPtr[i] = p0*xPtr[i];
        }
        for (auto j = 1; j < order; j++)
        {
            const auto pj = polynomialCoefficients[j];
            for (auto i = 0; i < nEvaluationPoints; i++)
            {
                yPtr[i] = (pj + yPtr[i])*xPtr[i];
            }
        }
        const auto pn = polynomialCoefficients[order];
        for (auto i = 0; i < nEvaluationPoints; i++)
        {
            yPtr[i] = pn + yPtr[i];
        }
    }
    return y;
}
*/

}
#endif
