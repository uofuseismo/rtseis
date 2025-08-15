#ifndef RTSEIS_UTILITIES_MATH_POLYNOMIAL_HPP
#define RTSEIS_UTILITIES_MATH_POLYNOMIAL_HPP
#include <cmath>
#include <memory>
#include <complex>
#include <algorithm>
#include "rtseis/vector.hpp"
#include "alignment.hpp"

namespace
{
#ifdef WITH_MKL
#include <mkl_lapacke.h>
#ifndef lapack_int
#define lapack_int MKL_INT
#endif
#ifndef lapack_complex_float
#define lapack_complex_float   MKL_Complex8
#endif
#ifndef lapack_complex_double
#define lapack_complex_double   MKL_Complex16
#endif
template<typename T>
RTSeis::Vector<std::complex<T>> computeEigenvalues(
    RTSeis::Vector<double> &A, const int nIn, const int ldaIn)
{
    auto n = static_cast<lapack_int> (nIn);
    auto lda = static_cast<lapack_int> (ldaIn);
    RTSeis::Vector<std::complex<T>> eigenvalues(n, 0 + 0i);
    std::vector<double> wr(n);
    std::vector<double> wi(n);
    double vl;
    constexpr lapack_int ldvl{1};
    double vr;
    constexpr lapack_int ldvr{1};
    auto info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, A.data(), lda,
                              wr.data(), wi.data(), &vl, ldvl, &vr, ldvr);
    if (info > 0)
    {
        throw std::runtime_error("Failed to compute eigenvalues");
    }
    for (int i = 0; i < n; ++i)
    {
        auto wri = static_cast<T> (wr[i]);
        auto wii = static_cast<T> (wi[i]);
        eigenvalues[i] = std::complex<T> (wri, wii);
    }
    return eigenvalues;
}
#endif
}

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
    auto yPtr = std::assume_aligned<ALIGNMENT> (y.data());
    if (order == 0) // Constant
    {
        const auto p0 = polynomialCoefficients[0];
        std::fill(yPtr, yPtr + nEvaluationPoints, p0);
    }
    else if (order == 1) // Linear
    {
        const auto xPtr
            = std::assume_aligned<ALIGNMENT> (evaluationPoints.data());
        const auto p0 = polynomialCoefficients[0];
        const auto p1 = polynomialCoefficients[1];
        for (int i = 0; i < nEvaluationPoints; ++i)
        {
            yPtr[i] = p0*xPtr[i] + p1;
        }
    }
    else if (order == 2) // Quadratic
    {
        const auto xPtr
            = std::assume_aligned<ALIGNMENT> (evaluationPoints.data());
        const auto p0 = polynomialCoefficients[0];
        const auto p1 = polynomialCoefficients[1];
        const auto p2 = polynomialCoefficients[2];
        for (int i = 0; i < nEvaluationPoints; ++i)
        {
            auto xi = xPtr[i];
            yPtr[i] = p2 + xi*(p1 + xi*p0);
        }
    }
    else // General case (Horner's rule)
    {
        const auto xPtr
            = std::assume_aligned<ALIGNMENT> (evaluationPoints.data());
        const auto p0 = polynomialCoefficients.at(0);
        for (auto i = 0; i < nEvaluationPoints; ++i)
        {
            yPtr[i] = p0*xPtr[i];
        }
        for (auto j = 1; j < order; j++)
        {
            const auto pj = polynomialCoefficients[j];
            for (auto i = 0; i < nEvaluationPoints; ++i)
            {
                yPtr[i] = (pj + yPtr[i])*xPtr[i];
            }
        }
        const auto pn = polynomialCoefficients[order];
        for (auto i = 0; i < nEvaluationPoints; ++i)
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

template<typename T>
RTSeis::Vector<std::complex<T>>
computeRoots(const RTSeis::Vector<T> &coefficients)
{
    auto nCoefficients = static_cast<int> (coefficients.size());
    if (nCoefficients < 1)
    {
        throw std::invalid_argument("No coefficeints");
    }
    auto order = nCoefficients - 1;
    constexpr T zero{0};
    if (coefficients[0] == zero)
    {
        throw std::invalid_argument("Highest order coefficient is zero");
    }
    // Set space for companion matrix
    int n   = order;
    int lda = std::max(8, order);
    RTSeis::Vector<double> A(lda*lda, 0.0);
    double ami = 1.0/coefficients[0]; //coefficient on highest order term
    // Fill out the non-zeros of the companion matrix
    for (int i = 1; i < order + 1; ++i)
    {
        int indx = (i - 1)*lda + 0;
        A[indx] =-ami*coefficients[i];
        // One in subdiagonal 
        if (i < order)
        {
            indx = (i - 1)*lda + i;
            A[indx] = 1.0;
        }
    }
    return ::computeEigenvalues<T>(A, n, lda);
}


}
#endif
