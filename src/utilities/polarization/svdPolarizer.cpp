#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <array>
#include <ipps.h>
#include <mkl.h>
#define RTSEIS_LOGGING 1
#include "rtseis/enums.h"
#include "rtseis/log.h"
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/polarization/svdPolarizer.hpp"

/// Worksapce size
#define LWORK 32

using namespace RTSeis::Utilities::Polarization;

namespace
{

/// Updates U
template<typename T>
inline void updateU(const bool lPositiveP,
                    const T *__restrict__ p, const T magp,
                    const T *__restrict__ A, // LDA = 3
                    T *__restrict__ U)
{
    T Uwork[9];
    if (lPositiveP)
    {
       // Rank reduction:
       Uwork[0] = U[0];
       Uwork[1] = U[1];
       Uwork[2] = p[0]/magp;
 
       Uwork[3] = U[2];
       Uwork[4] = U[3];
       Uwork[5] = p[1]/magp;
 
       Uwork[6] = U[4];
       Uwork[7] = U[5];
       Uwork[8] = p[2]/magp;
       // First row 
       U[0] = Uwork[0]*A[0] + Uwork[1]*A[3] + Uwork[2]*A[6];
       U[1] = Uwork[0]*A[1] + Uwork[1]*A[4] + Uwork[2]*A[7];
       // Second row
       U[2] = Uwork[3]*A[0] + Uwork[4]*A[3] + Uwork[5]*A[6];
       U[3] = Uwork[3]*A[1] + Uwork[4]*A[4] + Uwork[5]*A[7];
       // Third row
       U[0] = Uwork[6]*A[0] + Uwork[7]*A[3] + Uwork[8]*A[6];
       U[1] = Uwork[6]*A[1] + Uwork[7]*A[4] + Uwork[8]*A[7];
    }
    else
    { 
       std::copy(U, U+6, Uwork);
       // Rank preserving: U_{n} = U_{n-1} A
       // First row
       U[0] = Uwork[0]*A[0] + Uwork[1]*A[3];
       U[1] = Uwork[0]*A[1] + Uwork[1]*A[4];
       // Second row
       U[2] = Uwork[2]*A[0] + Uwork[3]*A[3];
       U[3] = Uwork[2]*A[1] + Uwork[3]*A[4];
       // Third row
       U[4] = Uwork[4]*A[0] + Uwork[5]*A[3];
       U[5] = Uwork[4]*A[1] + Uwork[5]*A[4];
    }
}

/// Computes the Euclidean norm of a length 3 vector
template<typename T>
inline T norm2(const T v[3])
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/// Computes \f$ U_{n-1}^T d_n \f$
template<typename T>
inline void computeUTd(const T *__restrict__ U,
                       const T *__restrict__ d,
                       T *__restrict__ m)
{
    m[0] = U[0]*d[0] + U[2]*d[1] + U[4]*d[2];
    m[1] = U[0]*d[0] + U[3]*d[1] + U[5]*d[2];
}

/// Computes p = (I - U_{n-1} U_{n-1}^T) d_n = d_n - U_{n-1} 
template<typename T>
inline void computeP(const T *__restrict__ d,
                     const T *__restrict__ U, 
                     const T *__restrict__ m,
                     T *__restrict__ p)
{
    p[0] = d[0] - (U[0]*m[0] + U[1]*m[1]);
    p[1] = d[1] - (U[2]*m[0] + U[3]*m[1]);
    p[2] = d[2] - (U[4]*m[0] + U[5]*m[1]);
}

/// Scales a length 3 vector y = alpha*x
template<typename T>
inline void scal3(const T alpha, const T x[3], T y[3])
{
    y[0] = alpha*x[0];
    y[1] = alpha*x[1];
    y[2] = alpha*x[2];
}

/// Computes the left singular vector and singular values of Q
inline void svd(double Q[], double S[], double U[])
{
    std::array<double, LWORK> work;
    constexpr lapack_int m = 3;
    constexpr lapack_int n = 3;
    auto info = LAPACKE_dgesvd_work(LAPACK_ROW_MAJOR,
                                    'A', 'N', m, n, Q, 3, 
                                    S, U, 3, NULL, 1, work.data(), LWORK);
#ifdef DEBUG
    assert(info == 0);
#endif
    if (info != 0){RTSEIS_ERRMSG("%s", "dgesvd failed");}
}

inline void svd(float Q[], float S[], float U[])
{
    std::array<float, LWORK> work;
    constexpr lapack_int m = 3;
    constexpr lapack_int n = 3;
    auto info = LAPACKE_sgesvd_work(LAPACK_ROW_MAJOR,
                                    'A', 'N', m, n, Q, 3,
                                    S, U, 3, NULL, 1, work.data(), LWORK);
#ifdef DEBUG
    assert(info == 0);
#endif
    if (info != 0){RTSEIS_ERRMSG("%s", "dgesvd failed");}
}

/// Fills the Q matrix ala Eqn 11
template<typename T>
inline void fillQ(const T lambda,
                  const T *__restrict__ S,
                  const T *__restrict__ m,
                  const T magp, T *__restrict__ Q)
{
    Q[0] = lambda*S[0];
    Q[1] = 0;
    Q[2] = m[0];

    Q[3] = 0;
    Q[4] = lambda*S[1];
    Q[5] = m[1];

    Q[6] = 0;
    Q[7] = 0;
    Q[8] = magp;
}

}

//============================================================================//
template<>
class SVDPolarizer<double>::SVDPolarizerImpl
{
public:
    /// Space for the 3 x 3 Q matrix
    std::array<double, 9> mQ;
    /// The 3 x 2 basis for the data polarization at the n'th sample. 
    std::array<double, 6> mU;
    /// The diagonal matrix of singular values
    std::array<double, 3> mS;
    /// The forgetting factor.
    double mLambda = 0.9995;
    /// The noise level
    double mNoise = 1.e-13;
    /// Flag indicating that the first sample has been set.
    bool mHaveFirstSample = false;
    /// Flag indicating that the class is initialized.
    bool mInitialized = false;
    /// The processing mode
    RTSeis::ProcessingMode mMode = RTSeis::ProcessingMode::POST_PROCESSING;
};

template<>
class SVDPolarizer<float>::SVDPolarizerImpl
{
public:
    /// Space for the 3 x 3 Q matrix
    std::array<float, 9> mQ;
    /// The 3 x 2 basis for the data polarization at the n'th sample.
    std::array<float, 6> mU;
    /// The diagonal matrix of singular values
    std::array<float, 3> mS;
    /// The forgetting factor.
    float mLambda = 0.9995;
    /// The noise level
    float mNoise = 1.e-5;
    /// Flag indicating that the first sample has been set.
    bool mHaveFirstSample = false;
    /// Flag indicating that the class is initialized.
    bool mInitialized = false;
    /// The processing mode
    RTSeis::ProcessingMode mMode = RTSeis::ProcessingMode::POST_PROCESSING;
};

//============================================================================//
/// Constructor
template<class T>
SVDPolarizer<T>::SVDPolarizer() :
    pImpl(std::make_unique<SVDPolarizerImpl> ())
{
/*
double U[9], q[9], s[3];
svd(3, 3, q, s, U);
svd(2, 3, q, s, U); 
*/
}

/// Destructor
template<class T>
SVDPolarizer<T>::~SVDPolarizer() = default;

/// Destroys the class
template<class T>
void SVDPolarizer<T>::clear() noexcept
{
    std::fill(pImpl->mU.begin(), pImpl->mU.end(), 0);
    std::fill(pImpl->mQ.begin(), pImpl->mQ.end(), 0);
    std::fill(pImpl->mS.begin(), pImpl->mS.end(), 0);
    pImpl->mLambda = 0.9995;
    pImpl->mHaveFirstSample = false;
    pImpl->mInitialized = false;
    pImpl->mMode = RTSeis::ProcessingMode::POST_PROCESSING;
}

/// Initialize
template<class T>
void SVDPolarizer<T>::initialize(const double decayFactor,
                                 const RTSeis::ProcessingMode mode)
{
    clear();
    if (decayFactor <= 0 || decayFactor >= 1)
    {
        RTSEIS_THROW_IA("decayFactor = %lf must be in range (0,1)",
                        decayFactor); 
    }
    pImpl->mLambda = decayFactor; 
    pImpl->mMode = mode;
    pImpl->mInitialized = true;
}

/// 
template<class T>
void SVDPolarizer<T>::polarize(const int npts, const T z[], const T n[], const T e[],
                               T incidenceAngle[], T rectilinearity[])
{
    if (npts < 1){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (z == nullptr){RTSEIS_THROW_RTE("%s", "z is NULL");}
    if (n == nullptr){RTSEIS_THROW_RTE("%s", "n is NULL");}
    if (e == nullptr){RTSEIS_THROW_RTE("%s", "e is NULL");}
    T *U = pImpl->mU.data();
    T *Q = pImpl->mQ.data();
    T *S = pImpl->mS.data();
    // Initialize?
    T d[3];
    int ibeg = 0;
    if (!pImpl->mHaveFirstSample)
    {
        // Extract (z, n, e)
        d[0] = z[0];
        d[1] = n[0];
        d[2] = e[0];
        // Eqn 7: X_0 = [ d_0 ] 
        // The decomposition is U_0 = [ d_0 / |d_0| ] with singular value |d_0|
        std::fill(pImpl->mU.begin(), pImpl->mU.end(), 0);
        auto magd = norm2(d);
        pImpl->mU[0] = d[0]/magd;
        pImpl->mU[2] = d[1]/magd;
        pImpl->mU[4] = d[2]/magd;
        std::fill(pImpl->mS.begin(), pImpl->mS.end(), 0);
        pImpl->mS[0] = magd;
        ibeg = 1;
        pImpl->mHaveFirstSample = true;
    }
    // Loop on samples
    T A[9], m[2], p[3];
    for (int i=ibeg; i<npts; ++i)
    {
        // Extract (z, n, e)
        d[0] = z[i];
        d[1] = n[i];
        d[2] = e[i]; 
        computeUTd(U, d, m);  // Eqn (5)
        computeP(d, U, m, p); // Eqn (6)
        T magp = norm2(p);    // |p| 
        // Instrument noise
        bool lPositiveP = true;
        if (magp < pImpl->mNoise)
        {
             magp = 0;
             lPositiveP = false;
        }
        // Eqn 11
        fillQ(pImpl->mLambda, S, m, magp, Q);
        // Diagonalize Q with SVD
        std::copy(Q, Q+9, A);
        svd(Q, S, A); // A now holds the left singular vectors of Q
        // Update U
        updateU(lPositiveP, p, magp, A, U);
    }
}

/// Initialized?
template<class T>
bool SVDPolarizer<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

///===========================================================================//
//                                 Utility Functions                          //
//============================================================================//

/// Polarizes a signal using p_n = I_n R_n r_n 
template<typename T>
void RTSeis::Utilities::Polarization::pPolarize(
    const int n, const T zne[],
    const T inc[], const T rect[], const T r[],
    T pPolarized[])
{
    #pragma omp simd
    for (int i=0; i<n; ++i)
    {
        auto ir = inc[i]*rect[i];
        pPolarized[3*i]   = ir;
        pPolarized[3*i+1] = ir;
        pPolarized[3*i+2] = ir;
    }
    #pragma omp simd
    for (int i=0; i<3*n; ++i){pPolarized[i] = pPolarized[i]*r[i];}
}
/// Polarizes a signal using s_n = (1 - I_n) R_n r_n
template<typename T>
void RTSeis::Utilities::Polarization::sPolarize(
    const int n, const T zne[],
    const T inc[], const T rect[], const T r[],
    T pPolarized[])
{
    #pragma omp simd
    for (int i=0; i<n; ++i)
    {
        auto ir = 1 - inc[i]*rect[i];
        pPolarized[3*i]   = ir;
        pPolarized[3*i+1] = ir;
        pPolarized[3*i+2] = ir;
    }
    #pragma omp simd
    for (int i=0; i<3*n; ++i){pPolarized[i] = pPolarized[i]*r[i];}
}

/// Template instantiation
template class RTSeis::Utilities::Polarization::SVDPolarizer<double>;
template class RTSeis::Utilities::Polarization::SVDPolarizer<float>;
