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

using namespace RTSeis::Utilities::Polarization;

namespace
{

/// Computes the Euclidean norm of a length 3 vector
template<typename T>
inline T norm2(const T v[3])
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
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
inline void svd(lapack_int m, lapack_int n, const double q[],
                double s[], double U[])
{
    std::array<double, 9> qwork;
    std::memcpy(qwork.data(), q, 9*sizeof(double));
    lapack_int lwork = 32;
    std::array<double, 32> work;
    auto info = LAPACKE_dgesvd_work(LAPACK_ROW_MAJOR,
                                    'A', 'N', m, n, qwork.data(), 3, 
                                    s, U, 3, NULL, 1, work.data(), lwork); 
#ifdef DEBUG
    assert(info == 0);
#endif
    if (info != 0){RTSEIS_ERRMSG("%s", "dgesvd failed");}
}

/// Computes the left singular vector and singular values of Q
inline void svd(lapack_int m, lapack_int n, const float q[],
                float s[], float U[])
{
    std::array<float, 9> qwork;
    std::memcpy(qwork.data(), q, 9*sizeof(double));
    lapack_int lwork = 32;
    std::array<float, 32> work;
    auto info = LAPACKE_sgesvd_work(LAPACK_ROW_MAJOR,
                                    'A', 'N', m, n, qwork.data(), 3,
                                    s, U, 3, NULL, 1, work.data(), lwork);
#ifdef DEBUG
    assert(info == 0);
#endif
    if (info != 0){RTSEIS_ERRMSG("%s", "dgesvd failed");}
}

}

//============================================================================//
template<>
class SVDPolarizer<double>::SVDPolarizerImpl
{
public:
    /// The 3 x 2 basis for the data polarization at the n'th sample. 
    double mU[6];
    /// The P-wave forgetting factor.
    double mLambdaP = 0.9995;
    /// The S-wave forgetting factor.
    double mLambdaS = 0.9990;
    /// Flag indicating that the first sample has been set.
    bool mHaveFirstSample = false;
    /// Flag indicating taht this is initialized.
    bool mInitialized = false; 
};

//============================================================================//
/// Constructor
template<class T>
SVDPolarizer<T>::SVDPolarizer() :
    pImpl(std::make_unique<SVDPolarizerImpl> ())
{
double U[9], q[9], s[3];
svd(3, 3, q, s, U);
svd(2, 3, q, s, U); 
}

/// Destructor
template<class T>
SVDPolarizer<T>::~SVDPolarizer() = default;

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
template <typename T>
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
