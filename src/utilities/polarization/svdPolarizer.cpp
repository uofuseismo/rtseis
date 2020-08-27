#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <array>
#include <ipps.h>
#include <mkl.h>
#define RTSEIS_LOGGING 1
#include "rtseis/enums.h"
#include "rtseis/log.h"
#include "private/throw.hpp"
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
                    const T *__restrict__ A,
                    T *__restrict__ U)
{
    T Uwork[9];
    if (lPositiveP)
    {
       // Rank reduction:
       Uwork[0] = U[0];
       Uwork[1] = U[1];
       Uwork[2] = U[2];

       Uwork[3] = U[3];
       Uwork[4] = U[4];
       Uwork[5] = U[5];

       Uwork[6] = p[0]/magp;
       Uwork[7] = p[1]/magp;
       Uwork[8] = p[2]/magp;
       // First column
       U[0] = Uwork[0]*A[0] + Uwork[3]*A[1] + Uwork[6]*A[2];
       U[1] = Uwork[1]*A[0] + Uwork[4]*A[1] + Uwork[7]*A[2];
       U[2] = Uwork[2]*A[0] + Uwork[5]*A[1] + Uwork[8]*A[2];
       // Second column 
       U[3] = Uwork[0]*A[3] + Uwork[3]*A[4] + Uwork[6]*A[5];
       U[4] = Uwork[1]*A[3] + Uwork[4]*A[4] + Uwork[7]*A[5];
       U[5] = Uwork[2]*A[3] + Uwork[5]*A[4] + Uwork[8]*A[5];
    }
    else
    { 
       std::copy(U, U+6, Uwork);
       // Rank preserving: U_{n} = U_{n-1} A
       // First column
       U[0] = Uwork[0]*A[0] + Uwork[3]*A[1];
       U[1] = Uwork[1]*A[0] + Uwork[4]*A[1];
       U[2] = Uwork[2]*A[0] + Uwork[5]*A[1];
       // Second column 
       U[3] = Uwork[0]*A[3] + Uwork[3]*A[4];
       U[4] = Uwork[1]*A[3] + Uwork[4]*A[4];
       U[5] = Uwork[2]*A[3] + Uwork[5]*A[4];
    }
}

/// Computes the KL transform by computing U_n U_n' d_n
template<typename T> 
void klTransform(const T *__restrict__ U,
                 const T *__restrict__ d,
                 T *dkz, T *dkn, T *dke)
{
    T y0 = U[0]*d[0] + U[1]*d[1] + U[2]*d[2];
    T y1 = U[3]*d[0] + U[4]*d[1] + U[5]*d[2];
    *dkz = U[0]*y0 + U[3]*y1;
    *dkn = U[1]*y0 + U[4]*y1;
    *dke = U[2]*y0 + U[5]*y1;
}

/// Computes the Euclidean norm of a length 3 vector
template<typename T>
inline T norm2(const T v[3])
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

template<typename T>
inline T norm2(const T a, const T b, const T c)
{
    return std::sqrt(a*a + b*b + c*c);
}

/// Computes \f$ U_{n-1}^T d_n \f$
template<typename T>
inline void computeUTd(const T *__restrict__ U,
                       const T *__restrict__ d,
                       T *__restrict__ m)
{
    m[0] = U[0]*d[0] + U[1]*d[1] + U[2]*d[2];
    m[1] = U[3]*d[0] + U[4]*d[1] + U[5]*d[2];
}

/// Computes p = (I - U_{n-1} U_{n-1}^T) d_n = d_n - U_{n-1} 
template<typename T>
inline void computeP(const T *__restrict__ d,
                     const T *__restrict__ U, 
                     const T *__restrict__ m,
                     T *__restrict__ p)
{
    p[0] = d[0] - (U[0]*m[0] + U[3]*m[1]);
    p[1] = d[1] - (U[1]*m[0] + U[4]*m[1]);
    p[2] = d[2] - (U[2]*m[0] + U[5]*m[1]);
}

/// Computes the left singular vector and singular values of Q
inline void svd(double Q[], double S[], double U[])
{
    std::array<double, LWORK> work;
    constexpr lapack_int m = 3;
    constexpr lapack_int n = 3;
    auto info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR,
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
    auto info = LAPACKE_sgesvd_work(LAPACK_COL_MAJOR,
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
    Q[2] = 0;

    Q[3] = 0;
    Q[4] = lambda*S[1];
    Q[5] = 0;

    Q[6] = m[0];
    Q[7] = m[1];
    Q[8] = magp;
}

/// Initializes the recurrence
template<typename T>
inline void initializeUS(const T z, const T n, const T e,
                         T U[], T S[])
{
    // The decomposition is U_0 = [ d_0 / |d_0| ] with singular value |d_0|
    auto magd = norm2(z, n, e);
    U[0] = z/magd;
    U[1] = n/magd;
    U[2] = e/magd;
    U[3] = 0;
    U[4] = 0;
    U[5] = 0;
    S[0] = magd;
    S[1] = 0;
}

/// Do the actual processing
template<class T>
void polarizeSignal(const int npts,
                    const T mLambda, const T mNoise,
                    const RTSeis::ProcessingMode mode,
                    const T *__restrict__ z,
                    const T *__restrict__ n,
                    const T *__restrict__ e,
                    bool *mHaveFirstSample,
                    T *__restrict__ U, T *__restrict__ Q, T *__restrict__ S,
                    T *__restrict__ klz,
                    T *__restrict__ kln,
                    T *__restrict__ kle,
                    T *__restrict__ cosIncidenceAngle,
                    T *__restrict__ rectilinearity)
{
    const T one = 1;
    // Do I want the KL transformed time series?
    bool lwantKL = false;
    if (klz && kln && kle){lwantKL = true;}
    // Rectilinearity?
    bool lwantRect = false;
    if (rectilinearity){lwantRect = true;}
    // Incidence angle?
    bool lwantInc = false;
    if (cosIncidenceAngle){lwantInc = true;}
    // Initialize?
    if (!*mHaveFirstSample)
    {
        initializeUS(z[0], n[0], e[0], U, S);
        *mHaveFirstSample = true;
    }
    // Loop on samples
    T d[3], A[9], m[2], p[3];
    for (int i=0; i<npts; ++i)
    {
        // Extract (z, n, e)
        d[0] = z[i];
        d[1] = n[i];
        d[2] = e[i];
        computeUTd(U, d, m);  // Eqn (5): m = U_{n-1}'d
        computeP(d, U, m, p); // Eqn (6): p = (I - U{n-1} U_{n-1}') d_n
        T magp = norm2(p);    // |p|
        // Instrument noise
        bool lPositiveP = true;
        if (magp < mNoise)
        {
             magp = 0;
             lPositiveP = false;
        }
        // Eqn 11
        fillQ(mLambda, S, m, magp, Q);
        // Diagonalize Q with SVD
        svd(Q, S, A); // A now holds the left singular vectors of Q
        // Update U
        updateU(lPositiveP, p, magp, A, U); 
        // Compute the KL transform.  Basically, we project the data into
        // the space defined by U, i.e., y = U' d, and then we project y
        // back to the original space, i.e., d_{kl} = U U' d.
        if (lwantKL){klTransform(U, d, &klz[i], &kln[i], &kle[i]);}
        // Save rectilinearity
        if (lwantRect)
        {
            rectilinearity[i] = 0;
            if (S[0] > 0)
            {
                rectilinearity[i] = 1 - S[1]/S[0];
            }
        }
        // Save cosine of incidence angle?
        if (lwantInc)
        {
            cosIncidenceAngle[i] = std::min(one, std::abs(U[0]));
        }
    }
    if (mode == RTSeis::ProcessingMode::POST_PROCESSING)
    {
        *mHaveFirstSample = false;
    }
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
    /// The initial conditions
    std::array<double, 3> mInitialConditions = {0.0, 0.0, 0.0};
    /// The forgetting factor.
    double mLambda = 0.9995;
    /// The noise level
    double mNoise = std::numeric_limits<double>::epsilon()*100; //1.e-13;
    /// Flag indicating that the first sample has been set.
    bool mHaveFirstSample = false;
    /// Flag indicating that the initial conditions were set.
    bool mHaveInitialConditions = false;
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
    /// The initial conditions
    std::array<float, 3> mInitialConditions = {0.0f, 0.0f, 0.0f};
    /// The forgetting factor.
    float mLambda = 0.9995;
    /// The noise level
    float mNoise = std::numeric_limits<float>::epsilon()*100; //1.e-5;
    /// Flag indicating that the first sample has been set.
    bool mHaveFirstSample = false;
    /// Flag indicating that the initial conditions were set.
    bool mHaveInitialConditions = false;
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
}

/// Copy c'tor
template<class T>
SVDPolarizer<T>::SVDPolarizer(const SVDPolarizer &polarizer)
{
    *this = polarizer;
}

/// Move c'tor
template<class T>
SVDPolarizer<T>::SVDPolarizer(SVDPolarizer &&polarizer) noexcept
{
    *this = std::move(polarizer);
}

/// Copy assignment operator
template<class T>
SVDPolarizer<T>& SVDPolarizer<T>::operator=(const SVDPolarizer &polarizer)
{
    if (&polarizer == this){return *this;}
    pImpl = std::make_unique<SVDPolarizerImpl> (*polarizer.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
SVDPolarizer<T>& SVDPolarizer<T>::operator=(SVDPolarizer &&polarizer) noexcept
{
    if (&polarizer == this){return *this;}
    pImpl = std::move(polarizer.pImpl);
    return *this;
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
void SVDPolarizer<T>::initialize(const T decayFactor,
                                 const T noise,
                                 const RTSeis::ProcessingMode mode)
{
    clear();
    if (decayFactor <= 0 || decayFactor >= 1)
    {   
        RTSEIS_THROW_IA("decayFactor = %lf must be in range (0,1)",
                        decayFactor); 
    }
    if (noise < 0)
    {
        RTSEIS_THROW_IA("noise = %lf must be at least 0", noise);
    }
    pImpl->mLambda = decayFactor; 
    pImpl->mNoise = noise;
    pImpl->mMode = mode;
    pImpl->mInitialized = true;
}

/// Initialize with default epsilon
template<class T>
void SVDPolarizer<T>::initialize(const T decayFactor,
                                 const RTSeis::ProcessingMode mode)
{
    initialize(decayFactor, std::numeric_limits<T>::epsilon()*100, mode);
}

/// Set the initial conditions
/*
template<class T>
void SVDPolarizer<T>::setInitialConditions(
    const T z, const T n, const T e)
{
    pImpl->mHaveInitialConditions = false;
    pImpl->mHaveFirstSample = false;
    pImpl->mInitialConditions[0] = 0;
    pImpl->mInitialConditions[1] = 0;
    pImpl->mInitialConditions[2] = 0;
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    if (z == 0 && n == 0 && e == 0)
    {
        throw std::invalid_argument("Initial samples cannot all be zero\n"); 
    } 
    pImpl->mInitialConditions[0] = z;
    pImpl->mInitialConditions[1] = n;
    pImpl->mInitialConditions[2] = e;
    initializeUS(pImpl->mInitialConditions[0],
                 pImpl->mInitialConditions[1],
                 pImpl->mInitialConditions[2],
                 pImpl->mU.data(), pImpl->mS.data());
    pImpl->mHaveFirstSample = true;
    pImpl->mHaveInitialConditions = true;
}
*/

/// Reset the initial conditions
template<class T>
void SVDPolarizer<T>::resetInitialConditions()
{
    pImpl->mHaveFirstSample = false;
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    std::fill(pImpl->mQ.begin(), pImpl->mQ.end(), 0);
    std::fill(pImpl->mU.begin(), pImpl->mU.end(), 0);
    std::fill(pImpl->mS.begin(), pImpl->mS.end(), 0);
}

/// Compute everything 
template<class T>
void SVDPolarizer<T>::polarize(const int npts,
                               const T z[], const T n[], const T e[],
                               T *klzIn[], T *klnIn[], T *kleIn[],
                               T *cosIncidenceAngleIn[], T *rectilinearityIn[])
{
    if (npts < 1){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (z == nullptr){RTSEIS_THROW_IA("%s", "z is NULL");}
    if (n == nullptr){RTSEIS_THROW_IA("%s", "n is NULL");}
    if (e == nullptr){RTSEIS_THROW_IA("%s", "e is NULL");}
    // Extract the pointers
    auto klz = *klzIn;
    auto kln = *klnIn;
    auto kle = *kleIn;
    auto cosIncidenceAngle = *cosIncidenceAngleIn;
    auto rectilinearity = *rectilinearityIn;
    if (klz == nullptr){RTSEIS_THROW_IA("%s", "klz is NULL");}
    if (kln == nullptr){RTSEIS_THROW_IA("%s", "kln is NULL");}
    if (kle == nullptr){RTSEIS_THROW_IA("%s", "kle is NULL");}
    if (cosIncidenceAngle == nullptr)
    {
         RTSEIS_THROW_IA("%s", "cosIncidenceAngle is NULL");
    }
    if (rectilinearity == nullptr)
    {
        RTSEIS_THROW_IA("%s", "rectilinearity is NULL");
    }
    // Extract the pointers
    polarizeSignal(npts,
                   pImpl->mLambda, pImpl->mNoise, pImpl->mMode,
                   z, n, e, 
                   &pImpl->mHaveFirstSample,
                   pImpl->mU.data(), pImpl->mQ.data(), pImpl->mS.data(),
                   klz, kln, kle,
                   cosIncidenceAngle, rectilinearity);
    
}

/// Compute KL transform only
template<class T>
void SVDPolarizer<T>::polarize(const int npts, 
                               const T z[], const T n[], const T e[],
                               T *klzIn[], T *klnIn[], T *kleIn[])
{
    if (npts < 1){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (z == nullptr){RTSEIS_THROW_IA("%s", "z is NULL");}
    if (n == nullptr){RTSEIS_THROW_IA("%s", "n is NULL");}
    if (e == nullptr){RTSEIS_THROW_IA("%s", "e is NULL");}
    // Extract the pointers
    auto klz = *klzIn;
    auto kln = *klnIn;
    auto kle = *kleIn;
    if (klz == nullptr){RTSEIS_THROW_IA("%s", "klz is NULL");}
    if (kln == nullptr){RTSEIS_THROW_IA("%s", "kln is NULL");}
    if (kle == nullptr){RTSEIS_THROW_IA("%s", "kle is NULL");}
    T *cosIncidenceAngle = nullptr;
    T *rectilinearity = nullptr;
    // Extract the pointers
    polarizeSignal(npts,
                   pImpl->mLambda, pImpl->mNoise, pImpl->mMode,
                   z, n, e,  
                   &pImpl->mHaveFirstSample,
                   pImpl->mU.data(), pImpl->mQ.data(), pImpl->mS.data(),
                   klz, kln, kle,
                   cosIncidenceAngle, rectilinearity);
}

/// Compute incidence angle and rectilinearity
template<class T>
void SVDPolarizer<T>::polarize(const int npts, 
                               const T z[], const T n[], const T e[],
                               T *cosIncidenceAngleIn[], T *rectilinearityIn[])
{
    if (npts < 1){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (z == nullptr){RTSEIS_THROW_IA("%s", "z is NULL");}
    if (n == nullptr){RTSEIS_THROW_IA("%s", "n is NULL");}
    if (e == nullptr){RTSEIS_THROW_IA("%s", "e is NULL");}
    // Extract the pointers
    auto cosIncidenceAngle = *cosIncidenceAngleIn;
    auto rectilinearity = *rectilinearityIn;
    if (cosIncidenceAngle == nullptr)
    {
         RTSEIS_THROW_IA("%s", "cosIncidenceAngle is NULL");
    }
    if (rectilinearity == nullptr)
    {
        RTSEIS_THROW_IA("%s", "rectilinearity is NULL");
    }
    T *klz = nullptr;
    T *kln = nullptr;
    T *kle = nullptr;
    // Extract the pointers
    polarizeSignal(npts,
                   pImpl->mLambda, pImpl->mNoise, pImpl->mMode,
                   z, n, e,  
                   &pImpl->mHaveFirstSample,
                   pImpl->mU.data(), pImpl->mQ.data(), pImpl->mS.data(),
                   klz, kln, kle,
                   cosIncidenceAngle, rectilinearity);
      
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

/// Modulates a signal using p_n = I_n R_n r_n 
template<typename T>
void RTSeis::Utilities::Polarization::modulateP(
    const int npts,
    const T z[], const T n[], const T e[],
    const T inc[], const T rect[],
    T *pzIn[], T *pnIn[], T *peIn[])
{
    if (npts < 1){return;}
    if (z == nullptr || n == nullptr || e == nullptr ||
        inc == nullptr || rect == nullptr)
    {
        if (z == nullptr){RTSEIS_THROW_IA("%s", "z is NULL");}
        if (n == nullptr){RTSEIS_THROW_IA("%s", "n is NULL");}
        if (e == nullptr){RTSEIS_THROW_IA("%s", "e is NULL");}
        if (inc == nullptr){RTSEIS_THROW_IA("%s", "cosIncidenceAngle is NULL");}
        RTSEIS_THROW_IA("%s", "rectilinearity is NULL");
    }
    auto pz = *pzIn;
    auto pn = *pnIn;
    auto pe = *peIn;
    if (pz == nullptr || pn == nullptr || pe == nullptr)
    {
        if (pz == nullptr){RTSEIS_THROW_IA("%s", "pz is NULL");}
        if (pn == nullptr){RTSEIS_THROW_IA("%s", "pn is NULL");}
        RTSEIS_THROW_IA("%s", "pe is NULL");
    }
    #pragma omp simd
    for (int i=0; i<npts; ++i)
    {
        auto r = inc[i]*rect[i];
        pz[i] = r*z[i];
        pn[i] = r*n[i];
        pe[i] = r*e[i];
    }
}

/// Modulates a signal using s_n = (1 - I_n) R_n r_n
template<typename T>
void RTSeis::Utilities::Polarization::modulateS(
    const int npts,
    const T z[], const T n[], const T e[],
    const T inc[], const T rect[],
    T *szIn[], T *snIn[], T *seIn[])
{
    if (npts < 1){return;}
    if (z == nullptr || n == nullptr || e == nullptr ||
        inc == nullptr || rect == nullptr)
    {
        if (z == nullptr){RTSEIS_THROW_IA("%s", "z is NULL");}
        if (n == nullptr){RTSEIS_THROW_IA("%s", "n is NULL");}
        if (e == nullptr){RTSEIS_THROW_IA("%s", "e is NULL");}
        if (inc == nullptr){RTSEIS_THROW_IA("%s", "cosIncidenceAngle is NULL");}
        RTSEIS_THROW_IA("%s", "rectilinearity is NULL");
    }
    auto sz = *szIn;
    auto sn = *snIn;
    auto se = *seIn;
    if (sz == nullptr || sn == nullptr || se == nullptr)
    {
        if (sz == nullptr){RTSEIS_THROW_IA("%s", "sz is NULL");}
        if (sn == nullptr){RTSEIS_THROW_IA("%s", "sn is NULL");}
        RTSEIS_THROW_IA("%s", "se is NULL");
    }
    #pragma omp simd
    for (int i=0; i<npts; ++i)
    {
        auto r = (1 - inc[i])*rect[i];
        sz[i] = r*z[i];
        sn[i] = r*n[i];
        se[i] = r*e[i];
    }
}

/// Template instantiation
template class RTSeis::Utilities::Polarization::SVDPolarizer<double>;
template class RTSeis::Utilities::Polarization::SVDPolarizer<float>;
template void RTSeis::Utilities::Polarization::modulateP(
    const int npts,
    const double z[], const double n[], const double e[],
    const double inc[], const double rect[],
    double *pzIn[], double *pnIn[], double *peIn[]);
template void RTSeis::Utilities::Polarization::modulateP(
    const int npts,
    const float z[], const float n[], const float e[],
    const float inc[], const float rect[],
    float *pzIn[], float *pnIn[], float *peIn[]);
template void RTSeis::Utilities::Polarization::modulateS(
    const int npts,
    const double z[], const double n[], const double e[],
    const double inc[], const double rect[],
    double *szIn[], double *snIn[], double *seIn[]);
template void RTSeis::Utilities::Polarization::modulateS(
    const int npts,
    const float z[], const float n[], const float e[],
    const float inc[], const float rect[],
    float *szIn[], float *snIn[], float *seIn[]);

namespace
{
template void polarizeSignal(const int npts,
                             const double mLambda, const double mNoise,
                             const RTSeis::ProcessingMode mode,
                             const double *__restrict__ z,
                             const double *__restrict__ n,
                             const double *__restrict__ e,
                             bool *mHaveFirstSample,
                             double *__restrict__ U,
                             double *__restrict__ Q,
                             double *__restrict__ S,
                             double *__restrict__ klz,
                             double *__restrict__ kln,
                             double *__restrict__ kle,
                             double *__restrict__ cosIncidenceAngle,
                             double *__restrict__ rectilinearity);
template void polarizeSignal(const int npts,
                             const float mLambda, const float mNoise,
                             const RTSeis::ProcessingMode mode,
                             const float *__restrict__ z,
                             const float *__restrict__ n,
                             const float *__restrict__ e,
                             bool *mHaveFirstSample,
                             float *__restrict__ U,
                             float *__restrict__ Q,
                             float *__restrict__ S,
                             float *__restrict__ klz,
                             float *__restrict__ kln,
                             float *__restrict__ kle,
                             float *__restrict__ cosIncidenceAngle,
                             float *__restrict__ rectilinearity);
}
