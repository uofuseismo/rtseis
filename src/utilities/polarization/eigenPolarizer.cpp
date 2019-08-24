#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <array>
#include <vector>
#include <mkl_lapacke.h>
#include <ipps.h>
#include "rtseis/enums.h"
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/polarization/eigenPolarizer.hpp"
#include "rtseis/utilities/rotate/utilities.hpp"

using namespace RTSeis::Utilities::Polarization;

namespace
{

/// Computes cross-variance matrix (Jurkevics, 1988)
void computeCrossVarianceMatrix(const int npts,
                                const double z[],
                                const double n[],
                                const double e[],
                                double crossVarianceMatrix[9])
{
    // Compute the means
    double zMean, nMean, eMean;
    ippsMean_64f(z, npts, &zMean);
    ippsMean_64f(n, npts, &nMean);
    ippsMean_64f(e, npts, &eMean);
    // Compute the cross power covariance
    double zz = 0;
    double zn = 0;
    double ze = 0;
    double nn = 0;
    double ne = 0;
    double ee = 0;
    #pragma omp simd reduction(+:zz, zn, ze, nn, ne, ee)
    for (auto i=0; i<npts; ++i)
    {
        // Strip out DC shift 
        auto zDemean = (z[i] - zMean);
        auto nDemean = (n[i] - nMean);
        auto eDemean = (e[i] - eMean);
        // row 1
        zz = zz + zDemean*zDemean;
        zn = zn + zDemean*nDemean;
        ze = ze + zDemean*eDemean;
        // row 2 
        nn = nn + nDemean*nDemean;
        ne = ne + nDemean*eDemean;
        // row 3
        ee = ee + eDemean*eDemean;
    }
    // Create the cross variance matrix (column major)
    crossVarianceMatrix[0] = zz/static_cast<double> (npts);
    crossVarianceMatrix[1] = zn/static_cast<double> (npts);
    crossVarianceMatrix[2] = ze/static_cast<double> (npts);
    crossVarianceMatrix[3] = zn/static_cast<double> (npts);
    crossVarianceMatrix[4] = nn/static_cast<double> (npts);
    crossVarianceMatrix[5] = ne/static_cast<double> (npts);
    crossVarianceMatrix[6] = ze/static_cast<double> (npts);
    crossVarianceMatrix[7] = ne/static_cast<double> (npts);
    crossVarianceMatrix[8] = ee/static_cast<double> (npts);
}

/// Computes cross-variance matrix (Jurkevics, 1988)
void computeCrossVarianceMatrix(const int npts,
                                const float z[],
                                const float n[],
                                const float e[],
                                double crossVarianceMatrix[9])
{
    // Compute the means
    float zMean, nMean, eMean;
    ippsMean_32f(z, npts, &zMean, ippAlgHintAccurate);
    ippsMean_32f(n, npts, &nMean, ippAlgHintAccurate);
    ippsMean_32f(e, npts, &eMean, ippAlgHintAccurate);
    // Compute the cross power covariance
    double zz = 0;
    double zn = 0;
    double ze = 0;
    double nn = 0;
    double ne = 0;
    double ee = 0;
    #pragma omp simd reduction(+:zz, zn, ze, nn, ne, ee)
    for (auto i=0; i<npts; ++i)
    {
        // Strip out DC shift
        auto zDemean = static_cast<double> (z[i] - zMean);
        auto nDemean = static_cast<double> (n[i] - nMean);
        auto eDemean = static_cast<double> (e[i] - eMean);
        // row 1
        zz = zz + zDemean*zDemean;
        zn = zn + zDemean*nDemean;
        ze = ze + zDemean*eDemean;
        // row 2
        nn = nn + nDemean*nDemean;
        ne = ne + nDemean*eDemean;
        // row 3
        ee = ee + eDemean*eDemean;
    }
    // Create the cross variance matrix (column major)
    crossVarianceMatrix[0] = zz/static_cast<double> (npts);
    crossVarianceMatrix[1] = zn/static_cast<double> (npts);
    crossVarianceMatrix[2] = ze/static_cast<double> (npts);
    crossVarianceMatrix[3] = zn/static_cast<double> (npts);
    crossVarianceMatrix[4] = nn/static_cast<double> (npts);
    crossVarianceMatrix[5] = ne/static_cast<double> (npts);
    crossVarianceMatrix[6] = ze/static_cast<double> (npts);
    crossVarianceMatrix[7] = ne/static_cast<double> (npts);
    crossVarianceMatrix[8] = ee/static_cast<double> (npts);
}

/// Computes eigendecomposition of 3 x 3 symmetric matrix
void computeEigendecomposition(const double A[9], double eig[3], double X[9])
{
    constexpr int n = 3;
    constexpr int lwork = 32;
    constexpr int lda = 3;
    std::array<double, 32> work;
    std::array<double, 9> Atemp;
    std::copy(A, A+9, Atemp.data());
    // Compute eigenvalues - convention is ascending order
#ifdef DEBUG
    auto info =
#endif
    LAPACKE_dsyev_work(LAPACK_COL_MAJOR, // Matrix of eigenvectors is col major
                       'V',              // Eigenvalues and eigenvectors
                       'U',              // Upper-triangle
                       n,                // A is [n x n]  = [3 x 3]
                       Atemp.data(),     // Will hold the eigenvectors
                       lda,              // Leading dimension is 3
                       eig,              // Eigenvalues
                       work.data(),      // Workspace
                       lwork);           // Size of workspace
    
#ifdef DEBUG
    assert(info == 0);
#endif
    // Reorder the eigenvalues so that they are in descending order
    std::reverse(eig, eig+3);
    std::copy(Atemp.begin()+6, Atemp.begin()+9, X);
    std::copy(Atemp.begin()+3, Atemp.begin()+6, X+3);
    std::copy(Atemp.begin(),   Atemp.begin()+3, X+6);
}

}

template<class T>
class EigenPolarizer<T>::EigenPolarizerImpl
{
public:
    /// Default constructor
    EigenPolarizerImpl() = default;
    /// Copy constructor
    EigenPolarizerImpl(const EigenPolarizerImpl &eigen)
    {
        *this = eigen;
    }
    /// Destructor
    ~EigenPolarizerImpl()
    {
        clear();
    }
    /// Copy assignment
    EigenPolarizerImpl& operator=(const EigenPolarizerImpl &eigen)
    {
        if (&eigen == this){return *this;}
        if (!eigen.mInitialized){return *this;}
        // Initialize class (sets nsamples + precision and allocates memory)
        initialize(eigen.mSamples, eigen.mPrecision);
        // Copy the basics 
        mCrossVarianceMatrix = eigen.mCrossVarianceMatrix;
        mEigenvalues = eigen.mEigenvalues;
        mEigenvectors = eigen.mEigenvectors;
        mIncidenceAngle = eigen.mIncidenceAngle;
        mRectilinearity = eigen.mRectilinearity;
        mAzimuth = eigen.mAzimuth;
        mBackAzimuth = eigen.mBackAzimuth;
        mHaveSignals = eigen.mHaveSignals;
        mHaveEigenvalues = eigen.mHaveEigenvalues;
        mInitialized = eigen.mInitialized;
        // Copy the arrays
        if (mPrecision == RTSeis::Precision::DOUBLE)
        {
            ippsCopy_64f(eigen.mVertical64f, mVertical64f, mSamples);
            ippsCopy_64f(eigen.mNorth64f, mNorth64f, mSamples);
            ippsCopy_64f(eigen.mEast64f, mEast64f, mSamples);
        }
        else
        {
            ippsCopy_32f(eigen.mVertical32f, mVertical32f, mSamples);
            ippsCopy_32f(eigen.mNorth32f, mNorth32f, mSamples);
            ippsCopy_32f(eigen.mEast32f, mEast32f, mSamples);
        }
        return *this;
    }
    /// Clears the class
    void clear()
    {
        if (mVertical64f){ippsFree(mVertical64f);}
        if (mNorth64f){ippsFree(mNorth64f);}
        if (mEast64f){ippsFree(mEast64f);}
        if (mVertical32f){ippsFree(mVertical32f);}
        if (mNorth32f){ippsFree(mNorth32f);}
        if (mEast32f){ippsFree(mEast32f);}
        mVertical64f = nullptr;
        mNorth64f = nullptr;
        mEast64f = nullptr;
        mVertical32f = nullptr;
        mNorth32f = nullptr;
        mEast32f = nullptr;
        mSamples = 0;
        mHaveSignals = false;
        mHaveEigenvalues = false;
        mInitialized = false;
    }
    /// Initialize the class
    void initialize(const int nSamples,
                    const RTSeis::Precision precision)
    {
        mSamples = nSamples;
        if (precision == RTSeis::Precision::DOUBLE)
        {
            mVertical64f = ippsMalloc_64f(mSamples);
            mNorth64f = ippsMalloc_64f(mSamples);
            mEast64f = ippsMalloc_64f(mSamples);
        }
        else
        {
            mVertical32f = ippsMalloc_32f(mSamples);
            mNorth32f = ippsMalloc_32f(mSamples);
            mEast32f = ippsMalloc_32f(mSamples);
        }
        mPrecision = precision;
        mInitialized = true;
    }
    /// Set the input signals and eigendecomposes crossvariance matrix
    void setSignalsAndDecompose(const int nSamples,
                                const double vertical[],
                                const double north[],
                                const double east[])
    {
        // Set the signals
        ippsCopy_64f(vertical, mVertical64f, nSamples);
        ippsCopy_64f(north, mNorth64f, nSamples);
        ippsCopy_64f(east, mEast64f, nSamples);
        // Compute the cross-variance matrix
        computeCrossVarianceMatrix(nSamples,
                                   mVertical64f,
                                   mNorth64f,
                                   mEast64f,
                                   mCrossVarianceMatrix.data());
        // Get the eigendecomposition and related metrics
        computeMetrics();
        mHaveSignals = true;
    }
    /// Set the input signals and eigendecomposes crossvariance matrix
    void setSignalsAndDecompose(const int nSamples,
                                const float vertical[],
                                const float north[],
                                const float east[])
    {
        // Set the signals
        ippsCopy_32f(vertical, mVertical32f, nSamples);
        ippsCopy_32f(north, mNorth32f, nSamples);
        ippsCopy_32f(east, mEast32f, nSamples);
        // Compute the cross-variance matrix
        computeCrossVarianceMatrix(nSamples,
                                   mVertical32f,
                                   mNorth32f,
                                   mEast32f,
                                   mCrossVarianceMatrix.data());
         // Get the eigendecomposition and related metrics
         computeMetrics();
         mHaveSignals = true;
    }
    /// Does the work - computes the eigendecomposition and computes the derived
    /// metrics.
    void computeMetrics()
    {
        // Get the eigendecomposition
        computeEigendecomposition(mCrossVarianceMatrix.data(),
                                  mEigenvalues.data(),
                                  mEigenvectors.data());
        // The frame is (Z,N,E).
        auto u11 = mEigenvectors[0]; // Z
        auto u21 = mEigenvectors[1]; // N
        auto u31 = mEigenvectors[2]; // E
        // Force ray to come up out of ground
        if (u11 < 0)
        {
            u11 =-u11;
            u21 =-u21;
            u31 =-u31;
        }
        mIncidenceAngle = std::acos(u11); // Absolute value already computed
        mRectilinearity = 1 - (mEigenvalues[1] + mEigenvalues[2])
                             /(2*mEigenvalues[0]);
        // atan2(N,E) = atan2(y,x).  This is measured positive counterclockwise
        // from x.  We measure counterclockwise (minus sign) from north which is
        // +90 degrees.   Note that atan2 is in range [-pi,pi).  Thus, if
        // mAzimuth is negative we enforce the convention it be between [0,360]
        // degerees by adding 360 degrees.
        constexpr double twopi = 2*M_PI;
        mAzimuth = M_PI/2 - std::atan2(u21, u31);
        if (mAzimuth < 0){mAzimuth = mAzimuth + twopi;}
        // To get to azimuth simply swing this around 180 degrees.  Note, that
        // since azimuth is in the range [0,360] that this can exceed the
        // convention that the backazimuth be in the range [0,360].  In that
        // case we add 360 degrees.
        mBackAzimuth = mAzimuth - M_PI;
        if (mBackAzimuth < 0){mBackAzimuth = mBackAzimuth + twopi;}
    }
//private:
    /// A copy of the input vertical channel.  This has dimension [mSamples].
    double *mVertical64f = nullptr;
    /// A copy of the input north channel.  This has dimension [mSamples].
    double *mNorth64f = nullptr;
    /// A copy of the input east channel.  This has dimension [mSamples].
    double *mEast64f = nullptr;
    /// A copy of the input vertical channel.  This has dimension [mSamples].
    float *mVertical32f = nullptr;
    /// A copy of the input north channel.  This has dimension [mSamples].
    float *mNorth32f = nullptr;
    /// A copy of the input east channel.  This has dimension [mSamples].
    float *mEast32f = nullptr;
    /// The cross-variance matrix.  This is a [3 x 3] symmetric matrix:
    /// [ zz zn ze ]
    /// [ zn nn ne ]
    /// [ ze ne ee ].
    std::array<double, 9> mCrossVarianceMatrix;
    /// The eigenavlues of the cross-variance matrix.  These are sorted in
    /// descending order and should all be positive.
    std::array<double, 3> mEigenvalues;
    /// The eigenvectors of the cross-variance matrix.  This is column major
    /// format and is ordered s.t. the first eigenvector (column) corresponds to
    /// the largest eigenvalue.
    /// [ u11 u21 u31 ]
    /// [ u12 u22 u32 ]
    /// [ u13 u23 u33 ]
    std::array<double, 9> mEigenvectors; 
    /// Angle of incidence in radians: acos(|u11|)
    double mIncidenceAngle = 0;
    /// Rectilinearity:  1 - (lambda2 + lambda3)/(2*lambda1)
    double mRectilinearity = 0;
    /// The apparent azimuth.  This is in the range [0,2*pi]
    double mAzimuth = 0;
    /// The apparent backazimuth.  This is in the range [0,2*pi]
    double mBackAzimuth = 0;
    //std::array<double, 32> mEigenWorkSpace;
    int mSamples = 0;
    /// Flag indicating that the signals are set
    bool mHaveSignals = false;
    /// Flag indicating that the eigenvalues are computed
    bool mHaveEigenvalues = false;
    /// Flag indicating that the class is initialized
    bool mInitialized = false;
    /// Precision of module
    RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
};

/// Default constructor
template<class T>
EigenPolarizer<T>::EigenPolarizer() :
    pImpl(std::make_unique<EigenPolarizerImpl> ())
{
}

/// Copy constructor
template<class T>
EigenPolarizer<T>::EigenPolarizer(const EigenPolarizer &polarizer)
{
    *this = polarizer;
}

/// Move constructor
template<class T>
EigenPolarizer<T>::EigenPolarizer(EigenPolarizer &&polarizer) noexcept
{
    *this = std::move(polarizer);
}

/// Copy assignment operator
template<class T>
EigenPolarizer<T>& EigenPolarizer<T>::operator=(const EigenPolarizer &polarizer)
{
    if (&polarizer == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<EigenPolarizerImpl> (*polarizer.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
EigenPolarizer<T>& 
EigenPolarizer<T>::operator=(EigenPolarizer &&polarizer) noexcept
{
    if (&polarizer == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(polarizer.pImpl);
    return *this;
}

/// Default destructor
template<class T>
EigenPolarizer<T>::~EigenPolarizer() = default;

/// Releases memory and resets class
template<class T>
void EigenPolarizer<T>::clear() noexcept
{
    pImpl->clear();
} 

/// Initialization
template<>
void EigenPolarizer<double>::initialize(const int nSamples)
{
    clear();
    if (nSamples < 1)
    {
        RTSEIS_THROW_IA("nSamples = %d must be positive", nSamples);
    }
    constexpr RTSeis::Precision precision = RTSeis::Precision::DOUBLE;
    pImpl->initialize(nSamples, precision);
}

template<>
void EigenPolarizer<float>::initialize(const int nSamples)
{
    clear();
    if (nSamples < 1)
    {
        RTSEIS_THROW_IA("nSamples = %d must be positive", nSamples);
    }
    constexpr RTSeis::Precision precision = RTSeis::Precision::FLOAT;
    pImpl->initialize(nSamples, precision);
}

/// Sets the input signals
template<class T>
void EigenPolarizer<T>::setSignals(const int nSamples,
                                   const T vertical[],
                                   const T north[],
                                   const T east[])
{
    pImpl->mHaveSignals = false;
    pImpl->mHaveEigenvalues = false;
    int nSamplesRef = getNumberOfSamples(); // Throws on initialization
    if (nSamplesRef != nSamples)
    {
        RTSEIS_THROW_IA("nSammples = %d must equal %d", nSamples, nSamplesRef);
    }
    if (vertical == nullptr || north == nullptr || east == nullptr)
    {
        if (vertical == nullptr){RTSEIS_THROW_IA("%s", "vertical is NULL");}
        if (north == nullptr){RTSEIS_THROW_IA("%s", "north is NULL");}
        RTSEIS_THROW_IA("%s", "east is NULL");
    }
    // Set signals, compute cross-power matrix, and its eigendecomposition
    pImpl->setSignalsAndDecompose(pImpl->mSamples, vertical, north, east);
}

/// Get the azimuth based on largest eigenvector
template<class T>
T EigenPolarizer<T>::getAzimuth(const bool wantRadians) const
{
    if (!haveInputSignals()){RTSEIS_THROW_RTE("%s", "Signals not yet set\n");}
    if (wantRadians){return pImpl->mAzimuth;}
    return pImpl->mAzimuth*(180/M_PI);
}

/// Get the backazimuth based on the largest eigenvector
template<class T>
T EigenPolarizer<T>::getBackAzimuth(const bool wantRadians) const
{
    if (!haveInputSignals()){RTSEIS_THROW_RTE("%s", "Signals not yet set\n");}
    if (wantRadians){return pImpl->mBackAzimuth;}
    return pImpl->mBackAzimuth*(180/M_PI);
}

/// Get the incidence angle based on the largest eigenvector
template<class T>
T EigenPolarizer<T>::getIncidenceAngle(const bool wantRadians) const
{
    if (!haveInputSignals()){RTSEIS_THROW_RTE("%s", "Signals not yet set\n");}
    if (wantRadians){return pImpl->mIncidenceAngle;}
    return pImpl->mIncidenceAngle*(180/M_PI);
}

/// Get the rectilinearity based on the eigenvalues
template<class T>
T EigenPolarizer<T>::getRectilinearity() const
{
    if (!haveInputSignals()){RTSEIS_THROW_RTE("%s", "Signals not yet set\n");}
    return pImpl->mRectilinearity;
}

/// Get number of samples
template<class T>
int EigenPolarizer<T>::getNumberOfSamples() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mSamples;
}

/// Are the input signals set?
template<class T>
bool EigenPolarizer<T>::haveInputSignals() const noexcept
{
    return pImpl->mHaveSignals;
}

/// Is the class is initialized?
template<class T>
bool EigenPolarizer<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Template instantiation
template class RTSeis::Utilities::Polarization::EigenPolarizer<double>;
template class RTSeis::Utilities::Polarization::EigenPolarizer<float>;
