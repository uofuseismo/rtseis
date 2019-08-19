#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <mkl/lapacke.h>
#include <ipps.h>
#include "rtseis/utilities/polarizer/eigenPolarizer.hpp"

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
    crossVarianceMatrix[0] = zz;
    crossVarianceMatrix[1] = zn;
    crossVarianceMatrix[2] = ze;
    crossVarianceMatrix[3] = zn;
    crossVarianceMatrix[4] = nn;
    crossVarianceMatrix[5] = ne;
    crossVarianceMatrix[6] = ze;
    crossVarianceMatrix[7] = ne;
    crossVarianceMatrix[8] = ee;
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
        double zDemean = static_cast<double> (z[i] - zMean);
        double nDemean = static_cast<double> (n[i] - nMean);
        double eDemean = static_cast<double> (e[i] - eMean);
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
    crossVarianceMatrix[0] = zz;
    crossVarianceMatrix[1] = zn;
    crossVarianceMatrix[2] = ze;
    crossVarianceMatrix[3] = zn;
    crossVarianceMatrix[4] = nn;
    crossVarianceMatrix[5] = ne;
    crossVarianceMatrix[6] = ze;
    crossVarianceMatrix[7] = ne;
    crossVarianceMatrix[8] = ee;
}

/// Computes eigendecomposition of 3 x 3 symmetric matrix
void computeEigendecomposition(const double A[9], double eig[3], double X[9])
{
    constexpr int lwork = 102;
    constexpr int lda = 3;
    std::array<double, lwork> work;
    std::copy(A, A+9, X);
#ifdef DEBUG
    auto info =
#endif
    LAPACKE_dsyev_work(LAPACK_COL_MAJOR, // Matrix of eigenvectors is col major
                       "V",              // Eigenvalues and eigenvectors
                        X,               // Will hold the eigenvectors
                        lda,             // Leading dimension is 3
                        eig,             // Eigenvalues
                        work.data(),     // Workspace
                        lwork);          // Size of workspace
    
#ifdef DEBUG
    assert(info == 0);
#endif             
}

}

class EigenPolarizer::EigenPolarizerImpl
{
    std::array<double, 9> mCrossVarianceMatrix; 
    std::array<double, 32> mEigenWorkSpace;
}


