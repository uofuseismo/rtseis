#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <algorithm>
#include <cfloat>
#include <mkl.h>
#include "rtseis/enums.hpp"
#include "rtseis/utilities/interpolation/weightedAverageSlopes.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"
#include "private/throw.hpp"
#include <ipps.h>

using namespace RTSeis::Utilities::Interpolation;

namespace
{

#pragma omp declare simd
template<class T>
inline void computeWiWiMi(const T mi, T *wi, T *wimi)
{
    constexpr T huge = std::numeric_limits<T>::max();
    constexpr T eps = std::numeric_limits<T>::epsilon();
    // The goal is to compute w_i and w_i m_i = 1/max(|m_i|, epsilon)*m_i
    // The trick is that to let w_i m_i safely go to 0 and let the
    // denominator become really big when the slope is 0.
    *wi = huge;
    *wimi = 0;
    auto ami = std::abs(mi);
    if (ami > eps)
    {
        *wi = 1/ami;
        // "If m_i and m_{i+1} have opposite signs then s_i is zero".
        // Realize, the copysign in the ensuing numerator will ensure
        // operations like 1 - 1, 1 + 1, -1 + 1, or -1 - 1.
        *wimi = std::copysign(1, mi);
    }
}

/// Implements first equation in Wiggins' Interpolation fo Digitized Curves
/// pg. 2077
template<typename T>
void computeUniformSlopes(const int npts, const T dx,
                          const T *__restrict__ y,
                          T *__restrict__ slopes,
                          T *__restrict__ splineCoeffs)
{
    const T one = 1;
    T dxi = one/dx;
    // Handle the initial conditions
    slopes[0] = (y[1] - y[0])/dx;
    #pragma omp simd
    for (int i=1; i<npts-1; ++i)
    {
        T mi  = (y[i] - y[i-1])*dxi;
        T mi1 = (y[i+1] - y[i])*dxi;
        // w_i and w_i m_i = 1/max(|m_i|, epsilon)*m_i
        T wi, wimi;
        computeWiWiMi(mi, &wi, &wimi);
        // w_{i+1} and w_{i+1} m_{i+1} = 1/(max(|m_{i+1}|, epsilon)*m_{i+1}
        T wi1, wi1mi1;
        computeWiWiMi(mi1, &wi1, &wi1mi1);
        // s_i = (w_i*m_i + w_{i+1}*m_{i+1})/(w_{i} + w_{i+1})
        slopes[i] = (wimi + wi1mi1)/(wi + wi1);
    }
    // Handle the final conditions
    slopes[npts-1] = (y[npts-1] - y[npts-2])/dx;
    // Compute the spline coefficients: Eqn 4 from
    // Monotone Piecewise Cubic Interpolation - Fritsch and Carlson 1980
    auto dxi2 = one/(dx*dx);
    for (int i=0; i<npts-1; ++i)
    {
        auto di = slopes[i];
        auto di1 = slopes[i+1];
        auto delta = (y[i+1] - y[i])*dxi;
        auto c1 = y[i];
        auto c2 = di;
        auto c3 = (-2*di - di1 + 3*delta)*dxi;
        auto c4 = (di + di1 - 2*delta)*dxi2;
        splineCoeffs[4*i+0] = c1; 
        splineCoeffs[4*i+1] = c2;
        splineCoeffs[4*i+2] = c3;
        splineCoeffs[4*i+3] = c4;
    }
}

template<typename T>
void computeUniformSlopes(const int npts,
                          const T *__restrict__ x,
                          const T *__restrict__ y,
                          T *__restrict__ slopes,
                          T *__restrict__ splineCoeffs)
{
    // Handle the initial conditions
    slopes[0] = (y[1] - y[0])/(x[1] - x[0]);
    #pragma omp simd
    for (int i=1; i<npts-1; ++i)
    {
        T mi  = (y[i] - y[i-1])/(x[i] - x[i-1]);
        T mi1 = (y[i+1] - y[i])/(x[i+1] - x[i]);
        // w_i and w_i m_i = 1/max(|m_i|, epsilon)*m_i
        T wi, wimi;
        computeWiWiMi(mi, &wi, &wimi);
        // w_{i+1} and w_{i+1} m_{i+1} = 1/(max(|m_{i+1}|, epsilon)*m_{i+1}
        T wi1, wi1mi1;
        computeWiWiMi(mi1, &wi1, &wi1mi1);
        // s_i = (w_i*m_i + w_{i+1}*m_{i+1})/(w_{i} + w_{i+1})
        slopes[i] = (wimi + wi1mi1)/(wi + wi1);
    }
    // Handle the final conditions
    slopes[npts-1] = (y[npts-1] - y[npts-2])/(x[npts-1] - x[npts-2]);
    // Compute the spline coefficients: Eqn 4 from
    // Monotone Piecewise Cubic Interpolation - Fritsch and Carlson 1980
    const T one = 1;
    for (int i=0; i<npts-1; ++i)
    {
        T dx = x[i+1] - x[i]; 
        T dxi = one/dxi; 
        auto dxi2 = dxi*dxi;
        auto di = slopes[i];
        auto di1 = slopes[i+1];
        auto delta = (y[i+1] - y[i])*dxi;
        auto c1 = y[i];
        auto c2 = di;
        auto c3 = (-2*di - di1 + 3*delta)*dxi;
        auto c4 = (di + di1 - 2*delta)*dxi2;
        splineCoeffs[4*i+0] = c1;
        splineCoeffs[4*i+1] = c2;
        splineCoeffs[4*i+2] = c3;
        splineCoeffs[4*i+3] = c4;
    }
}

}

template<class T>
class WeightedAverageSlopes<T>::WeightedAverageSlopesImpl
{
public:
    /// Copy assignment operator
    WeightedAverageSlopesImpl&
    operator=(const WeightedAverageSlopesImpl &slopes)
    {
        if (&slopes == this){return *this;}
        clear();
        mXiEqual64f[0] = slopes.mXiEqual64f[0];
        mXiEqual64f[1] = slopes.mXiEqual64f[1];
        mXiEqual32f[0] = slopes.mXiEqual32f[0];
        mXiEqual32f[1] = slopes.mXiEqual32f[1];
        mRange = slopes.mRange;
        mCoeffs = slopes.mCoeffs;
        mSites = slopes.mSites;
        mUniformPartition = slopes.mUniformPartition;
        mPrecision = slopes.mPrecision;
        mInitialized = slopes.mInitialized; 
        if (mCoeffs > 0)
        {
            if (mPrecision == RTSeis::Precision::DOUBLE) 
            {
                mSplineCoeffs64f = ippsMalloc_64f(mCoeffs);
                ippsCopy_64f(slopes.mSplineCoeffs64f, mSplineCoeffs64f,
                             mCoeffs);
            }
            else
            {
                mSplineCoeffs32f = ippsMalloc_32f(mCoeffs);
                ippsCopy_32f(slopes.mSplineCoeffs32f, mSplineCoeffs32f,
                             mCoeffs);
            }
        }
        if (!mUniformPartition && mSites > 0)
        {
            if (mPrecision == RTSeis::Precision::DOUBLE)
            {
                mXi64f = ippsMalloc_64f(mSites);
                ippsCopy_64f(slopes.mXi64f, mXi64f, mSites);
            } 
            else
            {
                mXi32f = ippsMalloc_32f(mSites);
                ippsCopy_32f(slopes.mXi32f, mXi32f, mSites);
            }
        }
        return *this;
    }
    /// Destructor
    ~WeightedAverageSlopesImpl()
    {
        clear(); 
    }
    /// Releases memory on the class
    void clear() noexcept
    {
        if (mSplineCoeffs64f){ippsFree(mSplineCoeffs64f);}
        if (mSplineCoeffs32f){ippsFree(mSplineCoeffs32f);}
        if (mXi64f){ippsFree(mXi64f);}
        if (mXi32f){ippsFree(mXi32f);}
        if (mHaveTask64f){dfDeleteTask(&mTask64f);}
        if (mHaveTask32f){dfDeleteTask(&mTask32f);}
        mSplineCoeffs64f = nullptr;
        mSplineCoeffs32f = nullptr;
        mXi64f = nullptr;
        mXi32f = nullptr;
        mTask64f = nullptr;
        mTask32f = nullptr;
        mXiEqual64f[0] = 0;
        mXiEqual64f[1] = 0;
        mXiEqual32f[0] = 0;
        mXiEqual32f[1] = 0;
        mRange = std::make_pair(0, 0);
        mCoeffs = 0;
        mSites = 0;
        mPrecision = RTSeis::Precision::DOUBLE;
        mUniformPartition = true;
        mHaveTask64f = false;
        mHaveTask32f = false;
        mInitialized = false;        
    }
//private: 
    /// The task pointers for the data fitting 
    DFTaskPtr mTask64f = nullptr; 
    DFTaskPtr mTask32f = nullptr; 
    /// The spline coefficients.  This has dimension [mCoeffs].
    double *mSplineCoeffs64f = nullptr;
    /// The abscissas for non-uniform inteprolation. This has dimension [mSites]
    double *mXi64f = nullptr;
    /// The spline coefficients.  This has dimension [mCoeffs].
    float *mSplineCoeffs32f = nullptr;
    /// The abscissas for non-uniform inteprolation. This has dimension [mSites]
    float *mXi32f = nullptr;
    /// The range for an equal interpolation
    double mXiEqual64f[2] = {0, 0};
    float mXiEqual32f[2] = {0, 0};
    /// The min/max x for interpolation
    std::pair<double, double> mRange = std::make_pair(0, 0);
    /// The number of spline coefficients.  This is splineOrder*(mSites - 1)
    int mCoeffs = 0;
    /// The number of interpolation sites
    int mSites = 0;
    /// The custom spline order
    const MKL_INT splineOrder = 4;
    /// Notes the precision of the module
    RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
    /// Using a uniform partition 
    bool mUniformPartition = true;
    /// Have the double MKL task.
    bool mHaveTask64f = false;
    /// Have the float MKL task.
    bool mHaveTask32f = false;
    /// Class is initialized.
    bool mInitialized = false;
};

/// Constructors
template<class T>
WeightedAverageSlopes<T>::WeightedAverageSlopes() :
    pImpl(std::make_unique<WeightedAverageSlopesImpl> ())
{
}

/// Copy constructor
template<class T>
WeightedAverageSlopes<T>::WeightedAverageSlopes(
    const WeightedAverageSlopes &slopes)
{
    *this = slopes;
}

/// Move constructor
template<class T>
WeightedAverageSlopes<T>::WeightedAverageSlopes(
    WeightedAverageSlopes &&slopes) noexcept
{
    *this = std::move(slopes);
}

/// Copy assignment operator
template<class T>
WeightedAverageSlopes<T>& WeightedAverageSlopes<T>::operator=(
    const WeightedAverageSlopes &slopes)
{
    if (&slopes == this){return *this;}
    pImpl = std::make_unique<WeightedAverageSlopesImpl> (*slopes.pImpl);
    return *this;
}

/// Move assignment operator
template<class T>
WeightedAverageSlopes<T>& WeightedAverageSlopes<T>::operator=(
    WeightedAverageSlopes &&slopes) noexcept
{
    if (&slopes == this){return *this;}
    pImpl = std::move(slopes.pImpl);
    return *this;
}

/// Destructors
template<class T>
WeightedAverageSlopes<T>::~WeightedAverageSlopes() = default;

/// Releases memory on the module
template<class T>
void WeightedAverageSlopes<T>::clear() noexcept
{
    pImpl->clear();
}

/// Determines if the class is initialized
template<class T>
bool WeightedAverageSlopes<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Initialize the class
template<>
void WeightedAverageSlopes<double>::initialize(
    const int npts,
    const std::pair<double, double> x,
    const double y[])
{
    clear();
    if (npts < 2)
    {
        throw std::invalid_argument("npts = " + std::to_string(npts)
                                  + " must be at least 2");
    }
    if (x.first >= x.second)
    {
        throw std::invalid_argument("x.first = " + std::to_string(x.first)
                                  + " must be less than x.second = "
                                  + std::to_string(x.second));
    }
    if (y == nullptr){throw std::invalid_argument("y is NULL");}
    // Compute the spline coefficients
    auto dx = (x.second - x.first)/static_cast<double> (npts - 1);
    pImpl->mSites = npts;
    pImpl->mCoeffs = pImpl->splineOrder*(pImpl->mSites - 1);
    pImpl->mSplineCoeffs64f = ippsMalloc_64f(pImpl->mCoeffs);
    auto slopes = ippsMalloc_64f(npts); // Workspace
    computeUniformSlopes(npts, dx, y, slopes, pImpl->mSplineCoeffs64f);
    ippsFree(slopes);
    // Create a custom piecewise 4th order spline
    pImpl->mTask64f = nullptr;
    pImpl->mRange.first = x.first;
    pImpl->mRange.second = x.second;
    pImpl->mXiEqual64f[0] = x.first;
    pImpl->mXiEqual64f[1] = x.second;
    auto status = dfdNewTask1D(&pImpl->mTask64f, npts, pImpl->mXiEqual64f,
                               DF_UNIFORM_PARTITION, 1, y, DF_NO_HINT); 
    if (status != DF_STATUS_OK)
    {
        dfDeleteTask(&pImpl->mTask64f);
        throw std::runtime_error("Failed to create task\n");
    }
    status = dfdEditPPSpline1D(pImpl->mTask64f, pImpl->splineOrder,
                               DF_PP_DEFAULT, DF_NO_BC, NULL, DF_NO_IC, NULL,
                               pImpl->mSplineCoeffs64f, DF_NO_HINT);
    if (status != DF_STATUS_OK)
    {
        dfDeleteTask(&pImpl->mTask64f);
        throw std::runtime_error("Failed to edit spline pipeline\n");
    }
    pImpl->mHaveTask64f = true;
    pImpl->mInitialized = true;
}

/// Initialize the class
template<>
void WeightedAverageSlopes<double>::initialize(
    const int npts,
    const double x[],
    const double y[])
{
    clear();
    if (npts < 2)
    {
        throw std::invalid_argument("npts = " + std::to_string(npts)
                                  + " must be at least 2");
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    if (!std::is_sorted(x, x+npts))
    {
        throw std::invalid_argument("x is not sorted");
    }
    pImpl->mSites = npts;
    pImpl->mCoeffs = pImpl->splineOrder*(pImpl->mSites - 1);
    pImpl->mSplineCoeffs64f = ippsMalloc_64f(pImpl->mCoeffs);
    auto slopes = ippsMalloc_64f(npts); // Workspace
    computeUniformSlopes(npts, x, y, slopes, pImpl->mSplineCoeffs64f);
    ippsFree(slopes);
    // Create a custom piecewise 4th order spline
    pImpl->mTask64f = nullptr;
    pImpl->mRange.first = x[0];
    pImpl->mRange.second = x[npts-1];
    pImpl->mXiEqual64f[0] = x[0];
    pImpl->mXiEqual64f[1] = x[npts-1];
    pImpl->mXi64f = ippsMalloc_64f(npts);
    ippsCopy_64f(x, pImpl->mXi64f, npts);
    auto status = dfdNewTask1D(&pImpl->mTask64f, npts, pImpl->mXi64f,
                               DF_NON_UNIFORM_PARTITION, 1, y, DF_NO_HINT);
    if (status != DF_STATUS_OK)
    {
        dfDeleteTask(&pImpl->mTask64f);
        throw std::runtime_error("Failed to create task\n");
    }
    status = dfdEditPPSpline1D(pImpl->mTask64f, pImpl->splineOrder,
                               DF_PP_DEFAULT, DF_NO_BC, NULL, DF_NO_IC, NULL,
                               pImpl->mSplineCoeffs64f, DF_NO_HINT);
    if (status != DF_STATUS_OK)
    {
        dfDeleteTask(&pImpl->mTask64f);
        throw std::runtime_error("Failed to edit spline pipeline\n");
    }
    pImpl->mHaveTask64f = true;
    pImpl->mInitialized = true;
}


// Interpolate
template<>
void WeightedAverageSlopes<double>::interpolate(
    const int nq, const double xq[], double *yqIn[]) const
{
    // Checks
    if (nq < 1){return;} // Nothing to do
    double xMin = getMinimumX(); // Throws on initialization
    double xMax = getMaximumX(); // Throws on initialization
    double *yq = *yqIn;
    if (xq == nullptr || yq == nullptr)
    {
        if (xq == nullptr){throw std::invalid_argument("xq is NULL");}
        throw std::invalid_argument("yq is NULL");
    }
    double xqMin, xqMax;
    ippsMinMax_64f(xq, nq, &xqMin, &xqMax);
    if (xqMin < xMin || xqMax > xMax)
    {
        throw std::invalid_argument("Min/max of xq = ("
                                  + std::to_string(xqMin) + ","
                                  + std::to_string(xqMax) 
                                  + ") Must be in range ["
                                  + std::to_string(xMin) + ","
                                  + std::to_string(xMax) + "]");
    }
    // Check this is sorted
    bool lsorted = Math::VectorMath::isSorted(nq, xq);
    // Interpolate
    const MKL_INT nsite = nq;
    MKL_INT sortedHint = DF_SORTED_DATA;
    if (!lsorted){sortedHint = DF_NO_HINT;}
    constexpr MKL_INT nOrder = 1;  // Length of dorder
    const MKL_INT dOrder[1] = {0}; // Order of derivatives
    auto status = dfdInterpolate1D(pImpl->mTask64f, DF_INTERP, DF_METHOD_PP,
                                   nsite, xq,
                                   sortedHint, nOrder, dOrder,
                                   DF_NO_APRIORI_INFO, yq,
                                   DF_MATRIX_STORAGE_ROWS, NULL);
    if (status != DF_STATUS_OK)
    {
        throw std::runtime_error("Interpolation failed");
    }
}

// Uniform interpolation
template<>
void WeightedAverageSlopes<double>::interpolate(
    const int nq, const std::pair<double, double> xInterval,
    double *yqIn[]) const
{
    // Checks
    if (nq < 1){return;} // Nothing to do
    double xMin = getMinimumX(); // Throws on initialization
    double xMax = getMaximumX(); // Throws on initialization
    double *yq = *yqIn;
    if (yq == nullptr){RTSEIS_THROW_IA("%s", "yq is NULL");}
    double xqMin = xInterval.first;
    double xqMax = xInterval.second;
    if (xqMin > xqMax)
    {
        RTSEIS_THROW_IA("xInterval.first = %lf > xInterval.second = %lf",
                         xqMin, xqMax);
    }
    if (xqMin < xMin || xqMax > xMax)
    {
        RTSEIS_THROW_IA("Min/max of xq = (%lf,%lf) must be in range [%lf,%lf]",
                        xqMin, xqMax, xMin, xMax);
    }
    // Interpolate
    const MKL_INT nsite = nq;
    MKL_INT sortedHint = DF_UNIFORM_PARTITION;
    constexpr MKL_INT nOrder = 1;  // Length of dorder
    const MKL_INT dOrder[1] = {0}; // Order of derivatives
    double xq[2] = {xqMin, xqMax};
    auto status = dfdInterpolate1D(pImpl->mTask64f, DF_INTERP, DF_METHOD_PP,
                                   nsite, xq,
                                   sortedHint, nOrder, dOrder,
                                   DF_NO_APRIORI_INFO, yq,
                                   DF_MATRIX_STORAGE_ROWS, NULL);
    if (status != DF_STATUS_OK)
    {
        throw std::runtime_error("Interpolation failed");
    }
}

/// Get minimum x
template<class T>
double WeightedAverageSlopes<T>::getMinimumX() const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    return pImpl->mRange.first;
}

/// Get maximum x
template<class T>
double WeightedAverageSlopes<T>::getMaximumX() const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Class is not initialized\n");
    }
    return pImpl->mRange.second;
}

/// Template class instantiation
template class RTSeis::Utilities::Interpolation::WeightedAverageSlopes<double>;
//template class RTSeis::Utilities::Interpolation::WeightedAverageSlopes<float>;
