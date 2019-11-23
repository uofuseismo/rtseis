#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <cmath>
#include <algorithm>
#include <exception>
#include <cassert>
#include <stdexcept>
#include <mkl.h>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/interpolation/linear.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::Interpolation;

class Linear::LinearImpl
{
public:
    /// Constructor
    LinearImpl()
    {
        std::memset(&mTask, 0, sizeof(DFTaskPtr));
    }
    /// Destructor
    ~LinearImpl()
    {
        clear();
    }
    /// Clears the memory
    void clear() noexcept
    {
        if (mHaveTask){dfDeleteTask(&mTask);}
        if (mSplineCoeffs){MKL_free(mSplineCoeffs);}
        if (mX){MKL_free(mX);}
        mSplineCoeffs = nullptr;
        mX = nullptr;
        mXMin = std::numeric_limits<double>::max();
        mXMax = std::numeric_limits<double>::min();
        mHaveTask = false;
        mInitialized = false;
    }
    /// Creates a new task for when x is uniform
    int initialize(const int npts,
                   const std::pair<double, double> x,
                   const double y[])
    {
        clear();
        MKL_INT nx = npts;
        mXMin = x.first;
        mXMax = x.second;
        constexpr MKL_INT ny = 1; // Function to interpolate is scalar
        mX = static_cast<double *> (MKL_calloc(npts, sizeof(double), 64));
        auto dx = (mXMax - mXMin)/static_cast<double> (npts - 1);
        #pragma omp simd
        for (auto i=0; i<npts; ++i)
        {
            mX[i] = mXMin + i*dx;
        }
        auto status = dfdNewTask1D(&mTask, nx, mX,
                                   DF_QUASI_UNIFORM_PARTITION,
                                   ny, y, DF_NO_HINT);
        mHaveTask = true;
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1;
        }
        // Set spline parameters in the data fitting task
        double *bc = nullptr;
        double *ic = nullptr;
        mSplineCoeffs = static_cast<double *>
                     (MKL_calloc(ny*mSplineOrder*(nx - 1), sizeof(double), 64));
        status = dfdEditPPSpline1D(mTask, mSplineOrder, mSplineType, mSplineBC,
                                   bc, mSplineIC, ic, mSplineCoeffs, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1;
        }
        // Construct the spline
        status = dfdConstruct1D(mTask, DF_PP_SPLINE, DF_METHOD_STD);
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1; 
        }
        mInitialized = true;
        return 0;
    }
    /// Creates a new task for when x is nonuniform
    int initialize(const int npts,
                   const double x[],
                   const double y[])
    {
        clear();
        MKL_INT nx = npts;
        mXMin = x[0];
        mXMax = x[npts-1];
        constexpr MKL_INT ny = 1; // Function to interpolate is scalar
        mX = static_cast<double *> (MKL_calloc(npts, sizeof(double), 64));
        std::copy(x, x+npts, mX);
        auto status = dfdNewTask1D(&mTask, nx, mX, DF_NO_HINT,
                                   ny, y, DF_NO_HINT);
        mHaveTask = true;
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1;
        }
        // Set spline parameters in the data fitting task
        double *bc = nullptr;
        double *ic = nullptr;
        mSplineCoeffs = static_cast<double *>
                     (MKL_calloc(ny*mSplineOrder*(nx - 1), sizeof(double), 64));
        status = dfdEditPPSpline1D(mTask, mSplineOrder, mSplineType, mSplineBC,
                                   bc, mSplineIC, ic, mSplineCoeffs, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1;
        }
        // Construct the spline
        status = dfdConstruct1D(mTask, DF_PP_SPLINE, DF_METHOD_STD);
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1;
        }
        mInitialized = true;
        return 0;
    }
    /// Interpolates over the uniform interval with many yq
    int interpolate(const int nq, 
                    const std::pair<double, double> xInterval,
                    double yq[])
    {
        double xq[2] = {xInterval.first, xInterval.second};
        const MKL_INT nsite = nq;
        constexpr MKL_INT nOrder = 1;  // Length of dorder
        const MKL_INT dOrder[1] = {0}; // Order of derivatives
        auto status = dfdInterpolate1D(mTask, DF_INTERP, DF_METHOD_PP,
                                       nsite, xq,
                                       DF_UNIFORM_PARTITION,
                                       nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yq, 
                                       DF_MATRIX_STORAGE_ROWS, NULL);
        if (status != DF_STATUS_OK){return -1;}
        return 0;  
    }
    /// Interpolates at select abscissa
    int interpolate(const int nq,
                    const double xq[],
                    double yq[],
                    const bool lsorted = false)
    {
        if (nq < 1){return 0;}
        const MKL_INT nsite = nq;
        const MKL_INT nOrder = 1;  // Length of dorder
        MKL_INT dOrder[1] = {0}; // Derivative order
        MKL_INT sortedHint = DF_SORTED_DATA;
        if (!lsorted){sortedHint = DF_NO_HINT;}
        auto status = dfdInterpolate1D(mTask, DF_INTERP, DF_METHOD_PP,
                                       nsite, xq,
                                       sortedHint, nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yq,
                                       DF_MATRIX_STORAGE_ROWS, NULL);
        if (status != DF_STATUS_OK){return -1;}
        return 0;
    }
//private:
    /// MKL data fitting task
    DFTaskPtr mTask;
    /// Copy of x data points.  Unclear if MKL copies x.  If not, then x
    /// can go out of scope or be deleted prior to usage.
    double *mX = nullptr;
    /// The spline coefficients.
    double *mSplineCoeffs = nullptr;
    /// The minimum value of x
    double mXMin = std::numeric_limits<double>::max();
    /// The maximum value of x 
    double mXMax = std::numeric_limits<double>::min();
    /// Appropriate parameters for MKL
    const MKL_INT mSplineOrder = DF_PP_LINEAR;
    const MKL_INT mSplineType  = DF_PP_LINEAR;
    const MKL_INT mSplineBC    = DF_NO_BC;
    const MKL_INT mSplineIC    = DF_NO_IC;
    /// Note if I have created the task and must therefore delete it.
    bool mHaveTask = false;
    /// Determines if the class is initialized.
    bool mInitialized = false;
};

/// Constructor
Linear::Linear() :
    pImpl(std::make_unique<LinearImpl> ())
{
}

/// Move constructor
Linear::Linear(Linear &&linear) noexcept
{
    *this = std::move(linear);
}

/// Move assignment operator
Linear& Linear::operator=(Linear &&linear) noexcept
{
    if (&linear == this){return *this;}
    if (pImpl){pImpl->clear();}
    pImpl = std::move(linear.pImpl);
    return *this;
}

/// Destructor
Linear::~Linear() = default;

/// Clears the memory
void Linear::clear() noexcept
{
    pImpl->clear();
}

/// Checks if the class is initialized 
bool Linear::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Initialize the linear interpolant for regularly spaced points
void Linear::initialize(const int npts,
                        const std::pair<double, double> xInterval,
                        const double y[])
{
    clear();
    if (npts < 2)
    {
        RTSEIS_THROW_IA("npts = %d must be at least 2", npts);
    }
    if (y == nullptr){RTSEIS_THROW_RTE("%s", "y is NULL");}
    if (xInterval.first >= xInterval.second)
    {
        RTSEIS_THROW_IA("x.first = %lf must be less than x.second = %lf",
                        xInterval.first, xInterval.second); 
    }
    // Initialize
    auto ierr = pImpl->initialize(npts, xInterval, y);
    if (ierr != 0)
    {
        RTSEIS_THROW_RTE("%s", "Failed to initialize interpolator");
    }
}

/// Initialize the linear interpolant for irregularly spaced points
void Linear::initialize(const int npts,
                        const double x[],
                        const double y[])
{
    clear();
    if (npts < 2)
    {   
        RTSEIS_THROW_IA("npts = %d must be at least 2", npts);
    }   
    if (x == nullptr){RTSEIS_THROW_RTE("%s", "x is NULL");}
    if (y == nullptr){RTSEIS_THROW_RTE("%s", "y is NULL");}
    // Assert that the data is increasing order
    int isInvalid = 0;
    #pragma omp simd reduction(+:isInvalid)
    for (auto i=0; i<npts-1; ++i)
    {   
        if (x[i+1] <= x[i]){isInvalid = isInvalid + 1;} 
    }   
    if (isInvalid > 0)
    {   
        RTSEIS_THROW_RTE("%s", "At least one x[i+1] <= x[i]\n");
    }    
    auto ierr = pImpl->initialize(npts, x, y);
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Interpolation failed\n");} 
}

/// Interpolate on regularly spaced points
void Linear::interpolate(const int nq,
                         const std::pair<double, double> xq,
                         double *yqIn[]) const
{
    if (nq <= 0){return;} // Nothing to do
    // Start checks
    auto xmin = getMinimumX(); // Throws on uninitialized
    auto xmax = getMaximumX();
    double *yq = *yqIn;
    if (yq == nullptr){RTSEIS_THROW_IA("%s", "yq is NULL");}
    if (xq.first >= xq.second)
    {
        RTSEIS_THROW_IA("xq.first = %lf must be less than xq.second  %lf\n",
                        xq.first, xq.second);
    }
    if (xq.first < xmin || xq.second > xmax)
    {
       RTSEIS_THROW_IA("Min/max of xq = (%lf,%lf) must be in range [%lf,%lf]",
                       xq.first, xq.second, xmin, xmax);
    }
    auto ierr = pImpl->interpolate(nq, xq, yq);
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Interpolation failed\n");}
}

/// Interpolate at a bunch of points
void Linear::interpolate(const int nq, const double xq[],
                         double *yqIn[]) const
{
    if (nq <= 0){return;} // Nothing to do
    // Now start performing checks
    auto xmin = getMinimumX(); // Throws on unitialized error
    auto xmax = getMaximumX();
    double *yq = *yqIn;
    if (xq == nullptr || yq == nullptr)
    {
        if (xq == nullptr){RTSEIS_THROW_IA("%s", "xq is NULL");}
        RTSEIS_THROW_IA("%s", "yq is NULL");
    }
    double xqMin, xqMax;
    ippsMinMax_64f(xq, nq, &xqMin, &xqMax);
    if (xqMin < xmin || xqMax > xmax)
    {
       RTSEIS_THROW_IA("Min/max of xq = (%lf,%lf) must be in range [%lf,%lf]",
                       xqMin, xqMax, xmin, xmax);
    }
    auto ierr = pImpl->interpolate(nq, xq, yq);
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Interpolation failed\n");}
}
/// Get minimum x for interpolation
double Linear::getMinimumX() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mXMin;
}

/// Get maximum x for interpolation
double Linear::getMaximumX() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mXMax;
}
