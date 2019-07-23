#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>
#include <mkl.h>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/interpolation/cubicSpline.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"

using namespace RTSeis::Utilities::Interpolation;

class CubicSpline::CubicSplineImpl
{
public:
    /*
    // TODO copy ctor requires me saving x and y.   After that, I can
    // recycle the spline coefficients and skip the "dfdConstruct" call.
    /// Copy assignment operator
    CubicSplineImpl& operator=(const CubicSplineImpl &spline)
    {
        if (&spline == this){return *this;}
        clear();
        if (!pImpl->mInititalized){return *this;} // Nothing to copy
        // Copy
        mSplineCoeffs = pImpl->mSplineCoeffs;
        mX = pImpl->mX;
        mXMin = pImpl->mXMin;
        mXMax = pImpl->mXMax;
        mSplineBC = pImpl->mSplineBC;
        mSplineIC = pImpl->mSplineIC;
        mBCType = pImpl->mBCType;
        mHaveTask = pImpl->mHaveTask;
        mInitialized = pImpl->mInitialized;
        // Set teh spline
        if (pImpl->mHaveTask)
        {
            
        }
        return *this;
    }
    */
    /// Destructor
    ~CubicSplineImpl()
    {
        clear();
    }
    /// Creates a new task for when x is uniform
    int createUniformXTask(const int npts,
                           const std::pair<double, double> xInterval,
                           const double y[])
    {
        deleteTask();
        mXMin = xInterval.first;
        mXMax = xInterval.second;
        mX.resize(npts);
        double dx = (mXMax - mXMin)/static_cast<double> (npts - 1);
        #pragma omp simd
        for (auto i=0; i<npts; ++i)
        {
            mX[i] = mXMin + i*dx;
        }
        const MKL_INT nx = npts; // Length of x
        const MKL_INT ny = 1;    // Dimension of vector valued function
        mSplineCoeffs.resize(ny*mSplineOrder*(nx - 1));
        std::fill(mSplineCoeffs.begin(), mSplineCoeffs.end(), 0);
        auto status = dfdNewTask1D(&mTask, nx, mX.data(),
                                   DF_QUASI_UNIFORM_PARTITION,
                                   ny, y, DF_NO_HINT);
        if (status != DF_STATUS_OK){return -1;}
        mHaveTask = true;
        return 0;
    }
    /// Creates a new task for when x is not necessarily uniofrm 
    int createTask(const int npts,
                   const double x[],
                   const double y[])
    {
        deleteTask();
        mX.resize(npts);
        std::copy(x, x+npts, mX.data());
        mXMin = x[0];
        mXMax = x[npts-1]; 
        const MKL_INT nx = npts;
        const MKL_INT ny = 1; // Dimension of vector valued function
        mSplineCoeffs.resize(ny*mSplineOrder*(nx - 1));
        auto status = dfdNewTask1D(&mTask, nx, mX.data(), DF_NO_HINT,
                                   ny, y, DF_NO_HINT);
        if (status != DF_STATUS_OK){return -1;}
        mHaveTask = true;
        return 0;
    }
    /// Sets the spline type and boundary conditions on the task
    int editPipeline(const CubicSplineBoundaryConditionType boundaryCondition)
    {
        // Figure out the boundary conditions first
        const double bcArray[2] = {0, 0};
        const double *bcs = NULL;
        const double *ics = NULL; 
        if (boundaryCondition == CubicSplineBoundaryConditionType::NATURAL)
        {
            mSplineBC = DF_BC_FREE_END;
            mSplineIC = DF_NO_IC;
        }
        else if (boundaryCondition == CubicSplineBoundaryConditionType::CLAMPED)
        {
            bcs = bcArray; // First derivatives are zero
            mSplineBC = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;
            mSplineIC = DF_NO_IC;
        }
        else if (boundaryCondition == CubicSplineBoundaryConditionType::NOT_A_KNOT)
        {
            mSplineBC = DF_BC_NOT_A_KNOT;
            mSplineIC = DF_NO_IC;
        }
        else if (boundaryCondition == CubicSplineBoundaryConditionType::PERIODIC)
        {
            mSplineBC = DF_BC_PERIODIC;
            mSplineIC = DF_NO_IC;
        }
        auto status = dfdEditPPSpline1D(mTask,
                                        mSplineOrder, // DF_PP_CUBIC
                                        mSplineType,  // DF_PP_NATURAL
                                        mSplineBC,
                                        bcs,
                                        mSplineIC,
                                        ics,
                                        mSplineCoeffs.data(),
                                        DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1;
        }
        return 0;
    }
    /// Makes the spline
    int constructSpline()
    {
        auto status = dfdConstruct1D(mTask, DF_PP_SPLINE, DF_METHOD_STD);
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1;
        }
        return 0;
    }
    /// Searches the cells for the interpolant points
    int searchCells(const int nq, const double xq[],
                    MKL_INT cell[],
                    const bool lsorted = false)
    {
        if (nq < 1){return 0;} // Nothing to do
        const MKL_INT nsite = nq;
        MKL_INT sortedHint = DF_SORTED_DATA;
        if (!lsorted){sortedHint = DF_NO_HINT;}
        auto status = dfdSearchCells1D(mTask, DF_METHOD_STD, nsite, xq,
                                       sortedHint, DF_NO_APRIORI_INFO,
                                       cell);
        if (status != DF_STATUS_OK){return -1;}
        return 0;
    }
    /// Interpolates over the uniform interval
    int interpolate(const int nq,
                    const std::pair<double, double> xInterval,
                    double yq[])
    {
        double xq[2] = {xInterval.first, xInterval.second};
        const MKL_INT nsite = nq;
        constexpr MKL_INT nOrder = 1;  // Length of dorder
        const MKL_INT dOrder[1] = {0}; // Order of derivatives
        auto status = dfdInterpolate1D(mTask, DF_INTERP, DF_METHOD_PP,
                                       nsite, xq, DF_UNIFORM_PARTITION,
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
    /// Integrates over the interval
    int integrate(const int nLimIn,
                  const double llim[], const double rlim[],
                  double integral[])
    {
        constexpr MKL_INT blockSize = 24; // > 48 seems to trigger bug
        const MKL_INT nLim = nLimIn;
        for (auto block=0; block<nLim; block=block+blockSize)
        {
            MKL_INT nInt = std::min(blockSize, nLim - block);
            auto status = dfdIntegrate1D(mTask, DF_METHOD_PP, nInt,
                                         &llim[block], DF_NO_HINT,
                                         &rlim[block], DF_NO_HINT,
                                         DF_NO_APRIORI_INFO,
                                         DF_NO_APRIORI_INFO,
                                         &integral[block], DF_NO_HINT);
            if (status != DF_STATUS_OK){return -1;}
        }
        return 0;
    }
    int integrate(const double lLim, const double rLim, double *integral)
    {
        constexpr MKL_INT nLim = 1;
        const double llim[1] = {lLim};
        const double rlim[1] = {rLim};
        auto status = dfdIntegrate1D(mTask, DF_METHOD_PP, nLim,
                                     llim, DF_UNIFORM_PARTITION,
                                     rlim, DF_UNIFORM_PARTITION,
                                     DF_NO_APRIORI_INFO, DF_NO_APRIORI_INFO,
                                     integral, DF_NO_HINT);
        if (status != DF_STATUS_OK){return -1;}
        return 0;
    }
         
    /// Deletes the task
    void deleteTask() noexcept
    {
        if (mHaveTask){dfDeleteTask(&mTask);}
        mHaveTask = false;
    }
    /// Clears the module
    void clear() noexcept
    {
        deleteTask();
        mSplineCoeffs.clear();
        mX.clear();
        mXMin =-DBL_MAX;
        mXMax = DBL_MAX;
        mSplineBC = DF_BC_FREE_END;
        mSplineIC = DF_NO_IC;
        mBCType = CubicSplineBoundaryConditionType::NATURAL;
        mHaveTask = false;
        mInitialized = false;
    }
///private:
    /// The data fitting task
    DFTaskPtr mTask;
    /// The spline coefficients
    std::vector<double> mSplineCoeffs;
    /// Copy of x data points.  Unclear if MKL copies x.  If not, then x
    /// can go out of scope or be deleted prior to usage.
    std::vector<double> mX;
    /// The minimum value of mX
    double mXMin =-DBL_MAX;
    /// The maximum value of mX
    double mXMax = DBL_MAX;    
    /// Default to natural cubic spline
    MKL_INT mSplineBC = DF_BC_FREE_END;
    /// Default to no initial conditions
    MKL_INT mSplineIC = DF_NO_IC;
    /// Some spline specific data fitting parameters
    const MKL_INT mSplineType = DF_PP_NATURAL;  // Deal with type in b.c.'s
    const MKL_INT mSplineOrder = DF_PP_CUBIC;   // This is for cubic splines
    const MKL_INT mSplineFormat = DF_PP_SPLINE; // Only supported option
    const MKL_INT mMethod = DF_METHOD_STD;      // Only supported option
    /// Defines my boundary conditions
    CubicSplineBoundaryConditionType mBCType
        = CubicSplineBoundaryConditionType::NATURAL;
    /// Lets me know if I can safely delete the task
    bool mHaveTask = false;
    /// Indicates that the class is initialized
    bool mInitialized = false;
};

/// Default constructor
CubicSpline::CubicSpline() :
    pImpl(std::make_unique<CubicSplineImpl> ())
{
}

/// Move constructor
CubicSpline::CubicSpline(CubicSpline &&spline) noexcept
{
    *this = std::move(spline);
}

/// Move assignment operator
CubicSpline& CubicSpline::operator=(CubicSpline &&spline) noexcept
{
    if (&spline == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(spline.pImpl);
    return *this;
}

/// Destructor
CubicSpline::~CubicSpline() = default;

/// Clear memory
void CubicSpline::clear() noexcept
{
    pImpl->clear();
}

/// Initialize the cubic spline for a uniform interval 
void CubicSpline::initialize(
    const int npts,
    const std::pair<double, double> xInterval,
    const double y[],
    const CubicSplineBoundaryConditionType boundaryConditionType)
{
    clear();
    // Check the inputs
    if (npts < 4)
    {
        RTSEIS_THROW_IA("npts = %d must be at least 4", npts);
    }
    if (y == nullptr){RTSEIS_THROW_RTE("%s", "y is NULL");}
    if (xInterval.first >= xInterval.second)
    {
        RTSEIS_THROW_IA("x.first = %lf must be less than x.second = %lf",
                        xInterval.first, xInterval.second); 
    }
    if (boundaryConditionType == CubicSplineBoundaryConditionType::PERIODIC)
    {
        // MKL is pretty stringent
        if (y[0] != y[npts-1])
        {
            RTSEIS_THROW_RTE("y[0] = %e != y[npts-1] = %e", y[0], y[npts-1]);
        }
    }
    // Create the pipeline
    auto ierr = pImpl->createUniformXTask(npts, xInterval, y);
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Failed to create task");}
    // Edit the pipeline to inform MKL which spline to create
    ierr = pImpl->editPipeline(boundaryConditionType);
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Failed to edit spline");}
    // Construct the task
    ierr = pImpl->constructSpline();
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Failed to construct spline");}
    pImpl->mInitialized = true;
}

/// Initialize the cubic spline for arbitrary points
void CubicSpline::initialize(
    const int npts,
    const double x[],
    const double y[],
    const CubicSplineBoundaryConditionType boundaryConditionType)
{
    clear();
    // Check the inputs
    pImpl->mBCType = boundaryConditionType;
    if (npts < 4)
    {
        RTSEIS_THROW_IA("npts = %d must be at least 4", npts);
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_RTE("%s", "y is NULL");
    }
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
    if (boundaryConditionType == CubicSplineBoundaryConditionType::PERIODIC)
    {
        // MKL is pretty stringent
        if (y[0] != y[npts-1])
        {
            RTSEIS_THROW_RTE("y[0] = %e != y[npts-1] = %e", y[0], y[npts-1]);
        }
    }
#ifdef DEBUG
    // Create the pipeline
    auto ierr = pImpl->createTask(npts, x, y);
    assert(ierr == 0);
    // Edit the pipeline to inform MKL which spline to create
    ierr = pImpl->editPipeline(boundaryConditionType);
    assert(ierr == 0);
    // Construct the task
    ierr = pImpl->constructSpline();
    assert(ierr == 0);
#else
    // Create the pipeline
    auto ierr = pImpl->createTask(npts, x, y);
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Failed to create task");}
    // Edit the pipeline to inform MKL which spline to create
    ierr = pImpl->editPipeline(boundaryConditionType);
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Failed to edit spline");}
    // Construct the task
    ierr = pImpl->constructSpline();
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Failed to construct spline");}
#endif
    pImpl->mInitialized = true;
}

/// Check if initialized
bool CubicSpline::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Get minimum x for interpolation
double CubicSpline::getMinimumX() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mXMin;
}

/// Get maximum x for interpolation
double CubicSpline::getMaximumX() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mXMax;
}

/// Interpolate at a bunch of points
void CubicSpline::interpolate(const int nq, const double xq[],
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
    // Check if this is sorted
    bool lsorted = Math::VectorMath::isSorted(nq, xq);
    // Interpolate
    int ierr = pImpl->interpolate(nq, xq, yq, lsorted); 
#ifdef DEBUG
    assert(ierr == 0);
#endif
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Interpolation failed");}
}

/// Integrate over an interval
double CubicSpline::integrate(const std::pair<double, double> &interval) const
{
    // Integral of a polynomial over interval of 0 length is 0.
    if (interval.first == interval.second){return 0;}
    // Now start performing checks
    auto xmin = getMinimumX(); // Throws on unitialized error
    auto xmax = getMaximumX();
    if (interval.first < xmin || interval.first > xmax)
    {
        RTSEIS_THROW_IA("interval.first = %lf must be in range [%lf,%lf]",
                        interval.first, xmin, xmax);
    }
    if (interval.second < xmin || interval.second > xmax)
    {
        RTSEIS_THROW_IA("interval.second = %lf must be in range [%lf,%lf]",
                        interval.second, xmin, xmax);
    }
    // Integrate 
    int ierr;
    double integral;
    if (interval.first > interval.second)
    {
        ierr = pImpl->integrate(interval.first, interval.second, &integral);
#ifdef DEBUG
        assert(ierr == 0);
#endif
    }
    else
    {
        // Reversing integration limits amounts to a sign change
        ierr = pImpl->integrate(interval.second, interval.first, &integral);
#ifdef DEBUG
        assert(ierr == 0);
#endif
        integral =-integral;
    }
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Integration failed");}
    return integral;
}

void CubicSpline::integrate(const int nIntervals,
                            const std::pair<double, double> intervals[],
                            double *integralsIn[]) const
{
    if (nIntervals < 1){return;}
    auto xmin = getMinimumX(); // Throws on unitialized error
    auto xmax = getMaximumX();
    // Check the inputs
    double *integrals = *integralsIn;
    if (intervals == nullptr || integrals == nullptr)
    {
        if (intervals == nullptr){RTSEIS_THROW_IA("%s", "intervals is NULL");}
        RTSEIS_THROW_IA("%s", "integrals is NULL");
    }
    std::vector<double> llim(nIntervals, 0);
    std::vector<double> rlim(nIntervals, 0);
    std::vector<double> swapSign(nIntervals, 0);
    double intervalMin = DBL_MAX;
    double intervalMax =-DBL_MAX;
    #pragma omp simd reduction(max:intervalMax) reduction(min:intervalMin)
    for (auto i=0; i<nIntervals; ++i)
    {
        llim[i] = intervals[i].first;
        rlim[i] = intervals[i].second;
        swapSign[i] = false;
        // Enforce left integrand as lower limit 
        if (intervals[i].first > intervals[i].second)
        {
            llim[i] = intervals[i].second;
            rlim[i] = intervals[i].first;
            swapSign[i] = true;
        }
        intervalMin = std::min(intervalMin, llim[i]);
        intervalMax = std::max(intervalMax, rlim[i]);
    }
    // Check the intervals
    if (intervalMin < xmin || intervalMax > xmax)
    {
        RTSEIS_THROW_IA("At least one interval is out of bounds [%lf,%lf]",
                        xmin, xmax);
    }
    // Integrate
    auto ierr = pImpl->integrate(nIntervals, llim.data(), rlim.data(),
                                 integrals); 
#ifdef DEBUG
    assert(ierr == 0):
#else
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Integration failed");}
#endif
    // When the integrand limits this manifests as the integral being off by
    // a sign factor.  This loop fixes that.
    //#pragma omp simd
    for (auto i=0; i<nIntervals; ++i)
    {
        if (swapSign[i]){integrals[i] =-integrals[i];}
    }
}
