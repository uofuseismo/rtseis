#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>
#include <mkl.h>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/math/cubicSpline.hpp"

using namespace RTSeis::Utilities::Math::Interpolation;

class CubicSpline::CubicSplineImpl
{
public:
    /// Destructor
    ~CubicSplineImpl()
    {
        clear();
    }
    /// Creates a new task for when x is uniform
    int createUniformXTask(const int npts,
                           const std::pair<double, double> xInterval,
                           const double y[],
                           const MKL_INT yHint = DF_NO_HINT)
    {
        deleteTask();
        mXMin = xInterval.first;
        mXMax = xInterval.second;
        const double x[2] = {mXMin, mXMax};
        constexpr MKL_INT nx = 2;
        const MKL_INT ny = npts;
        constexpr MKL_INT xHint = DF_UNIFORM_PARTITION;
        mSplineCoeffs.resize(ny*mSplineOrder*(nx - 1));
        auto status = dfdNewTask1D(&mTask, nx, x, xHint, ny, y, yHint); 
        if (status != DF_STATUS_OK){return -1;}
        mHaveTask = true;
        return 0;
    }
    /// Creates a new task for when x is not necessarily uniofrm 
    int createTask(const int npts,
                   const double x[],
                   const double y[],
                   const MKL_INT yHint = DF_NO_HINT)
    {
        deleteTask();
        mXMin = x[0];
        mXMax = x[npts-1]; 
        const MKL_INT nx = npts;
        const MKL_INT ny = npts;
        constexpr MKL_INT xHint = DF_NON_UNIFORM_PARTITION;
        mSplineCoeffs.resize(ny*mSplineOrder*(nx - 1));
        auto status = dfdNewTask1D(&mTask, nx, x, xHint, ny, y, yHint);
        if (status != DF_STATUS_OK){return -1;}
        mHaveTask = true;
        return 0;
    }
    /// Sets the spline type and boundary conditions on the task
    int editPipeline(const CubicSplineBoundaryConditionType boundaryCondition)
    {
        // Figure out the boundary conditions first
        const double bcArray[2] = {0, 0};
        const double *bcs = nullptr;
        const double *ics = nullptr; 
        if (boundaryCondition == CubicSplineBoundaryConditionType::NATURAL)
        {
            mSplineBC = DF_BC_FREE_END;
            mSplineIC = DF_NO_IC;
        }
        else if (boundaryCondition == CubicSplineBoundaryConditionType::CLAMPED)
        {
            bcs = bcArray;
            mSplineBC = DF_BC_1ST_LEFT_DER | DF_BC_2ND_RIGHT_DER;
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
        auto status = dfdConstruct1D(mTask, mSplineFormat, mMethod);
        if (status != DF_STATUS_OK)
        {
            clear();
            return -1; 
        }
        return 0;
    }
    /// Searches the cells for the interpolant points
    std::vector<MKL_INT> searchCells(const int nq, const double xq[],
                                     const bool lsorted = false)
    {
        const MKL_INT nsite = nq;
        std::vector<MKL_INT> cell(nsite);
        if (nq < 1){return cell;}
        MKL_INT sortedHint = DF_SORTED_DATA;
        if (!lsorted){sortedHint = DF_NO_HINT;}
        dfdSearchCells1D(mTask, DF_METHOD_STD, nsite, xq,
                         sortedHint, DF_NO_APRIORI_INFO,
                         cell.data());
        return cell;
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
                    double yq[])
    {
        const MKL_INT nsite = nq;
        std::vector<MKL_INT> cell = searchCells(nq, xq);
        dfdSearchCells1D(mTask, DF_METHOD_STD, nsite, xq,
                         DF_SORTED_DATA, DF_NO_APRIORI_INFO,
                         cell.data()); 
        constexpr MKL_INT nOrder = 1;  // Length of dorder
        const MKL_INT dOrder[1] = {0}; // Derivative order
        auto status = dfdInterpolate1D(mTask, DF_INTERP, DF_METHOD_PP,
                                       nsite, xq, DF_NON_UNIFORM_PARTITION,
                                       nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yq,
                                       DF_MATRIX_STORAGE_ROWS, cell.data());
        if (status != DF_STATUS_OK){return -1;}
        return 0;
    }
    /// Integrates over the interval
    int integrate(const double lLim, const double rLim, double *integral)
    {
        constexpr int nlim = 1;
        double llim[1] = {lLim};
        double rlim[1] = {rLim};
        auto status = dfdIntegrate1D(mTask, DF_METHOD_PP, nlim,
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
        mXMin =-DBL_MAX;
        mXMax = DBL_MAX;
        mSplineBC = DF_BC_FREE_END;
        mSplineIC = DF_NO_IC;
        mBCType = CubicSplineBoundaryConditionType::NATURAL;
        mHaveTask = false;
        mInitialized = false;
    }
///private:
    DFTaskPtr mTask;
    std::vector<double> mSplineCoeffs;
    double mXMin =-DBL_MAX;
    double mXMax = DBL_MAX;    
    /// Default to natural cubic spline
    MKL_INT mSplineBC = DF_BC_FREE_END;
    MKL_INT mSplineIC = DF_NO_IC;
    const MKL_INT mSplineType = DF_PP_NATURAL;  // Deal with type in b.c.'s
    const MKL_INT mSplineOrder = DF_PP_CUBIC;   // This is for cubic splines
    const MKL_INT mSplineFormat = DF_PP_SPLINE; // Only supported option
    const MKL_INT mMethod = DF_METHOD_STD;      // Only supported option
    CubicSplineBoundaryConditionType mBCType
        = CubicSplineBoundaryConditionType::NATURAL;
    bool mHaveTask = false;
    bool mInitialized = false;
};

CubicSpline::CubicSpline() :
    pImpl(std::make_unique<CubicSplineImpl> ())
{
}

CubicSpline::CubicSpline(CubicSpline &&spline) noexcept
{
    *this = std::move(spline);
}

CubicSpline& CubicSpline::operator=(CubicSpline &&spline) noexcept
{
    if (&spline == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(spline.pImpl);
    return *this;
}

CubicSpline::~CubicSpline() = default;

void CubicSpline::clear() noexcept
{
    pImpl->clear();
}

/*
void CubicSpline::initialize(
    const int npts,
    const std::pair<double, double> xInterval,
    const double y[],
    const CubicSplineBoundaryConditionType boundaryCondition)
{
    clear();
    // Check the inputs
    if (npts < 4)
    {
        RTSEIS_THROW_IA("npts = %d must be at least 4", npts);
    }
    if (y == nullptr){RTSEIS_THROW_RTE("%s", "y is NULL");}
    if (x.first >= x.second)
    {
        RTSEIS_THROW_IA("x.first = %lf must be less than x.second = %lf",
                        x.first, x.second); 
    }
}
*/

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
        if (x[i] >= x[i+1]){isInvalid = isInvalid + 1;}
    }
    if (isInvalid > 0)
    {
        RTSEIS_THROW_RTE("%s", "At least one x[i+1] <= x[i]\n");
    }
    if (boundaryConditionType == CubicSplineBoundaryConditionType::PERIODIC)
    {
        if (y[0] != y[npts-1])
        {
            RTSEIS_THROW_RTE("y[0] = %e != y[npts-1] = %e", y[0], y[npts-1]);
        }
    }
#ifdef DEBUG
    // Create the pipeline
    auto ierr = pImpl->createTask(npts, x, y, DF_NO_HINT);
    assert(ierr == 0);
    // Edit the pipeline to inform MKL which spline to create
    ierr = pImpl->editPipeline(boundaryConditionType);
    assert(ierr == 0);
    // Construct the task
    ierr = pImpl->constructSpline();
    assert(ierr == 0);
#else
    // Create the pipeline
    auto ierr = pImpl->createTask(npts, x, y, DF_NO_HINT);
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

bool CubicSpline::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

double CubicSpline::getMinimumX() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mXMin;
}

double CubicSpline::getMaximumX() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mXMax;
}

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
    // Interpolate
    int ierr = pImpl->interpolate(nq, xq, yq); 
#ifdef DEBUG
    assert(ierr == 0);
#endif
    if (ierr != 0){RTSEIS_THROW_RTE("%s", "Interpolation failed");}
}

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


