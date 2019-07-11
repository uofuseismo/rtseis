#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>
#include <mkl.h>
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
        constexpr MKL_INT scoeffHint = DF_NO_HINT;
        auto status = dfdEditPPSpline1D(mTask,
                                        mSplineOrder,
                                        mSplineType,
                                        mSplineBC,
                                        bcs,
                                        mSplineIC,
                                        ics,
                                        mSplineCoeffs.data(),
                                        scoeffHint);
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
    /// Interpolates over the uniform interval
    int interpolate(const int nq,
                    const std::pair<double, double> xInterval,
                    double yq[])
    {
        double x[2] = {xInterval.first, xInterval.second};
        constexpr MKL_INT nOrder = 1;  // Length of dorder
        const MKL_INT dOrder[1] = {0}; // Order of derivatives
        auto status = dfdInterpolate1D(mTask, DF_INTERP, DF_METHOD_PP,
                                       nq, x, DF_UNIFORM_PARTITION,
                                       nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yq,
                                       DF_NO_HINT, NULL);
        if (status != DF_STATUS_OK){return -1;}
        return 0; 
    }
    /// Interpolates at select abscissa
    int interpolate(const int nq,
                    const double x[],
                    double yq[])
    {
        constexpr MKL_INT nOrder = 1;  // Length of dorder
        const MKL_INT dOrder[1] = {0}; // Order of derivatives
        auto status = dfdInterpolate1D(mTask, DF_INTERP, DF_METHOD_PP,
                                       nq, x, DF_NON_UNIFORM_PARTITION,
                                       nOrder, dOrder,
                                       DF_NO_APRIORI_INFO, yq,
                                       DF_NO_HINT, NULL);
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
    pImpl->createTask(npts, x, y, DF_NO_HINT);
    // Edit the pipeline to inform MKL which spline to create
    pImpl->editPipeline(boundaryConditionType);
    // Construct the task
    pImpl->constructSpline();
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
/*
void CubicSpline::integrate(std::vector<std::pair<double,double>> intervals)
{

}

double CubicSpline::integrate(const double x1, const double x2)
{

}
*/
