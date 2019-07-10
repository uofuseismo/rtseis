#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <mkl.h>

using namespace RTSeis::Utilities::Math::Interpolation::CubicSpline;

/*!
 */
enum class BoundaryConditionType
{
    NOT_A_KNOT, /*!< The first and second segment at a curve end are the
                     same.  This is useful when there is no information on
                     the boundary conditions. */
    NATURAL,    /*!< The second derivative at the curve ends are set to 0, i.e.,
                     \f$ y''[0] = y''[ny-1] \f$. */
    CLAMPED,    /*!< The first derivative at the curve ends are set to 0, i.e.,
                     \f$  y'[0] = y'[ny-1] \f$. */
    PERIODIC    /*!< The interpolated function is assumed to be periodic.
                     In this case, the first y value must equal the last
                     y value.  The resulting boundary condition will be
                     \f$  y'[0] = y'[ny-1] \f$ and 
                     \f$ y''[0] = y''[ny-1] \f$. */
};

class CubicSpline::CubicSplineImpl
{
public:
    /// Creates a new task for when x is uniform
    int createUniformXTask(const int npts,
                           const std::pair<double, double> xPair,
                           const double y[],
                           const int MKL_INT yHint = DF_NO_HINT)
    {
         deleteTask();
         mXMin = x.first;
         mXMax = x.second;
         const double x[2] = {mXMin, mXMax};
         constexpr MKL_INT nx = 2;
         const MKL_INT ny = npts;
         constexpr MKL_INT xHint = DF_UNIFORM_PARTITION;
         mSplineCoeffs.resize(ny*splineOrder*(nx - 1));
         auto status = dfdNewTask1D(&mTask, nx, x, xHint, ny, y, yhint); 
         if (status != MKL_STATUS_OK){return -1;}
         mhaveTask = true;
         return 0;
    }
    /// Creates a new task for when x is not necessarily uniofrm 
    int createTask(const int npts,
                   const double x[],
                   const double y[],
                   const int MKL_INT yHint = DF_NO_HINT)
    {
         deleteTask();
         mXMin = x[0];
         mXMax = x[npts-1]; 
         const MKL_INT nx = npts;
         const MKL_INT ny = npts;
         constexpr MKL_INT xHint = DF_UNIFORM_PARTITION;
         mSplineCoeffs.resize(ny*splineOrder*(nx - 1));
         auto status = dfdNewTask1D(&mTask, nx, x, xHint, ny, y, yhint);
         if (status != MKL_STATUS_OK){return -1;}
         mhaveTask = true;
         return 0;
    }
    /// Sets the spline type and boundary conditions on the task
    int editPipeline(const BoundaryConditionType boundaryConditionType)
    {
        // Figure out the boundary conditions first
        const double bcArray[2] = {0, 0};
        const double *bcs = nullptr;
        const double *ics = nullptr; 
        if (boundaryConditionType == BoundaryConditionType::NATURAL)
        {
            pImpl->mSplineBC = DF_BC_FREE_END;
            pImpl->mSplineIC = DF_NO_IC;
        }
        else if (boundaryConditionType == BoundaryConditiontype::CLAMPED)
        {
            bcs = bcArray;
            pImpl->mSplineBC = DF_BC_1ST_LEFT_DER | DF_BC_2ND_RIGHT_DER;
            pImpl->mSplineIC = DF_NO_IC;
        }
        else if (boundaryConditionType == BoundaryConditionType::NOT_A_KNOT)
        {
            pImpl->mSplineBC = DF_BC_NOT_A_KNOT;
            pImpl->mSplineIC = DF_NO_IC;
        }
        else if (boundaryConditionType == BoundaryConditionType::PERIODIC)
        {
            pImpl->mSplineBC = DF_BC_PERIODIC;
            pImpl->mSplineIC = DF_NO_IC;
        }
        constexpr MKL_INT scoeffHint = DF_HINT_HINT;
        auto status = dfdEditPPSpline1D(mTask,
                                        pImpl->mSplineOrder,
                                        pImpl->mSplineType,
                                        pImpl->mSplineBC,
                                        bcs,
                                        pImpl->mSplineIC,
                                        ics,
                                        pImpl->mSplineCoeffs.data(),
                                        scoeffHint);
        if (status != MKL_SUCCESS)
        {
            clear();
            return -1;
        }
        return 0;
    }
    /// Makes the spline
    int constructSpline()
    {
        auto status = dfdConstruct1D(pImpl->mTask,
                                     pImpl->mStatusFormat,
                                     pImpl->mMethod);
        if (status != MKL_SUCCESS)
        {
            clear();
            return -1; 
        }
        return 0;
    }
    /// Deletes the task
    void deleteTask() noexcept
    {
        if (mHaveTask){dfDeleteTask(&mTask);}
        mHaveTask = false;
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
    bool mHaveTask = false;
    bool mInitialized = false;
};

CubicSpline::CubicSpline() :
    pImpl(std::make_unique<CubicSplineImpl> ())
{
}

CubicSpline::~CubicSpline()
{
    clear();
}

void CubicSpline::clear() noexcept
{
    pImpl->deleteTask();
    pImpl.mSplineCoeffs.clear();
    pImpl->mXMin =-DBL_MAX;
    pImpl->mXMax = DBL_MAX; 
    pImpl->mSplineBC = DF_BC_FREE_END;
    pImpl->mSplineIC = DF_NO_IC;
    pImpl->mHaveTask = false;
    pImpl->mInitialized = false;
}

void CubicSpline::initialize(
    const int npts,
    const std::pair<double, double> x,
    const double y[],
    const BoundaryConditionType boundaryConditionType)
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
    // Create the data fitting task
    pImpl->creatNewUniformXTask(npts, x, y);
    // Edit the pipeline to inform MKL which spline to create
    pImpl->editPipeline(boundaryConditionType);
    // Construct the task
    pImpl->constructSpline();
}

void CubicSpline::initialize(
    const int npts,
    const double xIn[],
    const double yIn[],
    const BoundaryConditionType boundaryConditionType)
{
    clear();
    // Check the inputs
    if (npts < 4)
    {
        RTSEIS_THROW_IA("npts = %d must be at least 4", npts);
    }
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_RTE("%s", "y is NULL");
    }
    // Assert that the data is sorted
    auto lsorted = RTSeis::Utilities::Math::VectorMath::isSorted(npts, x);
    const double *x = xIn;
    const double *y = yIn;
    if (!lsorted)
    {
        // Sort
        std::vector<int> permutation(npts);
        ippsSortIndexAscend_64f_I(x, permuation.data(), npts);
        std::vector<double> xSorted(npts);
        std::vector<double> ySorted(npts);
        #pragma omp simd
        for (auto i=0; i<npts; ++i)
        {
           xSorted[i] = x[permutation[i]];
           ySorted[i] = y[permutation[i]];
        }
        if (xSorted[0] >= xSorted[npts-1])
        {
           RTSEIS_THROW_IA("min(x) = %lf must be less than max(x) = %lf",
                           xSorted[0], xSorted[npts-1]);
        }
        // Create the data fitting task
        pImpl->creatNewTask(npts, x, y);
    }
    else
    {

    }

    // Edit the pipeline to inform MKL which spline to create
    pImpl->editPipeline(boundaryConditionType);
    // Construct the task
    pImpl->constructSpline();
}

void CubicSpline::integrate(std::vector<std::pair<double,double>> intervals)
{

}

double CubicSpline::integrate(const double x1, const double x2)
{

}
