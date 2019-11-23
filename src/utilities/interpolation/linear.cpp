#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <cmath>
#include <algorithm>
#include <exception>
#include <assert>
#include <stdexcept>
#include <mkl.h>
#include <ipps.h>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/interpolation/linear.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"
#include "rtseis/log.h"

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
        if (mHaveTask){dfDeleteTask(&task);}
        if (mSplineCoeffs){MKL_free(mSplineCoeffs);}
        mXMin = std::numeric_limits<double>::max();
        mXMax = std::numeric_limits<double>::min();
        mHaveTask = false;
        mInitialized = false;
    }
    /// Builds the interpolator
    void initialize(const int npts,
                    const std::pair<double, double> x,
                    const double v[])
    {
        clear();
        constexpr int ny = 1; // Function to interpolate is scalar
        nwork = std::max(8, ny *mSplineOrder*(npts - 1));
        double xpts[2] = {x.first, x.second};
        auto status = dfdNewTask1D(&task, nx, xpts, DF_UNIFORM_PARTITION,
                                   ny, v, DF_NO_HINT);
#ifdef DEBUG
        assert(status == DF_STATUS_OK);
#else
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Task creation failed\n");
        }
#endif
        mHaveTask = true;
        // Set spline parameters in the data fitting task
        double *bc = nullptr;
        double *ic = nullptr;
        mSplineCoeffs = static_cast<double *>
                        (MKL_calloc(static_cast<size_t>(nwork), sizeof(double));
        status = dfdEditPPSpline1D(task, splineOrder, splineType, splineBC,
                                   bc, splineIC, ic, mSplineCoeffs, DF_NO_HINT);
#ifdef DEBUG
        assert(status == DF_STATUS_OK);
#else
        if (status != DF_STATUS_OK)
        {
            throw std::runtime_error("Spline editing failed\n");
        }
#endif
        // Construct the spline
        mXMin = x.first;
        mXMax = x.second;
    }
 
    DFTaskPtr mTask;
    const MKL_INT mSplineOrder = DF_PP_LINEAR;
    const MKL_INT mSplineType  = DF_PP_LINEAR;
    const MKL_INT mSplineBC    = DF_NO_BC;
    const MKL_INT mSplineIC    = DF_NO_IC;
    double mXMin = std::numeric_limits<double>::max();
    double mXMax = std::numeric_limits<double>::min();
    double *mSplineCoeffs = nullptr;
    bool mHaveTask = false;
    bool mInitialized = false;
};

/// Constructor
Linear::Linear() :
    pImpl(std::make_unique<LinearImpl> ())
{
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
