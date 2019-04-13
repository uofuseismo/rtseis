#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <mkl.h>
#include <ipps.h>
#define RTSEIS_LOGGING 1
#include "rtseis/utilities/math/interpolate.hpp"
#include "rtseis/utilities/math/vectorMath.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utilities::Math::Interpolate;
using namespace RTSeis::Utilities::Math;

void Interpolate::interpft(
    const std::vector<double> &x,
    const int npnew,
    std::vector<double> &yint)
{
    yint.resize(0);
    if (npnew < 1)
    {
        throw std::invalid_argument("No points at which to intepolate");
    }
    int npts = static_cast<int> (x.size());
    if (npts < 1)
    {
        throw std::invalid_argument("x is empty");
    }
    // Straight copy case
    yint.resize(npnew);
    if (npts == npnew)
    {
        ippsCopy_64f(x.data(), yint.data(), npnew);
        return;
    }
    // Figure out the size of the forward and inverse transforms
    IppStatus status;
    int bufferSizeI, bufferSizeF, specSizeF, specSizeI, sizeInitF, sizeInitI;
    status = ippsDFTGetSize_R_64f(npts,
                                  IPP_FFT_DIV_INV_BY_N,
                                  ippAlgHintNone,
                                  &specSizeF,
                                  &sizeInitF,
                                  &bufferSizeF);
    if (status != ippStsNoErr)
    {
        throw std::invalid_argument("Forward transform inquiry failed");
    }
    status = ippsDFTGetSize_R_64f(npnew,
                                  IPP_FFT_DIV_INV_BY_N,
                                  ippAlgHintNone,
                                  &specSizeI,
                                  &sizeInitI,
                                  &bufferSizeI);
    if (status != ippStsNoErr)
    {   
        throw std::invalid_argument("Inverse transform inquiry failed");
    }
    // Initialize the forward transform
    auto *pDFTForwardSpec = (IppsDFTSpec_R_64f *) ippsMalloc_8u(specSizeF);
    Ipp8u *pDFTInitBuf = ippsMalloc_8u(std::max(sizeInitF, sizeInitI));
    status = ippsDFTInit_R_64f(npts, 
                               IPP_FFT_DIV_FWD_BY_N,
                               ippAlgHintNone,
                               pDFTForwardSpec,
                               pDFTInitBuf);
    if (status != ippStsNoErr)
    {
        ippsFree(pDFTForwardSpec);
        ippsFree(pDFTInitBuf);
        throw std::invalid_argument("Forward transform init failed");
    }
    // Initialize the inverse transform
    auto *pDFTInverseSpec = (IppsDFTSpec_R_64f *) ippsMalloc_8u(specSizeI);
    status = ippsDFTInit_R_64f(npnew,
                               IPP_FFT_DIV_FWD_BY_N,
                               ippAlgHintNone,
                               pDFTInverseSpec,
                               pDFTInitBuf);
    if (pDFTInitBuf){ippsFree(pDFTInitBuf);}
    if (status != ippStsNoErr)
    {
        ippsFree(pDFTForwardSpec);
        ippsFree(pDFTInverseSpec);
        ippsFree(pDFTInitBuf);
        throw std::invalid_argument("Inverse transform init failed");
    }
    // Set the workspace
    Ipp8u *pBuf = ippsMalloc_8u(std::max(bufferSizeF, bufferSizeI));
    int maxDFTLen = std::max(npts/2+1, npnew/2+1);
    Ipp64f *pDst = ippsMalloc_64f(2*maxDFTLen); // Hold real and complex
    ippsZero_64f(pDst, 2*maxDFTLen); // Pre-zero-pad in frequency domain
    // Forward transform
    ippsDFTFwd_RToCCS_64f(x.data(), pDst, pDFTForwardSpec, pBuf);
    // It's pre-zero padded so inverse transform
    ippsDFTInv_CCSToR_64f(pDst, yint.data(), pDFTInverseSpec, pBuf);
    // Clean up
    ippsFree(pBuf);
    ippsFree(pDst);
    ippsFree(pDFTForwardSpec);
    ippsFree(pDFTInverseSpec);
}

class Interp1D::Interp1DImpl
{
    public:
        /// Default constructor
        Interp1DImpl()
        {
            std::memset(&task, 0, sizeof(DFTaskPtr));
        }
        /// Can enable copy construction once task is re-initialized
        Interp1DImpl(const Interp1DImpl &interp1d) = delete;
        /// Destructor
        ~Interp1DImpl()
        {
            clear();
        }
        /// Copy assignment operator
        Interp1DImpl &operator=(const Interp1DImpl &interp1d)
        {
            if (&interp1d == this){return *this;}
            if (!interp1d.linit){return *this;}
            initialize(interp1d.xin, interp1d.vin, interp1d.method);
            return *this; 
        }
        /// Releases the memory on the interpolator
        void clear()
        {
            xin.clear();
            vin.clear();
            if (lhaveTask){dfDeleteTask(&task);}
            delete[] scoeff;
            scoeff = nullptr;
            xmin = 0;
            xmax = 0;
            splineOrder = DF_PP_STD;
            splineType  = DF_CR_STEPWISE_CONST_INTERPOLANT;
            splineBC    = DF_NO_BC;
            splineIC    = DF_NO_IC;
            nwork = 0;
            method = Interp1D::Method::NEAREST;
            //luniform = false;
            lhaveTask = false;
            linit = false;
        }
        /// Builds the interpolator
        int initialize(const std::vector<double> &x,
                       const std::vector<double> &v,
                       const Interp1D::Method)
        {
            clear();
            int nx = static_cast<int> (x.size());
            setSplineInfo(method);
            // Initialize space
            constexpr int ny = 1; // Function to interpolate is scalar
            nwork = std::max(8, ny*splineOrder*(nx - 1));
            // Create the data fitting task
            //DFTaskPtr task;
            MKL_INT status = dfdNewTask1D(&task, nx, x.data(), DF_NO_HINT,
                                          ny, v.data(), DF_NO_HINT);
            lhaveTask = true;
            if (status != DF_STATUS_OK)
            {
                RTSEIS_ERRMSG("%s", "Failed to create task");
                clear();
                return -1;
            }
            // Set spline parameters in the data fitting task
            double *bc = nullptr;
            double *ic = nullptr;
            scoeff = new double[nwork]; // Holds
            status = dfdEditPPSpline1D(task, splineOrder, splineType, splineBC,
                                       bc, splineIC, ic, scoeff, DF_NO_HINT);
            if (status != DF_STATUS_OK)
            {
                clear();
                return -1;
            }
            // Construct the spline
            if (method == Interp1D::Method::NEAREST)
            {
                status = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
                if (status != DF_STATUS_OK)
                {
                    RTSEIS_ERRMSG("%s", "Failed to compute spline");
                    clear();
                    return -1;
                }
            }
            // Get mins and maxes
            ippsMin_64f(x.data(), nx, &xmin);
            ippsMax_64f(x.data(), nx, &xmax);
            // Save data
            xin = x;
            vin = v;
            linit = true;
            return 0;
        }
        /// Builds the interpolator for uniform points 
        int initialize(const int npts,
                       const std::pair<double, double> x,
                       const std::vector<double> &v, 
                       const Interp1D::Method)
        {
            clear();
            int nx = npts;
            setSplineInfo(method);
            // Initialize space
            constexpr int ny = 1; // Function to interpolate is scalar
            nwork = std::max(8, ny*splineOrder*(nx - 1));
            double xpts[2] = {x.first, x.second};
            // Create the data fitting task
            //DFTaskPtr task;
            MKL_INT status = dfdNewTask1D(&task, nx, xpts, DF_UNIFORM_PARTITION,
                                          ny, v.data(), DF_NO_HINT);
            lhaveTask = true;
            if (status != DF_STATUS_OK)
            {
                RTSEIS_ERRMSG("%s", "Failed to create task");
                clear();
                return -1; 
            }
            // Set spline parameters in the data fitting task
            double *bc = nullptr;
            double *ic = nullptr;
            scoeff = new double[nwork]; // Holds
            status = dfdEditPPSpline1D(task, splineOrder, splineType, splineBC,
                                       bc, splineIC, ic, scoeff, DF_NO_HINT);
            if (status != DF_STATUS_OK)
            {
                clear();
                return -1; 
            }
            // Construct the spline
            if (method == Interp1D::Method::NEAREST)
            {
                status = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
                if (status != DF_STATUS_OK)
                {
                    RTSEIS_ERRMSG("%s", "Failed to compute spline");
                    clear();
                    return -1; 
                }
            }
            // Get mins and maxes
            xmin = x.first;
            xmax = x.second;
            // Save data
            //luniform = true;
            double dx = 0;
            if (nx > 1){dx = (xmax - xmin)/static_cast<double> (nx - 1);}
            xin.resize(nx); 
            #pragma omp simd
            for (int i=0; i<nx; i++)
            {
                xin[i] = xmin + static_cast<double> (i)*dx;
            }
            vin = v;
            linit = true;
            return 0;
        }
        /// Sets the spline information
        void setSplineInfo(const Interp1D::Method)
        {
            if (method == Interp1D::Method::NEAREST)
            {
                splineOrder = DF_PP_STD;
                splineType  = DF_CR_STEPWISE_CONST_INTERPOLANT;
                splineBC    = DF_NO_BC;
                splineIC    = DF_NO_IC;
            }
            else if (method == Interp1D::Method::LINEAR)
            {
                splineOrder = DF_PP_LINEAR;
                splineType  = DF_PP_LINEAR;
                splineBC    = DF_NO_BC;
                splineIC    = DF_NO_IC;
            }
            else if (method == Interp1D::Method::CSPLINE_NATURAL)
            {
                splineOrder = DF_PP_CUBIC;
                splineType  = DF_PP_NATURAL;
                splineBC    = DF_BC_FREE_END;
                splineIC    = DF_NO_IC;
            }
            else if (method == Interp1D::Method::CSPLINE_PERIODIC)
            {
                splineOrder = DF_PP_CUBIC;
                splineType  = DF_PP_NATURAL;
                splineBC    = DF_BC_PERIODIC;
                splineIC    = DF_NO_IC;
            }
            else if (method == Interp1D::Method::AKIMA)
            {
                splineOrder = DF_PP_CUBIC;
                splineType  = DF_PP_AKIMA;
                splineBC    = DF_BC_FREE_END;
                splineIC    = DF_NO_IC;
            }
            else //if (method == Interp1D::Method::AKIMA_PERIODIC)
            {
                splineOrder = DF_PP_CUBIC;
                splineType  = DF_PP_AKIMA;
                splineBC    = DF_BC_PERIODIC;
                splineIC    = DF_NO_IC;
            }
        }
        /// Interpolate
        int interpolate(const std::vector<double> &xq,
                        std::vector<double> &vq,
                        const bool lsorted)
        {
           int nx = static_cast<int> (xq.size());
           vq.resize(nx);
           if (nx == 0){return 0;}
           // Check mins/maxes
           double xqMin, xqMax;
           ippsMin_64f(xq.data(), nx, &xqMin);
           ippsMax_64f(xq.data(), nx, &xqMax);
           const double *xqData = nullptr;
           if (xqMin >= xmin && xqMax <= xmax)
           {
               xqData = xq.data(); 
           }   
           else
           {
               // Need to threshold
               std::vector<double> xqThresh(nx);
               #pragma omp simd
               for (int i=0; i<nx; i++)
               {
                   xqThresh[i] = std::min(xmax, std::max(xmin, xq[i]));
               }
               xqData = xqThresh.data(); 
           }
           // Compute cells of interpolant points
           MKL_INT status;
           auto *cell = static_cast<MKL_INT *>
                        (mkl_malloc(xq.size()*sizeof(MKL_INT), 64));
           auto nsite = static_cast<MKL_INT> (xq.size());
           if (lsorted)
           {
               status = dfdSearchCells1D(task, DF_METHOD_STD, nsite, xqData,
                                         DF_SORTED_DATA, DF_NO_APRIORI_INFO,
                                         cell);
           }
           else
           {
               status = dfdSearchCells1D(task, DF_METHOD_STD, nsite, xqData,
                                         DF_NO_HINT, DF_NO_APRIORI_INFO,
                                          cell);
           }
           if (status != DF_STATUS_OK)
           {
               delete[] cell;
               RTSEIS_ERRMSG("%s", "Failed searching cells");
               return -1;
           }
           // Interpolate
           int ierr = 0;
           if (method != Interp1D::Method::NEAREST)
           {
               constexpr MKL_INT norder = 1;
               const MKL_INT dorder[1] = {0}; // Order of derivatives
               if (lsorted)
               {
                   status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
                                             nsite, xqData,
                                             DF_SORTED_DATA, norder, dorder,
                                             DF_NO_APRIORI_INFO, vq.data(),
                                             DF_MATRIX_STORAGE_ROWS, cell);
               }
               else
               {
                   status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
                                             nsite, xqData,
                                             DF_NO_HINT, norder, dorder,
                                             DF_NO_APRIORI_INFO, vq.data(),
                                             DF_MATRIX_STORAGE_ROWS, cell);
               }
               if (status != 0)
               {
                   RTSEIS_ERRMSG("%s", "Interpolation failed");
                   ierr = 1; 
               }
           }
           else
           {
               for (MKL_INT i=0; i<nsite; ++i)
               {
                   MKL_INT indx = cell[i];
                   vq[i] = vin[indx];
                   if (std::abs(xqData[i] - xin[indx]) >
                       std::abs(xqData[i] - xin[indx-1]))
                   {
                       vq[i] = vin[indx-1];
                   }
               }
           }
           delete[] cell;
           return ierr;
        }
        /// Checks if the class is initialized
        bool isInitialized() const
        {
            return linit;
        }
        /// Gets the interolation method
        Interp1D::Method getMethod() const
        {
            return method;
        }
 
    private:
        DFTaskPtr task;
        std::vector<double> xin;
        std::vector<double> vin;
        double *scoeff = nullptr; // has dimension [nwork]
        double xmin = 0;
        double xmax = 0;
        MKL_INT splineOrder = DF_PP_STD;
        MKL_INT splineType  = DF_CR_STEPWISE_CONST_INTERPOLANT;
        MKL_INT splineBC    = DF_NO_BC;
        MKL_INT splineIC    = DF_NO_IC;
        int nwork = 0;
        Interp1D::Method method = Interp1D::Method::NEAREST;
        //bool luniform = false;
        bool lhaveTask = false;
        bool linit = false;
};

Interp1D::Interp1D() :
    pImpl(std::unique_ptr<Interp1DImpl>())
{
}

Interp1D::~Interp1D() = default;

void Interp1D::initialize(const std::vector<double> &x,
                          const std::vector<double> &v,
                          const Interp1D::Method method)
{
    pImpl->clear();
    if (x.size() != v.size())
    {
        throw std::invalid_argument("x.size() != v.size()");
    }
    if (method == Method::NEAREST)
    {
        if (x.empty())
        {
            throw std::invalid_argument("At least 1 point for nearest");
        }
    }
    else if (method == Method::LINEAR)
    {
        if (x.size() < 3)
        {
            throw std::invalid_argument("At least 2 points for linear");
        }
    }
    else if (method == Method::CSPLINE_NATURAL || 
             method == Method::AKIMA)
    {
        if (x.size() < 4)
        {
            throw std::invalid_argument("At least 3 points for cspline/akima");
        }
    }
    else if (method == Method::CSPLINE_PERIODIC ||
             method == Method::AKIMA_PERIODIC)
    {
        if (x.size() < 4)
        {
            throw std::invalid_argument("At least 3 points for cspline/akima");
        }
        if (v[0] != v[v.size()-1])
        {
            throw std::invalid_argument("v[0] != v[v.size()-1]");
        }
    } 
    int ierr = pImpl->initialize(x, v, method);
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Initialization failed");
        throw std::invalid_argument("Initialization failed");
    }
}

void Interp1D::initialize(const int npts,
                          const std::pair<double,double> x,
                          const std::vector<double> &v, 
                          const Interp1D::Method method)
{
    pImpl->clear();
    if (static_cast<size_t> (npts) != v.size())
    {   
        throw std::invalid_argument("x.size() != v.size()");
    }   
    if (x.first >= x.second)
    {
        throw std::invalid_argument("x.first >= x.second");
    }
    if (method == Method::NEAREST)
    {
        if (npts < 1)
        {
            throw std::invalid_argument("At least 1 point for nearest");
        }
    }   
    else if (method == Method::LINEAR)
    {   
        if (npts < 3)
        {
            throw std::invalid_argument("At least 2 points for linear");
        }
    }   
    else if (method == Method::CSPLINE_NATURAL ||  
             method == Method::AKIMA)
    {   
        if (npts < 4)
        {
            throw std::invalid_argument("At least 3 points for cspline/akima");
        }
    }   
    else if (method == Method::CSPLINE_PERIODIC ||
             method == Method::AKIMA_PERIODIC)
    {   
        if (npts < 4)
        {
            throw std::invalid_argument("At least 3 points for cspline/akima");
        }
        if (v[0] != v[v.size()-1])
        {
            throw std::invalid_argument("v[0] != v[v.size()-1]");
        }
    }   
    int ierr = pImpl->initialize(npts, x, v, method);
    if (ierr != 0)
    {   
        RTSEIS_ERRMSG("%s", "Initialization failed");
        throw std::invalid_argument("Initialization failed");
    }   
}



void Interp1D::apply(const std::vector<double> &xq,
                     std::vector<double> &vq,
                     const bool lxqSorted=false)
{
    if (!pImpl->isInitialized())
    {
        throw std::invalid_argument("Interp1D not yet inititalized");
    }
    if (xq.size() != vq.size())
    {
        throw std::invalid_argument("xq.size() != vq.size()");
    }
    bool lsorted = false;
    if (!lxqSorted){lsorted = VectorMath::isSorted(xq);}
    int ierr = pImpl->interpolate(xq, vq, lsorted); 
    if (ierr != 0)
    {
        throw std::invalid_argument("Internal failure");
    }
}

/*
void RTSeis::Utilities::Math::Interpolate::interp1d(
    const std::vector<double> &x,
    const std::vector<double> &v,
    const std::vector<double> &xq,
    std::vector<double> &vq,
    const Interpolate::Method method, 
    const bool lxSorted,
    const bool lxqSorted)
        if (nx < 4)
        {
            throw std::invalid_argument("Splines require at least 3 points");
        }
    }
    else if (method == Interpolate::Method::CSPLINE_PERIODIC)
    {
        splineOrder = DF_PP_CUBIC;
        splineType  = DF_PP_NATURAL;
        splineBC    = DF_BC_PERIODIC;
        splineIC    = DF_NO_IC;
        if (nx < 4)
        {
            throw std::invalid_argument("Splines require at least 3 points");
        } 
        if (std::abs(v[0] - v[nx-1]) > 1.e-14)
        {
            throw std::invalid_argument("v[0] != v[n-1]");
        }
    }
    else if (method == Interpolate::Method::AKIMA)
    {
        splineOrder = DF_PP_CUBIC;
        splineType  = DF_PP_AKIMA;
        splineBC    = DF_BC_FREE_END;
        splineIC    = DF_NO_IC;
        if (nx < 4)
        {
            throw std::invalid_argument("Splines require at least 3 points");
        }
    }
    else if (method == Interpolate::Method::AKIMA_PERIODIC)
    {
        splineOrder = DF_PP_CUBIC;
        splineType  = DF_PP_AKIMA;
        splineBC    = DF_BC_PERIODIC;
        splineIC    = DF_NO_IC;
        if (nx < 4)
        {   
            throw std::invalid_argument("Splines require at least 3 points");
        }
        if (std::abs(v[0] - v[nx-1]) > 1.e-14)
        {
            throw std::invalid_argument("v[0] != v[n-1]");
        }
    }
    else
    {
        throw std::invalid_argument("Unsupported interpolation type"); 
    }
    // Check the interpolation points are sorted
    bool lintSorted = true;
    if (!lxqSorted){lintSorted = VectorMath::isSorted(xq);}
    // Initialize space
    constexpr int ny = 1; // Function to interpolate is scalar
    int nwork = std::max(8, ny*splineOrder*(nx - 1));
    // Create the data fitting task
    DFTaskPtr task;
    MKL_INT status = dfdNewTask1D(&task, nx, x.data(), DF_NO_HINT,
                                  ny, v.data(), DF_NO_HINT);
    if (status != DF_STATUS_OK)
    {
        dfDeleteTask(&task);
        std::invalid_argument("Invalid inputs to dfdNewTask1D");
    }
    // Set spline parameters in the data fitting task
    double *bc = NULL;
    double *ic = NULL;
    double *scoeff = new double[nwork]; // Holds 
    status = dfdEditPPSpline1D(task, splineOrder, splineType, splineBC, bc,
                               splineIC, ic, scoeff, DF_NO_HINT);
    if (status != DF_STATUS_OK)
    {
        dfDeleteTask(&task);
        delete[] scoeff;
        std::invalid_argument("Invalid spline parameters");
    }
    // Construct the spline
    if (method == Interpolate::Method::NEAREST)
    {

    }
*/
/*
    // Create the data fitting task
    DFTaskPtr task;
    status = dfdNewTask1D(&task, nx, x.data(), DF_NO_HINT,
                          ny, v.data(), DF_NO_HINT);
    if (status != DF_STATUS_OK)
    {
        isclPrintError("%s", "Error creating spline task");
        ierr = ISCL_MKL_FAILURE;
        goto ERROR;
    }
    // Set spline parameters in the data fitting task
    status = dfdEditPPSpline1D(task, splineOrder, splineType, splineBC, bc,
                               splineIC, ic, scoeff, DF_NO_HINT);
    if (status != DF_STATUS_OK)
    {
        isclPrintError("%s", "Error setting spline parameters");
        ierr = ISCL_MKL_FAILURE;
        goto ERROR;
    }
    // Check if outputs are sorted
    bool lsorted = true; // If I am to assume its sorted then default to true
    if (!assumeSorted){lsorted = VectorMath::isSorted(xq);}
    // Compute the cell indices
    MKL_INT *cell = static_cast<MKL_INT *>
                    (mkl_malloc((size_t) nq*sizeof(MKL_INT), 64));
    if (lsorted)
    {
        status = dfdSearchCells1D(task, DF_METHOD_STD, nsite, xqwork,
                                 DF_NO_HINT, DF_NO_APRIORI_INFO, cell);
    }
    else
    {
        status = dfdSearchCells1D(task, DF_METHOD_STD, nsite, xqwork,
                                 DF_SORTED_DATA, DF_NO_APRIORI_INFO, cell);
    }
    if (status != DF_STATUS_OK)
    {
        delete[] cell; 
    }
    // Do interpolation
    vq.resize(nq);
    if (method != Interpolate::Method::NEAREST)
    {
        const MKL_INT dorder[1] = {0};  // order of derivatives to be computed
        const int nsite = static_cast<MKL_INT> (nq);
        if (lsorted)
        {
            status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
                                      nsite, xqwork,
                                      DF_SORTED_DATA, norder, dorder,
                                      DF_NO_APRIORI_INFO, vq.data(),
                                      DF_MATRIX_STORAGE_ROWS, cell);
        }
        else
        {
            status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
                                      nsite, xqwork,
                                      DF_NO_HINT, norder, dorder,
                                      DF_NO_APRIORI_INFO, vq.data(),
                                      DF_MATRIX_STORAGE_ROWS, cell);
        }
        if (status != DF_STATUS_OK)
        {
            vq.resize(0);
            delete[] cell;
            dfDeleteTask(&task);
            throw std::invalid_argument("This shouldn't happen");
        }
    }
    else
    {
        for (i=0; i<nq; i++)
        {
            indx = static_cast<int> (cell[i]);
            vq[i] = v[indx];
            if (std::abs(xqwork[i] - x[indx]) > std::abs(xqwork[i] - x[indx-1]))
            {
                vq[i] = v[indx-1];
            }
        }
    }
    delete[] cell;
    dfDeleteTask(&task);
*/
/*
    return;
}
*/
