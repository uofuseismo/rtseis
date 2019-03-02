#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <mkl.h>
#include <ipps.h>
#include "rtseis/utilities/math/interpolate.hpp"


void RTSeis::Utilities::Math::Interpolate::interpft(
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
    IppsDFTSpec_R_64f *pDFTForwardSpec
        = (IppsDFTSpec_R_64f *) ippsMalloc_8u(specSizeF);
    Ipp8u *pDFTInitBuf = ippsMalloc_8u(std::max(sizeInitF, sizeInitI));
    status = ippsDFTInit_R_64f(npts, 
                               IPP_FFT_DIV_INV_BY_N,
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
    IppsDFTSpec_R_64f *pDFTInverseSpec
        = (IppsDFTSpec_R_64f *) ippsMalloc_8u(specSizeI);
    status = ippsDFTInit_R_64f(npnew,
                               IPP_FFT_DIV_INV_BY_N,
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
    // Remove scaling from FFT
    double xscal = static_cast<double> (npnew)/static_cast<double> (npts);
    ippsMulC_64f_I(xscal, yint.data(), npnew);
    // Clean up
    ippsFree(pBuf);
    ippsFree(pDst);
    ippsFree(pDFTForwardSpec);
    ippsFree(pDFTInverseSpec);
    return;
}

/*

void Interpolate::interp1d(const std::vector<double> &x,
                           const std::vector<double> &v,
                           const std::vector<double> &xq,
                           std::vector<double> &vq,
                           const Interpolate::Method method, 
                           const bool lxSorted,
                           const bool lxqSorted)
{
    // Check that there's actually something to do
    int nq = xq.size();
    vq.resize(0);
    if (nq <= 0){return 0;} // Nothing to do
    // Make sure v = f(x) makes sense in terms of vector size
    int n = x.size();
    if (x.size() != v.size())
    {
        throw std::invalid_argment("x size != v size");
    }
    // Figure out the interpolation options for MKL
    if (method == Interpolate::Method::NEAREST)
    {
        splineOrder = DF_PP_STD;
        splineType  = DF_CR_STEPWISE_CONST_INTERPOLANT;
        splineBC    = DF_NO_BC;
        splineIC    = DF_NO_IC;
    }
    else if (method == Interpolate::Method::LINEAR)
    {
        splineOrder = DF_PP_LINEAR;
        splineType  = DF_PP_LINEAR;
        splineBC    = DF_NO_BC;
        splineIC    = DF_NO_IC;
        if (n < 3)
        {
            throw std::invalid_argument("LinInt requires at least 2 points");
        }
    }
    else if (method == Interpolate::Method::CSPLINE_NATURAL)
    {
        splineOrder = DF_PP_CUBIC;
        splineType  = DF_PP_NATURAL;
        splineBC    = DF_BC_FREE_END;
        splineIC    = DF_NO_IC;
        if (n < 4)
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
        if (n < 4)
        {
            throw std::invalid_argument("Splines require at least 3 points");
        } 
        if (std::abs(v[0] - v[n-1]) > 1.e-14)
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
        if (n < 4)
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
        if (n < 4)
        {   
            throw std::invalid_argument("Splines require at least 3 points");
        }
        if (std::abs(v[0] - v[n-1]) > 1.e-14)
        {
            throw std::invalid_argument("v[0] != v[n-1]");
        }
    }
    else
    {
        throw std::invalid_argument("Unsupported interpolation type"); 
    }
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
    return;
}
*/
