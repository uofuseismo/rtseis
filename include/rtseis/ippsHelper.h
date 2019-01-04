#ifndef RTSEIS_IPPSHELPER_H__
#define RTSEIS_IPPSHELPER_H__
#include <stdbool.h>
#include "rtseis/enums.h"

/*!
 * @defgroup rtseis_ipps IPP Interface
 * @brief Interfaces to the Intel Performance Primitives library.  IPP is
 *        the workhorse that drives ISPL.  Most users will likely never
 *        call the functions in this module directly.
 * @copyright ISTI distributed under the Apache 2 license.
 */

struct ippsDFTR2C_struct
{
    /*!
     * @struct ippsDFTR2C_struct
     * @brief DFT structure with variables to be used by IPP.
     * @ingroup rtseis_ipps_dft
     */
    void *ftHandle;  /*!< Holds the DFT or FFT handle. */
    void *xwork;     /*!< Workspace to hold the signal to transform.  This is an
                          Ipp64f or Ipp32f array of dimension [2*length]. */
    void *ywork;     /*!< Workspace to hold the transform.   This is an Ipp64f
                          or Ipp32f array of dimension [2*length]. */ 
    void *fftBuffer; /*!< An additional buffer for the FFT. */
    void *pBuf;      /*!< Ipp8u buffer for DFT.  This is an array of dimension
                          [bufferSize]. */
    int length;      /*!< Length of the signal to transform. */
    int lenft;       /*!< Length of the Fourier transform: length/2 + 1. */
    int specSize;    /*!< Size of the FFT handle. */ 
    int order;       /*!< For the base 2 FFT the length is \f$ 2^{order} \f$. */
    int bufferSize;  /*!< Size of pBuf. */
    enum rtseisPrecision_enum precision; /*!< Precision. */
    bool luseFFT;    /*!< If true then use the FFT.  Otherwise, use the DFT. */
    bool linit;      /*!< If true then the structure is initialized. */
};

struct ippsIIRFilter_struct
{
    /*!
     * @struct ippsIIRFilter_struct
     * @brief IIR filter structure with variables to be used by IPP.
     * @ingroup rtseis_ipps_iir
     */
    void *bufsrc;   /*!< This is an array of dimension [2 x 1024]. */
    void *bufdst;   /*!< This is an array of dimension [2 x 1024]. */
    void *bufdly;   /*!< This is an array of dimension [nbDly]. */
    void *pTaps;    /*!< This is an Ipp64f or Ipp32f array of length
                         [2*(order+1)] holding the feedforward then 
                         feedbackward coefficients. */
    void *pBuf;     /*!< This is an Ipp8u array of length bufferSize. */
    void *state;    /*!< IIR state structure. */
    void *zi;       /*!< Initial conditions. */
    double *bRef;   /*!< Reference numerator coefficients.  This has
                         dimension [nbRef]. */
    double *aRef;   /*!< Reference denominator coefficients.  This has
                         dimension [naRef]. */
    int bufferSize; /*!< Length of pBuf. */
    int nbDly;      /*!< Length of bufldy.  This is equal to
                         MAX(128, order+1). */
    int nbRef;      /*!< Reference number of input numerator coefficients. */
    int naRef;      /*!< Input reference number of denominator coefficients. */
    int order;      /*!< Filter order.  This is an array of dimension
                         [order]. */
    int accumLen;   /*!< Number samples mod 1024 accumulated. */
    enum rtseisPrecision_enum precision; /*!< Precision. */
    bool lnextBlock; /*!< If true then the initial conditions need to be
                          set in IPP >= v.9. */
    bool lrt;   /*!< True if the filter is for real-time applications. */
    bool linit; /*!< True if the structure has been initialized. */
    bool useLegacy; /*!< IPP made a pretty weird switch in 2019 with regards
                         to delay lines.  This is a flag indicating one
                         should use the legacy source or the new source. */
};

struct ippsZPIIRFilter_struct
{
    /*!
     * @struct ippsZPIIRFilter_struct 
     * @brief Zero-phase IIR filter structure with variables to be used by 
IPP.
     * ingroup rtseis_ipps_zpiir
     */
    void *pTaps;    /*!< This is an Ipp64f or Ipp32f array of length
                         [2*(order+1)] holding the feedforward then 
                         feedbackward coefficients. */
    void *pBuf;     /*!< This is an Ipp8u array of length bufferSize. */
    void *state;    /*!< IIR state structure. */
    void *zi;       /*!< Initial conditions. */
    double *bRef;   /*!< Reference numerator coefficients.  This is an
                         array of dimension [nbRef]. */
    double *aRef;   /*!< Reference denominator ceofficients.  This is an
                         array of dimension [naRef]. */
    int bufferSize; /*!< Length of pBuf. */
    int nbDly;      /*!< Length of bufldy.  This is equal to
                         MAX(128, order+1). */
    int order;      /*!< Filter order.  This is an array of dimension
                         [order]. */
    int nbRef;      /*!< Number of reference numerator coefficients. */
    int naRef;      /*!< Number of reference denominator coefficients. */
    enum rtseisPrecision_enum precision; /*!< Precision. */
    bool lhaveZI;   /*!< If true then the filter has initial conditions. 
                         Otherwise, zero initial conditions are assumed. */
    bool linit; /*!< True if the structure has been initialized. */
    char pad[2];
};

struct ippsFIRFilter_struct
{
    /*!
     * @struct ippsFIRFilter_struct
     * @brief FIR filter structure with variables to be used by IPP.
     * @ingroup rtseis_ipps_fir
     */
    void *dlysrc;   /*!< Input delay line (i.e., the initial conditions). 
                         This is an array of dimension [nbDly]. */
    void *dlydst;   /*!< Output delay line (i.e., the final conditions). 
                         This is an array of dimension [nbDly]. */ 
    void *pBuf;     /*!< This is an Ipp8u array of length bufferSize. */
    void *spec;     /*!< FIR state structure. */
    void *pTaps;    /*!< This is an Ipp64f or Ipp32f array of length
                         [order + 1] holding the numerator coefficients. */
    double *zi;     /*!< Initial conditions.  This is an array of dimension
                         [order]. */
    double *tapsRef;/*!< Input taps.  This is an array of dimension of
                         order [order + 1] and allows for save copying
                         of filter taps without loss of precision. */
    int bufferSize; /*!< Length of pBuf. */
    int specSize;   /*!< Size of state. */
    int order;      /*!< Filter order. */
    int nbDly;      /*!< Length of delay lines.  This is equal to 
                         MAX(128, order+1). */
    enum rtseisPrecision_enum precision; /*!< Precision. */
    bool lrt;       /*!< True if the filter is for real-time applications. */ 
    bool linit;     /*!< True if the structure has been initialized. */
    char pad[2];
};

struct ippsFIRMRFilter_struct
{
    /*!
     * @struct ippsFIRMRFilter_struct
     * @brief FIR multi-rate filter structure with variables to be used by IPP.
     * @ingroup rtseis_ipps_firmr
     */
    void *work;     /*!< Ipp64 or Ipp32 workspace for holding excess samples
                         in input signal. This is an array of dimension
                         [downFactor]. */
    void *pSrc;     /*!< Ipp64 or Ipp32 workspace for holding input
                         signal.  This is an array of dimension 
                         [nwork]. */
    void *pDst;     /*!< ipps64 or ipps32 workspace for holding output
                         signal.  This is an array of dimension
                         [nwork]. */ 
    void *dlysrc;   /*!< Input delay line (i.e., the initial conditions). 
                         This is an array of dimension [mbDly]. */
    void *dlydst;   /*!< Output delay line (i.e., the final conditions). 
                         This is an array of dimension [mbDly]. */ 
    void *pBuf;     /*!< This is an Ipp8u array of length bufferSize. */
    void *spec;     /*!< FIR state structure. */
    void *pTaps;    /*!< This is an Ipp64f or Ipp32f array of length
                         [order + 1] holding the numerator coefficients. */
    double *zi;      /*!< Initial conditions.  This is an array of dimension
                         [mbDly]. */
    /*!< A copy of the filter taps.  This has dimension [tapsLen]. */
    double *tapsRef;
    int bufferSize; /*!< Length of pBuf. */
    int specSize;   /*!< Size of state. */
    int order;      /*!< Filter order. */
    int nbDly;      /*!< Length of delay lines.  This is equal to 
                         (tapsLen + upFactor - 1)/upFactor. */
    int mbDly;      /*!< Max workspace for delay lines.  This is equal to
                         MAX(128, nbDly). */
    int upFactor;   /*!< Upsampling factor. */
    int downFactor; /*!< Downsampling factor. */
    int upPhase;    /*!< Phase of up-sampler. */
    int downPhase;  /*!< Phase of down-sampler. */
    int tapsLen;    /*!< Number of taps.  This is order + 1. */
    int nwork;      /*!< Size of pSrc and pDst. */
    int nexcess;    /*!< Number of excess samples retained from previous
                         application. */
    enum rtseisPrecision_enum precision; /*!< Precision. */
    bool lrt;       /*!< True if the filter is for real-time applications. */ 
    bool linit;     /*!< True if the structure has been initialized. */
    char pad[2];
};

struct ippsSOSFilter_struct
{
    /*!
     * @struct ippsSOSFilter_struct
     * @brief SOS filter structure with variables to be used by IPP.
     * @ingroup rtseis_ipps_sos
     */
    void *dlysrc;   /*!< Input delay line (i.e., the initial condtions).
                         This is an array of dimension [nbDly]. */
    void *dlydst;   /*!< Output delay line (i.e., the final conditions).
                         This is an array of dimension [nbDly]. */
    void *pBuf;     /*!< This is an Ipp8u workspace array of dimension
                         [bufferSize]. */
    void *state;    /*!< SOS filter state. */
    void *pTaps;    /*!< This is an Ipp64f or Ipp32f array containing the
                         filter taps.  It has length [6*ns]. */
    void *zi;       /*!< Initial conditions.  This is an array of dimension
                         [nbDly]*/ 
    double *bsRef;  /*!< Reference numerator coefficients.  This is an array
                         of dimension [3 x ns]. */
    double *asRef;  /*!< Reference denominator coefficients.  This is an array
                         of dimension [3 x ns]. */
    int ns;         /*!< Number of sections. */
    int nbDly;      /*!< This is MAX(2*ns, 128). */ 
    int bufferSize; /*!< Size of pBuf. */
    enum rtseisPrecision_enum precision; /*!< Precision. */
    bool lrt;       /*!< True if the filter is for real-time applications. */
    bool linit;     /*!< True if the structure has been initialized. */
};

struct ippsDownsample_struct
{
    /*!
     * @struct ippsDownsample_struct
     * @brief Downsampling structure with variables to be used by IPP.
     * @ingroup rtseis_ipps_downsample
     */
    int downFactor; /*!< Downsampling factor. */
    int phase;      /*!< Phase of downsampler. */
    int phase0;     /*!< Initial phase of downsampler. */
    enum rtseisPrecision_enum precision; /*!< Precision. */
    bool lrt;       /*!< True if the filter is for real-time applications. */
    bool linit;     /*!< True if the structure has been initialized. */
    char pad[6];
};

struct ippsSparseFIRFilter_struct
{
    /*!
     * @struct ippsSparseFIRFilter_struct
     * @brief Sparse FIR filter structure with variables to be used by IPP.
 *   * @ingroup rtseis_ipps_sparsefir
     */
    void *dlysrc;    /*!< Input delay line (i.e., the initial conditions).
                          This is an array of dimension [nbDly]. */
    void *dlydst;    /*!< Output delay line (i.e., the final conditions).
                          This is an array of dimension [nbDly]. */
    void *pState;    /*!< Sparse FIR state structure. */
    double *zi;      /*!< Initial conditions in the sparse delay line.  This is
                          an array of dimension [order]. */
    void *pNZTaps;   /*!< Non-zero filter tap coefficients.  This is an array
                          of dimension [nnzCoeffs]. */
    void *pNZTapPos; /*!< Indices corresponding to the non-zero taps.  This is
                         an array of dimension [nnzCoeffs]. */ 
    void *pBuf;      /*!< Workspace for sparse FIR filtering.  This is an 
                          ipps8u array of dimension [bufferSize]. */
    double *tapsRef; /*!< Copy of the input taps.  This is an array of dimension
                          [nnzCoeffs]. */
    int bufferSize;  /*!< Length of pBuf. */
    int order;       /*!< Order of the sparse filter. */
    int nbDly;       /*!< This is MAX(256, order+1). */
    int nnzCoeffs;   /*!< Number of non-zero filter taps. */
    enum rtseisPrecision_enum precision; /*!< Precision. */
    bool lrt;        /*!< True if the filter is for real-time applications. */
    bool linit;      /*!< True if the structure has been initialized. */
};

struct ippsMedianFilter_struct
{
    /*!
     * @struct ippsMedianFilter_struct
     * @brief Median filter structure with variables to be used by IPP.
     * @ingroup rtseis_ipps_median
     */
    void *dlysrc;   /*!< The filter initial conditions.  This is a float or
                         double array of dimension [MAX(1, maskSize-1)]. */
    void *dlydst;   /*!< The filter final conditions.  This is a float or double
                         array of dimension [MAX(1, maskSize-1)]. */
    void *pBuf;     /*!< Workspace for median.  This is an array of dimension
                         [bufferSize]. */
    void *zi;       /*!< Initial conditions.  This is an array of dimension
                         [MAX(1, maksSize-1)] but only the first maskSize-1
                         elements are accessed. */
    int bufferSize; /*!< The length of pBuf. */
    int maskSize;   /*!< Length of window over which the median is computed. */
    enum rtseisPrecision_enum precision; /*!< Precision of filter. */
    bool lrt;       /*!< Flag indicating whether the filter is for real-time
                         or post-processing. */
    bool linit;     /*!< If true the structure is initialized. */
};

#ifdef __cplusplus
extern "C"
{
#endif

/*----------------------------------------------------------------------------*/
/*                                  Downsampling                              */
/*----------------------------------------------------------------------------*/
/* Initialization function */
int rtseis_ippsDownsample_initialize(const int downFactor,
                                     const bool lisRealTime,
                                     const enum rtseisPrecision_enum precision,
                                     struct ippsDownsample_struct *ippsDS);
/* Apply the downsampler */
int rtseis_ippsDownsample_apply64f(const int n,
                                   const double x[],
                                   const int ny, int *len,
                                   double y[],
                                   struct ippsDownsample_struct *ippsDS);
int rtseis_ippsDownsample_apply32f(const int n,
                                   const float x[],
                                   const int ny, int *len,
                                   float y[],
                                   struct ippsDownsample_struct *ippsDS);
/* Space estimate. */
int rtseis_ippsDownsample_estimateSpace(
    const int n,
    const struct ippsDownsample_struct ippsDS);
/* Toggle the filter to be for real-time processing or for post-processing */
int rtseis_ippsDownsample_toggleRealTime(const bool lisRealTime,
                                         struct ippsDownsample_struct *ippsDS);
/* Finalize the downsampler */
int rtseis_ippsDownsample_finalize(struct ippsDownsample_struct *ippsDS);
/* Resets initial conditions */
int rtseis_ippsDownsample_resetInitialConditions(
    struct ippsDownsample_struct *ippsDS);
int rtseis_ippsDownsample_setInitialConditions(
    const int phase, struct ippsDownsample_struct *ippsDS);
/*----------------------------------------------------------------------------*/
/*                                 IIR Filtering                              */
/*----------------------------------------------------------------------------*/
/* Initialization function */
int rtseis_ippsIIRFilter_initialize(const int nb, const double b[],
                                   const int na, const double a[],
                                   const bool lisRealTime,
                                   const enum rtseisPrecision_enum precision,
                                   struct ippsIIRFilter_struct *ippIIR);
/* Apply IIR filter with IPP */
int rtseis_ippsIIRFilter_apply64f(const int n, const double x[],
                                  double y[],
                                  struct ippsIIRFilter_struct *ippIIR);
int rtseis_ippsIIRFilter_apply32f(const int n, const float x[],
                                 float y[],
                                 struct ippsIIRFilter_struct *ippsIIR);
/* Copies an IIR filter. */
int rtseis_ippsIIRFIlter_copy(const struct ippsIIRFilter_struct ippsIIRIn,
                              struct ippsIIRFilter_struct *ippsIIROut);
/* Toggle the filter to be for real-time processing or for post-processing */
int rtseis_ippsIIRFilter_toggleRealTime(const bool lisRealTime,
                                        struct ippsIIRFilter_struct *ippsIIR);
/* Finalize IPP IIR filter */
int rtseis_ippsIIRFilter_finalize(struct ippsIIRFilter_struct *ippsIIR);
/* Resets initial/final conditions */
int rtseis_ippsIIRFilter_setInitialConditions(
    const int nz, const double zi[], struct ippsIIRFilter_struct *ippsIIR);
int rtseis_ippsIIRFilter_resetInitialConditions(
    struct ippsIIRFilter_struct *ippsIIR);
/*----------------------------------------------------------------------------*/
/*                            Zero-Phase IIR Filtering                        */
/*----------------------------------------------------------------------------*/
/* Initialization function */
int rtseis_ippsZPIIRFilter_initialize(const int nb, const double b[],
                                      const int na, const double a[],
                                      const enum rtseisPrecision_enum precision,
                                      struct ippsZPIIRFilter_struct *ippsZPIIR);
/* Finalize IPP zero-phase IIR filter */
int rtseis_ippsZPIIRFilter_finalize(struct ippsZPIIRFilter_struct *ippsZPIIR);
/* Apply zero-phase IIR filter */
int rtseis_ippsZPIIRFilter_apply64f(const int n, const double x[],
                                    double y[],
                                    struct ippsZPIIRFilter_struct *ippsZPIIR);
int rtseis_ippsZPIIRFilter_apply32f(const int n, const float x[],
                                    float y[],
                                    struct ippsZPIIRFilter_struct *ippsZPIIR);
/* Copies a zero-phase IIR filter structure. */
int rtseis_ippsZPIIRFilter_copy(const struct ippsZPIIRFilter_struct ippsZPIIRIn,
                                struct ippsZPIIRFilter_struct *ippsZPIIROut);
/* Finalize the IPP zero-phase IIR filter */
int rtseis_ippsZPIIRFilter_finalize(struct ippsZPIIRFilter_struct *ippsZPIIR);
/*----------------------------------------------------------------------------*/
/*                               FIR Filtering                                */
/*----------------------------------------------------------------------------*/
/* Initialization function */
int rtseis_ippsFIRFilter_initialize(const int nb, const double b[],
                                    const bool lisRealTime,
                                    const enum rtseisPrecision_enum precision,
                                    struct ippsFIRFilter_struct *ippFIR);
/* Copies an FIR structure */
int rtseis_ippsFIRFilter_copy(const struct ippsFIRFilter_struct ippsFIRIn,
                              struct ippsFIRFilter_struct *ippsFIROut);
/* Resets the initial conditions */
int rtseis_ippsFIRFilter_resetInitialConditions(
    struct ippsFIRFilter_struct *ippsFIR);
/* Sets the initial conditions */
int rtseis_ippsFIRFilter_setInitialConditions(
    const int nz, const double zi[],
    struct ippsFIRFilter_struct *ippsFIR);
/* Apply the filter */
int rtseis_ippsFIRFilter_apply64f(const int n,
                                  const double x[],
                                  double y[],
                                  struct ippsFIRFilter_struct *ippsFIR);
int rtseis_ippsFIRFilter_apply32f(const int n,
                                  const float x[],
                                  float y[],
                                  struct ippsFIRFilter_struct *ippsFIR);
/* Determines if the filter is for real-time or post-processing. */
bool rtseis_ippsFIRFilter_isRealTime(const struct ippsFIRFilter_struct ippsFIR);
/* Toggle the filter to be for real-time processing or for post-processing */
int rtseis_ippsFIRFilter_toggleRealTime(const bool lisRealTime,
                                        struct ippsFIRFilter_struct *ippsFIR);
/* Finalize IPP FIR filter */
int rtseis_ippsFIRFilter_finalize(struct ippsFIRFilter_struct *ippsFIR);
/*----------------------------------------------------------------------------*/
/*                          Multirate FIR Filtering                           */
/*----------------------------------------------------------------------------*/
/* Initialization function */
int rtseis_ippsFIRFilterMR_initialize(
    const int upFactor, const int downFactor,
    const int nb, const double b[],
    const bool lisRealTime,
    const enum rtseisPrecision_enum precision,
    struct ippsFIRMRFilter_struct *ippsFIRMR);
/* Copies the multirate FIR filter. */
int rtseis_ippsFIRFilterMR_copy(
    const struct ippsFIRMRFilter_struct ippsFIRMRIn,
    struct ippsFIRMRFilter_struct *ippsFIRMROut);
/* Resets the initial conditions */
int rtseis_ippsFIRFilterMR_resetInitialConditions(
    struct ippsFIRMRFilter_struct *ippsFIRMR);
/* Sets the initial conditions */
int rtseis_ippsFIRFilterMR_setInitialConditions(
    const int nz, const double zi[],
    struct ippsFIRMRFilter_struct *ippsFIRMR);
/* Apply the filter */
int rtseis_ippsFIRFilterMR_apply64f(const int n,
                                    const double x[],
                                    const int ny, int *len,
                                    double y[],
                                    struct ippsFIRMRFilter_struct *ippsFIRMR);
int rtseis_ippsFIRFilterMR_apply32f(const int n,
                                    const float x[],
                                    const int ny, int *len,
                                    float y[],
                                    struct ippsFIRMRFilter_struct *ippsFIRMR);
/* Toggle the filter to be for real-time processing or for post-processing */
int rtseis_ippsFIRFilterMR_toggleRealTime(
    const bool lisRealTime,
    struct ippsFIRMRFilter_struct *ippsFIRMR);
/* Finalize the IPP FIRMR filter */
int rtseis_ippsFIRFilterMR_finalize(struct ippsFIRMRFilter_struct *ippsFIRMR);
/*----------------------------------------------------------------------------*/
/*                           Sparse FIR Filtering                             */
/*----------------------------------------------------------------------------*/
/* Initializes the sparse FIR filtering structure */
int rtseis_ippsSparseFIR_initialize(
    const int nnzCoeffs,
    const int nzIndices[],
    const double nzCoeffs[],
    const bool lisRealTime,
    const enum rtseisPrecision_enum precision,
    struct ippsSparseFIRFilter_struct *ippsSparseFIR);
/* Toggle the filter to be for real-time processing or for post-processing */
int rtseis_ippsSparseFIR_toggleRealTime(
    const bool lisRealTime,
    struct ippsSparseFIRFilter_struct *ippsSparseFIR);
/* Applies the sparse FIR filtering */
int rtseis_ippsSparseFIR_apply64f(
    const int n, const double x[], double y[],
    struct ippsSparseFIRFilter_struct *ippsSparseFIR);
int rtseis_ippsSparseFIR_apply32f(
    const int n, const float x[], float y[],
    struct ippsSparseFIRFilter_struct *ippsSparseFIR);
/* Sets the initial conditions */
int rtseis_ippsSparseFIR_setInitialConditions(
    const int nz, const double zi[],
    struct ippsSparseFIRFilter_struct *ippsSparseFIR);
/* Resets the filter initial conditions */
int rtseis_ippsSparseFIR_resetInitialConditions(
    struct ippsSparseFIRFilter_struct *ippsSparseFIR);
/* Frees the sparse FIR filtering structure */
int rtseis_ippsSparseFIR_finalize(
    struct ippsSparseFIRFilter_struct *ippsSparseFIR);
/*----------------------------------------------------------------------------*/
/*                          SOS (BiQuad) Filtering                            */
/*----------------------------------------------------------------------------*/
/* Initialization function. */
int rtseis_ippsSOSFilter_initialize(const int ns, 
                                    const double bs[],
                                    const double as[],
                                    const bool lisRealTime,
                                    const enum rtseisPrecision_enum precision,
                                    struct ippsSOSFilter_struct *ippsSOS);
/* Resets the initial conditions. */
int rtseis_ippsSOSFilter_resetInitialConditions(
    struct ippsSOSFilter_struct *ippsSOS);
/* Sets the initial conditions. */
int rtseis_ippsSOSFilter_setInitialConditions(
    const int nz, const double zi[],
    struct ippsSOSFilter_struct *ippsSOS);
/* Toggle the filter to be for real-time processing or for post-processing. */
int rtseis_ippsSOSFilter_toggleRealTime(const bool lisRealTime,
                                        struct ippsSOSFilter_struct *ippsSOS);
/* Finalize the filter. */
int rtseis_ippsSOSFilter_finalize(struct ippsSOSFilter_struct *ippsSOS);
/* Apply the filter. */
int rtseis_ippsSOSFilter_apply64f(const int n,
                                  const double x[],
                                  double y[],
                                  struct ippsSOSFilter_struct *ippsSOS);
int rtseis_ippsSOSFilter_apply32f(const int n,
                                const float x[],
                                float y[],
                                struct ippsSOSFilter_struct *ippsSOS);

/*----------------------------------------------------------------------------*/
/*                              Median Filtering                              */
/*----------------------------------------------------------------------------*/
/* Initialize the filter */
int rtseis_ippsMedianFilter_initialize(
    const int n,
    const bool lisRealTime,
    const enum rtseisPrecision_enum precision,
    struct ippsMedianFilter_struct *ippsMedian);
/* Toggle the filter for real-time or post-processing. */
int rtseis_ippsMedianFilter_toggleRealTime(
    const bool lisRealTime,
    struct ippsMedianFilter_struct *ippsMedian);
/* Resets the initial conditions on the median filter. */
int rtseis_ippsMedianFilter_resetInitialConditions(
    struct ippsMedianFilter_struct *ippsMedian);
/* Sets the initial conditions on the median filter. */
int rtseis_ippsMedianFilter_setInitialConditions(
    const int nz, const double zi[],
    struct ippsMedianFilter_struct *ippsMedian);
/* Apply the median filter */
int rtseis_ippsMedianFilter_apply64f(
    const int n,
    const double x[],
    double y[],
    struct ippsMedianFilter_struct *ippsMedian);
int rtseis_ippsMedianFilter_apply32f(
    const int n,
    const float x[],
    float y[],
    struct ippsMedianFilter_struct *ippsMedian);
/* Free the median filter */
int rtseis_ippsMedianFilter_finalize(
    struct ippsMedianFilter_struct *ippsMedian);
/*----------------------------------------------------------------------------*/
/*                           Discrete Fourier Transform                       */
/*----------------------------------------------------------------------------*/
/* Initialize the DFT. */
int rtseis_ippsDFTR2C_initialize(const bool luseFFT,
                                 const int length,
                                 const enum rtseisPrecision_enum precision,
                                 struct ippsDFTR2C_struct *ippsDFTR2C);
/* Frees the DFT. */
int rtseis_ippsDFTR2C_finalize(struct ippsDFTR2C_struct *ippsDFTR2C);
/* Transform from real to complex */
int rtseis_ippsDFTR2C_applyForwardDFT64f(const int n,
                                         const double x[],
                                         const int maxy,
                                         int *ny,
                                         void *y, // double complex y[]
                                         struct ippsDFTR2C_struct *ippsDFTR2C);
int rtseis_ippsDFTR2C_applyForwardDFT32f(const int n,
                                         const float x[],
                                         const int maxy,
                                         int *ny,
                                         void *y, //float complex y[],
                                         struct ippsDFTR2C_struct *ippsDFTR2C);
/* Inverse transform from complex to real. */
int rtseis_ippsDFTR2C_applyInverseDFT64f(const int lenft,
                                         const void *x, //double complex x[],
                                         const int ny, 
                                         double y[],
                                         struct ippsDFTR2C_struct *ippsDFTR2C);
int rtseis_ippsDFTR2C_applyInverseDFT32f(const int lenft,
                                         const void *x, //float complex x[],
                                         const int ny, 
                                         float y[],
                                         struct ippsDFTR2C_struct *ippsDFTR2C);


#ifdef __cplusplus
}
#endif
#endif
