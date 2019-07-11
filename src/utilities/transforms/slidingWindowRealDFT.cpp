#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <complex> // Put this before fftw
#include <fftw/fftw3.h>
#include <ipps.h>
#include <rtseis/private/throw.hpp>
#include "rtseis/utilities/transforms/slidingWindowRealDFT.hpp"
#include "rtseis/utilities/transforms/slidingWindowRealDFTParameters.hpp"
#include "rtseis/utilities/filterImplementations/detrend.hpp"


using namespace RTSeis::Utilities::Transforms;

namespace
{

inline int padLength64f(const int n, const int alignment=64)
{
    auto size = static_cast<int> (sizeof(double));
    int padLength = 0;
    auto xmod = (n*size)%alignment;
    if (xmod != 0){padLength = (alignment - xmod)/size;}
    auto nptsPadded = n + padLength;
    return nptsPadded;
}

inline int padLength32f(const int n, const int alignment=64)
{
    auto size = static_cast<int> (sizeof(float));
    int padLength = 0;
    auto xmod = (n*size)%alignment;
    if (xmod != 0){padLength = (alignment - xmod)/size;}
    auto nptsPadded = n + padLength;
    return nptsPadded;
}

}

class SlidingWindowRealDFT::SlidingWindowRealDFTImpl
{
public:
    /// Constructor
    SlidingWindowRealDFTImpl() = default;
    /// Destructor
    ~SlidingWindowRealDFTImpl()
    {
        clear();
    }
    /// Releases memory
    void clear()
    {
        mParameters.clear();
        if (mHaveDoublePlan){fftw_destroy_plan(mDoublePlan);}
        if (mHaveFloatPlan){fftwf_destroy_plan(mFloatPlan);}
        if (mOutData64f != nullptr){fftw_free(mOutData64f);}
        if (mOutData32f != nullptr){fftwf_free(mOutData32f);}
        if (mWindow64f != nullptr){ippsFree(mWindow64f);}
        if (mWindow32f != nullptr){ippsFree(mWindow32f);}
        if (mInData64f != nullptr){fftw_free(mInData64f);}
        if (mInData32f != nullptr){fftwf_free(mInData32f);}
        mOutData64f = nullptr;
        mOutData32f = nullptr;
        mWindow64f = nullptr;
        mWindow32f = nullptr;
        mInData64f = nullptr;
        mInData32f = nullptr;
        mInDataLength = 0;
        mOutDataLength = 0;
        mSamples = 0;
        mSamplesPerSegment = 0;
        mSamplesInOverlap = 0;
        mDFTLength = 0;
        mNumberOfFrequencies = 0;
        mNumberOfColumns = 0;
        mDataOffset = 0;
        mFTOffset = 0;
        mPrecision = RTSeis::Precision::DOUBLE;
        mDetrendType = SlidingWindowDetrendType::REMOVE_NONE;
        mHaveDoublePlan = false;
        mHaveFloatPlan = false;
        mApplyWindow = false;
        mHaveTransform = false;
        mInitialized = false;
    }

//private:
    /// The parameters that went into initialization
    class SlidingWindowRealDFTParameters mParameters;
    /// FFTw plan
    fftw_plan  mDoublePlan;
    /// Holds the data to Fourier transform.  This is an array of dimension
    /// [mInDataOffset x mNumberOfColumns]
    double *mInData64f = nullptr;
    /// Holds the Fourier transformed data.  This is an array of dimension
    /// [mOutDataOffset x mNumberOfColumns]
    fftw_complex *mOutData64f = nullptr;
    /// Holds the window function.  This is an array of dimension
    /// [mSamplesPerSegment].  This is used with mApplyWindow.
    double *mWindow64f = nullptr;
    /// Holds the FFTw floating arithmetic plan
    fftwf_plan mFloatPlan;
    /// Holds the data to Fourier transform.  This is an array of dimension
    /// [mInDataOffset x mNumberOfColumns]
    float *mInData32f = nullptr;
    /// Holds the Fourier transformed data.  This is an array of dimension
    /// [mOutDataOffset x mNumberOfColumns]
    fftwf_complex *mOutData32f = nullptr;
    /// Holds the window function.  This is an array of dimension
    /// [mSamples].  This is used with mApplyWindow.
    float *mWindow32f = nullptr;
    /// The number of samples in the input signal
    int mSamples = 0;
    /// The number of samples in an overlap
    int mSamplesInOverlap = 0;
    /// The length of array mInData = mInDataOffset x mNumberOfColumns
    int mInDataLength = 0;
    /// The length of array mOutData = mOutDataOffset x mNumberOfColumns
    int mOutDataLength = 0;
    /// The length of the DFT.  This is for padding purposes.
    int mDFTLength = 0;
    /// The number of frequencies. = mDFTLength/2 = 1
    int mNumberOfFrequencies = 0;
    /// The number of time domain samples where the sliding DFT is tabulated.
    int mNumberOfColumns = 0;
    /// This is a padded variant of mDFTLength.  The padding ensures 64 bit
    /// cache alignment for each row of the input data to transform.
    int mDataOffset = 0;
    /// This is a padded variant of mNumberOfFrequencies.  The padding ensures
    /// 64 bit cache alignment ofr each row of the output data to transform.
    int mFTOffset = 0;
    /// The number of samples in a segment
    int mSamplesPerSegment = 0; 
    /// The precision of the module
    RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
    /// The detrend strategy
    SlidingWindowDetrendType mDetrendType
       = SlidingWindowDetrendType::REMOVE_NONE;
    /// Flag indicating that I have a plan (for double precision)
    bool mHaveDoublePlan = false;
    /// Flag indicating that I have a plan (for single precision)
    bool mHaveFloatPlan = false;
    /// Flag indicating whether or not I will apply the window function.
    bool mApplyWindow = false;
    /// Flag indicating the transform was applied
    bool mHaveTransform = false;
    /// Flag indicating the class is inititalized.
    bool mInitialized = false;
};

/// Constructor
SlidingWindowRealDFT::SlidingWindowRealDFT() :
    pImpl(std::make_unique<SlidingWindowRealDFTImpl> ())
{
}

/// Copy constructor
SlidingWindowRealDFT::SlidingWindowRealDFT(const SlidingWindowRealDFT &swdft)
{
    *this = swdft;
}

/// Move constructor
SlidingWindowRealDFT::SlidingWindowRealDFT(
    SlidingWindowRealDFT &&swdft) noexcept
{
    *this = std::move(swdft);
}

/// Move assignment operator
SlidingWindowRealDFT&
SlidingWindowRealDFT::operator=(SlidingWindowRealDFT &&swdft) noexcept
{
    if (&swdft == this){return *this;}
    pImpl = std::move(swdft.pImpl);
    return *this;
}

/// Copy assignment operator
SlidingWindowRealDFT&
SlidingWindowRealDFT::operator=(const SlidingWindowRealDFT &swdft)
{
    if (&swdft == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<SlidingWindowRealDFTImpl> ();
    if (!swdft.pImpl->mInitialized){return *this;}
    // Call initialize so that we get a fresh FFTw context
    try
    {
        initialize(swdft.pImpl->mParameters);
    }
    catch (const std::exception &e)
    {
        clear();
        RTSEIS_THROW_RTE("%s", "Failed to initialize class");
    }
    // Copy the workspace
    if (pImpl->mPrecision == RTSeis::Precision::DOUBLE)
    {
        auto nbytes = static_cast<size_t> (pImpl->mInDataLength)
                     *sizeof(double);
        std::memcpy(pImpl->mInData64f, swdft.pImpl->mInData64f, nbytes);
        nbytes = static_cast<size_t> (pImpl->mOutDataLength)
                *sizeof(fftw_complex);
        std::memcpy(pImpl->mOutData64f, swdft.pImpl->mOutData64f, nbytes);
    }
    else
    {
        auto nbytes = static_cast<size_t> (pImpl->mInDataLength)
                     *sizeof(float);
        std::memcpy(pImpl->mInData32f, swdft.pImpl->mInData32f, nbytes);
        nbytes = static_cast<size_t> (pImpl->mOutDataLength)
                *sizeof(fftwf_complex);
        std::memcpy(pImpl->mOutData32f, swdft.pImpl->mOutData32f, nbytes);
    }
    return *this;
}

/// Destructor
SlidingWindowRealDFT::~SlidingWindowRealDFT() = default;

void SlidingWindowRealDFT::clear() noexcept
{
    pImpl->clear();
}

/// Initialize - this will fire up the FFTw engine
void SlidingWindowRealDFT::initialize(
    const SlidingWindowRealDFTParameters &parameters)
{
    clear();
    if (!parameters.isValid())
    {
        RTSEIS_THROW_IA("%s", "parameters are not valid");
    }
    // Extract the parameters 
    int nSamples = parameters.getNumberOfSamples();
    int nSamplesPerSegment = parameters.getWindowLength();
    int nSamplesInOverlap = parameters.getNumberOfSamplesInOverlap();
    int dftLength = parameters.getDFTLength();
    int windowLength = parameters.getWindowLength();
    bool luseWindow = false;
    if (parameters.getWindowType() != SlidingWindowWindowType::BOXCAR)
    {
        luseWindow = true; 
    }
    // Compute the sizes
    auto cols = static_cast<double> (nSamples - nSamplesInOverlap)
               /static_cast<double> (nSamplesPerSegment - nSamplesInOverlap);
    auto ncols = static_cast<int> (cols);
    pImpl->mSamples = nSamples;
    pImpl->mSamplesInOverlap = nSamplesInOverlap;
    pImpl->mSamplesPerSegment = nSamplesPerSegment;
    pImpl->mDFTLength = dftLength;
    pImpl->mNumberOfFrequencies = dftLength/2 + 1;
    pImpl->mNumberOfColumns = ncols;
    pImpl->mPrecision = parameters.getPrecision();
    pImpl->mDetrendType = parameters.getDetrendType();
    if (luseWindow)
    {
        auto window = parameters.getWindow();
        // Always copy the window for the copy constructor
        pImpl->mWindow64f = ippsMalloc_64f(windowLength);
        ippsCopy_64f(window.data(), pImpl->mWindow64f, windowLength);
        if (pImpl->mPrecision == RTSeis::Precision::FLOAT)
        {
            pImpl->mWindow32f = ippsMalloc_32f(windowLength);
            ippsConvert_64f32f(window.data(), pImpl->mWindow32f, windowLength);
        }
        pImpl->mApplyWindow = true;
    }
    // Initialize the Fourier transform
    constexpr int rank = 1;
    int nForward[1] = {pImpl->mDFTLength};
    int howMany = pImpl->mNumberOfColumns;
    constexpr int istride = 1;
    constexpr int ostride = 1;
    constexpr int inembed[1] = {0}; //NULL;
    constexpr int onembed[1] = {0}; //NULL;
    // Figure out padding for cache alignment
    if (pImpl->mPrecision == RTSeis::Precision::DOUBLE)
    {
        pImpl->mDataOffset = padLength64f(pImpl->mDFTLength, 64);
        pImpl->mFTOffset = padLength64f(pImpl->mNumberOfFrequencies, 64);
    }
    else
    {
        pImpl->mDataOffset = padLength32f(pImpl->mDFTLength, 64);
        pImpl->mFTOffset = padLength32f(pImpl->mNumberOfFrequencies, 64);
    }
    // Make the real-to-complex Fourier transform plans
    pImpl->mInDataLength  = pImpl->mDataOffset*pImpl->mNumberOfColumns;
    pImpl->mOutDataLength = pImpl->mFTOffset*pImpl->mNumberOfColumns;
    if (pImpl->mPrecision == RTSeis::Precision::DOUBLE)
    {
        auto nbytes = static_cast<size_t> (pImpl->mInDataLength)
                     *sizeof(double);
        pImpl->mInData64f = static_cast<double *> (fftw_malloc(nbytes));
        memset(pImpl->mInData64f, 0, nbytes);
        nbytes = static_cast<size_t> (pImpl->mOutDataLength)
                *sizeof(fftw_complex);
        pImpl->mOutData64f
            = reinterpret_cast<fftw_complex *> (fftw_malloc(nbytes));
        memset(pImpl->mOutData64f, 0, nbytes);
        pImpl->mDoublePlan = fftw_plan_many_dft_r2c(rank, nForward, howMany,
                                                    pImpl->mInData64f, inembed,
                                                    istride, pImpl->mDataOffset,
                                                    pImpl->mOutData64f, onembed,
                                                    ostride, pImpl->mFTOffset,
                                                    FFTW_PATIENT);
        pImpl->mHaveDoublePlan = true;
    }
    else
    {
        auto nbytes = static_cast<size_t> (pImpl->mInDataLength)
                     *sizeof(float);
        pImpl->mInData32f = static_cast<float *> (fftw_malloc(nbytes));
        memset(pImpl->mInData32f, 0, nbytes);
        nbytes = static_cast<size_t> (pImpl->mOutDataLength)
                *sizeof(fftwf_complex);
        pImpl->mOutData32f
            = reinterpret_cast<fftwf_complex *> (fftw_malloc(nbytes));
        memset(pImpl->mOutData32f, 0, nbytes);
        pImpl->mFloatPlan = fftwf_plan_many_dft_r2c(rank, nForward, howMany,
                                                    pImpl->mInData32f, inembed,
                                                    istride, pImpl->mDataOffset,
                                                    pImpl->mOutData32f, onembed,
                                                    ostride, pImpl->mFTOffset,
                                                    FFTW_PATIENT);
        pImpl->mHaveFloatPlan = true;
    }
    pImpl->mHaveTransform = false;
    pImpl->mInitialized = true;
}

/*
/// Initialize DFT engine
void SlidingWindowRealDFT::initialize(const int nSamples,
                                      const int nSamplesPerSegment,
                                      const int dftLength,
                                      const int nSamplesInOverlap,
                                      const int windowLength,
                                      const double window[],
                                      const SlidingWindowDetrendType detrendType,
                                      const RTSeis::Precision precision)
{
    clear();
    if (nSamples < 1)
    {
        RTSEIS_THROW_IA("nSamples = %d must be positive", nSamples);
    }
    if (nSamplesPerSegment < 1)
    {
        RTSEIS_THROW_IA("Samples per segement = %d must be positive",
                        nSamplesPerSegment);
    }
    if (dftLength < nSamplesPerSegment)
    {
        RTSEIS_THROW_IA("DFT length = %d must be at least = %d",
                        dftLength, nSamplesPerSegment);
    }
    if (nSamplesInOverlap < 0 || nSamplesInOverlap >= nSamplesPerSegment)
    {
        RTSEIS_THROW_IA("Overlap size = %d must be in range [0,%d]",
                        nSamplesInOverlap, nSamplesPerSegment-1);
    }
    bool luseWindow = false;
    if (windowLength > 0)
    {
        if (windowLength != nSamplesPerSegment)
        {
            RTSEIS_THROW_IA("Window length = %d must equal %d",
                            windowLength, nSamplesPerSegment);
        }
        if (window == nullptr){RTSEIS_THROW_IA("%s", "Window is NULL");}
        luseWindow = true;
    }
    // Compute the sizes
    auto cols = static_cast<double> (nSamples - nSamplesInOverlap)
               /static_cast<double> (nSamplesPerSegment - nSamplesInOverlap);
    auto ncols = static_cast<int> (cols);
    pImpl->mSamples = nSamples;
    pImpl->mSamplesInOverlap = nSamplesInOverlap;
    pImpl->mSamplesPerSegment = nSamplesPerSegment;
    pImpl->mDFTLength = dftLength;
    pImpl->mNumberOfFrequencies = dftLength/2 + 1;
    pImpl->mNumberOfColumns = ncols;
    pImpl->mPrecision = precision;
    pImpl->mDetrendType = detrendType;
    if (luseWindow)
    {
        // Always copy the window for the copy constructor
        pImpl->mWindow64f = ippsMalloc_64f(windowLength);
        ippsCopy_64f(window, pImpl->mWindow64f, windowLength);
        if (precision == RTSeis::Precision::FLOAT)
        {
            pImpl->mWindow32f = ippsMalloc_32f(windowLength);
            ippsConvert_64f32f(window, pImpl->mWindow32f, windowLength);
        }
        pImpl->mApplyWindow = true;
    }
    // Initialize the Fourier transform
    constexpr int rank = 1;
    int nForward[1] = {pImpl->mDFTLength};
    int howMany = pImpl->mNumberOfColumns;
    constexpr int istride = 1;
    constexpr int ostride = 1;
    constexpr int inembed[1] = {0}; //NULL;
    constexpr int onembed[1] = {0}; //NULL;
    // Figure out padding for cache alignment 
    if (precision == RTSeis::Precision::DOUBLE)
    {
        pImpl->mDataOffset = padLength64f(pImpl->mDFTLength, 64);
        pImpl->mFTOffset = padLength64f(pImpl->mNumberOfFrequencies, 64);
    }
    else
    {
        pImpl->mDataOffset = padLength32f(pImpl->mDFTLength, 64);
        pImpl->mFTOffset = padLength32f(pImpl->mNumberOfFrequencies, 64);
    }
    // Make the real-to-complex Fourier transform plans
    pImpl->mInDataLength  = pImpl->mDataOffset*pImpl->mNumberOfColumns;
    pImpl->mOutDataLength = pImpl->mFTOffset*pImpl->mNumberOfColumns;
    if (precision == RTSeis::Precision::DOUBLE)
    {
        auto nbytes = static_cast<size_t> (pImpl->mInDataLength)
                     *sizeof(double);
        pImpl->mInData64f = static_cast<double *> (fftw_malloc(nbytes));
        memset(pImpl->mInData64f, 0, nbytes);
        nbytes = static_cast<size_t> (pImpl->mOutDataLength)
                *sizeof(fftw_complex);
        pImpl->mOutData64f 
            = reinterpret_cast<fftw_complex *> (fftw_malloc(nbytes));
        memset(pImpl->mOutData64f, 0, nbytes);
        pImpl->mDoublePlan = fftw_plan_many_dft_r2c(rank, nForward, howMany,
                                                    pImpl->mInData64f, inembed,
                                                    istride, pImpl->mDataOffset,
                                                    pImpl->mOutData64f, onembed,
                                                    ostride, pImpl->mFTOffset,
                                                    FFTW_PATIENT);
        pImpl->mHaveDoublePlan = true;
    }
    else
    {
        auto nbytes = static_cast<size_t> (pImpl->mInDataLength)
                     *sizeof(float);
        pImpl->mInData32f = static_cast<float *> (fftw_malloc(nbytes));
        memset(pImpl->mInData32f, 0, nbytes);
        nbytes = static_cast<size_t> (pImpl->mOutDataLength)
                *sizeof(fftwf_complex);
        pImpl->mOutData32f 
            = reinterpret_cast<fftwf_complex *> (fftw_malloc(nbytes));
        memset(pImpl->mOutData32f, 0, nbytes);
        pImpl->mFloatPlan = fftwf_plan_many_dft_r2c(rank, nForward, howMany,
                                                    pImpl->mInData32f, inembed,
                                                    istride, pImpl->mDataOffset,
                                                    pImpl->mOutData32f, onembed,
                                                    ostride, pImpl->mFTOffset,
                                                    FFTW_PATIENT);
        pImpl->mHaveFloatPlan = true;
    }
    pImpl->mHaveTransform = false;
    pImpl->mInitialized = true;
}
*/

/// Gets number of frequencies
int SlidingWindowRealDFT::getNumberOfFrequencies() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");} 
    return pImpl->mNumberOfFrequencies;
}

/// Gets number of time samples
int SlidingWindowRealDFT::getNumberOfTransformWindows() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mNumberOfColumns;
}

/// Gets number of samples
int SlidingWindowRealDFT::getNumberOfSamples() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mSamples;
}

/// Determines if the class is intitialized
bool SlidingWindowRealDFT::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Gets the precision
RTSeis::Precision SlidingWindowRealDFT::getPrecision() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mPrecision;
}

/// Actually perform the transform
void SlidingWindowRealDFT::transform(const int nSamples, const double x[])
{
    pImpl->mHaveTransform = false;
    // Check the class is initialized and that the inputs are as expected
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (nSamples != getNumberOfSamples())
    {
        RTSEIS_THROW_IA("Number of samples = %d must equal %d",
                        nSamples, getNumberOfSamples());
    }
    if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
    // Handle the float case
    if (pImpl->mPrecision != RTSeis::Precision::DOUBLE)
    {
        RTSEIS_THROW_RTE("%s", "Float precision not yet implemented");
    }
    // Loop over the windows
    //#pragma omp parallel shared(inData, window) default(None)
    {
    auto *window = pImpl->mWindow64f;
    auto *inData = pImpl->mInData64f;
    int nDataOffset = pImpl->mDataOffset;
    int nPtsPerSeg = pImpl->mSamplesPerSegment;
    int nOverlap   = pImpl->mSamplesInOverlap;
    int shift = nPtsPerSeg - nOverlap;
    auto lwindow = pImpl->mApplyWindow;
    SlidingWindowDetrendType detrendType = pImpl->mDetrendType;
    //#pragma omp for
    for (auto icol=0; icol<pImpl->mNumberOfColumns; ++icol)
    {
        auto xIndex = icol*shift; // Extract x
        auto dataIndex = icol*nDataOffset;
        auto xptr = &x[xIndex];
        auto dptr = &inData[dataIndex];
#ifdef DEBUG
        assert(xIndex < nSamples);
#endif
        auto ncopy = std::min(nPtsPerSeg, nSamples - xIndex);
        // Zero and copy
        std::memset(dptr, 0, static_cast<size_t> (nDataOffset)*sizeof(double));
        std::copy(xptr, xptr + ncopy, dptr); //dptr, xptr, static_cast<size_t> (ncopy)*sizeof(double));
        // Demean?
        if (detrendType == SlidingWindowDetrendType::REMOVE_MEAN)
        {
            double mean;
            FilterImplementations::removeMean(ncopy, xptr, &dptr, &mean);
        }
        else if (detrendType == SlidingWindowDetrendType::REMOVE_TREND)
        {
            double intercept;
            double slope;
            FilterImplementations::removeTrend(ncopy, xptr, &dptr,
                                               &intercept, &slope);
        }
        // Window
        if (lwindow){ippsMul_64f_I(window, dptr, nPtsPerSeg);}
    }
    } // End parallel
    // Transform
    fftw_execute(pImpl->mDoublePlan);
    pImpl->mHaveTransform = true;
}

/// Returns a pointer to the transform in the i'th window
const std::complex<double> *
SlidingWindowRealDFT::getTransform64f(const int iWindow) const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (!pImpl->mHaveTransform)
    {
        RTSEIS_THROW_RTE("%s", "Transform not yet applied");
    }
    if (pImpl->mPrecision != RTSeis::Precision::DOUBLE)
    {
        RTSEIS_THROW_RTE("%s", "Precision is FLOAT - call getTransform32f");
    }
    if (iWindow < 0 || iWindow >= pImpl->mNumberOfColumns)
    {
        RTSEIS_THROW_IA("iWindow = %d must be in range [0,%d]",
                        iWindow, pImpl->mNumberOfColumns);
    }
    int indx = pImpl->mFTOffset*iWindow;
    auto ptr = reinterpret_cast<const std::complex<double> *>
               (pImpl->mOutData64f + indx);
    return ptr; 
} 

const std::complex<float> *
SlidingWindowRealDFT::getTransform32f(const int iWindow) const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    if (!pImpl->mHaveTransform)
    {
        RTSEIS_THROW_RTE("%s", "Transform not yet applied");
    }
    if (pImpl->mPrecision != RTSeis::Precision::FLOAT)
    {
        RTSEIS_THROW_RTE("%s", "Precision is DOUBLE - call getTransform64f");
    }
    if (iWindow < 0 || iWindow >= pImpl->mNumberOfColumns)
    {
        RTSEIS_THROW_IA("iWindow = %d must be in range [0,%d]",
                        iWindow, pImpl->mNumberOfColumns);
    }
    int indx = pImpl->mFTOffset*iWindow;
    auto ptr = reinterpret_cast<const std::complex<float> *>
               (pImpl->mOutData32f + indx);
    return ptr;
}
