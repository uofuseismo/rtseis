#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fftw/fftw3.h>
#include <ipps.h>
#include <rtseis/private/throw.hpp>
#include "rtseis/utilities/transforms/enums.hpp"

#include <complex>
#include "rtseis/enums.h"
#include <memory>


using namespace RTSeis::Utilities::Transforms;

/*!
 * @brief This is a utility class that slides a window over a real time series
 *        and computes the DFT in each window.  This is the basis of many 
 *        higher order spectral analysis methods such as the STFT, spectrogram,
 *        and Welch's method.  
 * @note This class is not intended for use by a large audience.  Because it is
 *       foundational to higher-order methods whose results are more desirable
 *       this interface emphasizes efficiency at the sake of clarity.
 * @author Ben Baker, University of Utah
 * @copyright Ben Baker distributed under the MIT license.
 * @date May 2019.
 */
class SlidingWindowRealDFT
{
public:
    /*!
     * @brief Default constructor.
     */
    SlidingWindowRealDFT();

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~SlidingWindowRealDFT();
    /*!
     * @brief Releases memory on the class.
     */
    void clear() noexcept;
    /*! @} */

    /*!
     * @brief Initializes the sliding-window real discrete Fourier Transform.
     * @param[in] nSamples            The number of samples in the signal.
     *                                This must be positive.
     * @param[in] nSamplesPerSegment  The number of samples in each DFT segment.
     * @param[in] dftLength           The length of the DFT.  If the number of
     *                                samples in each segment results in an 
     *                                inefficient transform length then this
     *                                will zero pad zero pad each window to 
     *                                the desired DFT length.  This must be
     *                                greater than or equal to
     *                                nSamplesPerSegment.
     * @param[in]
     * @param[in] windowLength        The length of the window.  If 0 then no
     *                                window function will be applied.
     *                                Otherwise, this must equal 
     *                                nSamplesPerSegment.
     * @param[in] window              The window to apply to each segment.
     *                                If windowLength is 0 then this is ignored.
     *                                Otherwise, it is an array of dimension
     *                                [windowLength].
     * @param[in] detrendType         Defines the detrend strategy to be 
     *                                applied to each signal prior to windowing.
     * @param[in] precision           The precision of the underlying DFT.
     * @throws std::invalid_argument if any parameters are incorrect.
     */
    void initialize(const int nSamples,
                    const int nSamplesPerSegment,
                    const int dftLength,
                    const int nSamplesInOverlap,
                    const int windowLength,
                    const double window[],
                    const SlidingWindowDetrendType detrendType,
                    const RTSeis::Precision precision = RTSeis::Precision::DOUBLE);

    /*!
     * @brief Flag indicating whether or not the class is initialized.
     * @result True indicates that the class is inititalized.
     */
    bool isInitialized() const noexcept;
    /*!
     * @brief Returns the number of frequencies.
     * @result The number of frequencies.
     * @throws std::runtime_error if the class is not intitialized.
     */
    int getNumberOfFrequencies() const;
    /*!
     * @brief Returns the number of time samples in the sliding window DFT.
     *        This is the number of columns in the output matrix.
     * @result The number of time samples or, equivalently, columsn in the
     *         output matrix.
     * @throws std::runtime_error if the class is not intitialized.
     */
    int getNumberOfTimeSamples() const;
    /*!
     * @brief Returns the expected number of samples in the time series.
     * @result The expected number of samples.
     * @throws std::runtime_error if the class is not initialized.
     */
    int getNumberOfSamples() const;

    /*!
     * @brief Computes the sliding window DFT of the real signal.
     * @param[in] nSamples  The number of samples in the signal.
     *                      This must match the result of 
     *                      \c getNumberOfSamples().
     * @param[in] x         The signal to transform.  This is an array whose
     *                      dimension is [nSamples].
     * @param[in] nWork     The workspace allocated to y.  This must be at
     *                      least \c getNumberOfFrequencies() x 
     *                            \c getNumberOfTimeSamples().
     * @param[out] y        This is a row-major matrix containing the
     *                      DFT's in each window.  Here, the rows correspond
     *                      to frequency and the columns correspond to time.
     *                      While this has dimension [nwork] only the first
     *                      [\c getNumberOfFrequencies() x
     *                       \c getNumberOfTimeSamples()] samples are valid.
     * @throws std::invalid_argument if any arguments are invalid.
     * @throws std::runtime_error if the class is not initalized.
     * @sa \c getNumberOfTimeSamples(), \c getNumberOfFrequencies()
     */
    void transform(const int nSamples, const double x[],
                   const int nWork, std::complex<double> y[]);
private:
    class SlidingWindowRealDFTImpl;
    std::unique_ptr<SlidingWindowRealDFTImpl> pImpl;
};

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
        if (mHaveDoublePlan){fftw_free(mDoublePlan);}
        if (mHaveRealPlan){fftwf_free(mFloatPlan);}
        if (mWindow64f != nullptr){ippsFree(mWindow64f);}
        if (mWindow32f != nullptr){ippsFree(mWindow32f);}
        if (mInData64f != nullptr){fftw_free(mInData64f);}
        if (mInData32f != nullptr){fftwf_free(mInData32f);}
        mWindow64f = nullptr;
        mWindow32f = nullptr;
        mInData64f = nullptr;
        mInData32f = nullptr;
        mNumberOfSamples = 0;
        mDFTLength = 0;
        mNumberOfFrequencies = 0;
        mNumberOfColumns = 0;
        mPrecision = RTSeis::Precision::DOUBLE;
        mHaveDoublePlan = false;
        mHaveRealPlan = false;
        mApplyWindow = false;
        mInitialized = false;
    }

//private:
    fftw_plan  mDoublePlan;
    fftwf_plan mFloatPlan;
    fftw_complex *mOutData64f = nullptr;
    fftwf_complex *mOutData32f = nullptr;
    double *mWindow64f = nullptr;
    double *mInData64f = nullptr;
    float *mWindow32f = nullptr;
    float *mInData32f = nullptr;
    int mNumberOfSamples = 0;
    int mDFTLength = 0;
    int mNumberOfFrequencies = 0;
    int mNumberOfColumns = 0;
    RTSeis::Precision mPrecision = RTSeis::Precision::DOUBLE;
    bool mHaveDoublePlan = false;
    bool mHaveRealPlan = false;
    bool mApplyWindow = false;
    bool mInitialized = false;
};

SlidingWindowRealDFT::SlidingWindowRealDFT() :
    pImpl(std::make_unique<SlidingWindowRealDFTImpl> ())
{
}

SlidingWindowRealDFT::~SlidingWindowRealDFT() = default;

void SlidingWindowRealDFT::clear() noexcept
{
    pImpl->clear();
}

/// Initialize DFT
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
    pImpl->mNumberOfSamples = nSamples;
    pImpl->mDFTLength = dftLength;
    pImpl->mNumberOfFrequencies = dftLength/2 + 1;
    pImpl->mNumberOfColumns = ncols;
    pImpl->mPrecision = precision;
    if (luseWindow)
    {
        if (precision == RTSeis::Precision::DOUBLE)
        {
            pImpl->mWindow64f = ippsMalloc_64f(windowLength);
            ippsCopy_64f(window, pImpl->mWindow64f, windowLength);
        }
        else
        {
            pImpl->mWindow32f = ippsMalloc_32f(windowLength);
            ippsConvert_64f32f(window, pImpl->mWindow32f, windowLength);
        }
    }
    // Initialize the Fourier transform
    constexpr int rank = 2; // 2D matrix
    // More natural to make the frequencies the columns and the times the rows
    int matrixDimensions[2] = {pImpl->mNumberOfColumns,
                               pImpl->mDFTLength};
    constexpr int howMany = 1;
    constexpr int *inembed = NULL;
    constexpr int *onembed = NULL;
    if (precision == RTSeis::Precision::DOUBLE)
    {
        auto nbytes = static_cast<size_t> (pImpl->mDFTLength)
                     *sizeof(double);
        pImpl->mInData64f = static_cast<double *> (fftw_malloc(nbytes));
        nbytes = static_cast<size_t> (pImpl->mNumberOfFrequencies)
                *sizeof(fftw_complex);
        //pImpl->mDoublePlan = fftw_plan_many_dft_r2c(rank, matrixDimensions, howMany, );
    }
    else
    {
        auto nbytes = static_cast<size_t> (pImpl->mDFTLength)
                     *sizeof(float);
        nbytes = static_cast<size_t> (pImpl->mNumberOfFrequencies)
                *sizeof(fftwf_complex);
        //pImpl->mFloatPlan = fftwf_plan_many_dft_r2c( );
    }
}

int SlidingWindowRealDFT::getNumberOfFrequencies() const
{
    if (!pImpl->mInitialized){RTSEIS_THROW_RTE("%s", "Class not initialized");} 
    return pImpl->mNumberOfFrequencies;
}

int SlidingWindowRealDFT::getNumberOfTimeSamples() const
{
    if (!pImpl->mInitialized){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mNumberOfColumns;
}

int SlidingWindowRealDFT::getNumberOfSamples() const
{
    if (!pImpl->mInitialized){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->mNumberOfSamples;
}

bool SlidingWindowRealDFT::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

