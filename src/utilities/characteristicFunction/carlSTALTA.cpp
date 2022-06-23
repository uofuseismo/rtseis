#include <stdexcept>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#include <algorithm>
#include <ipps.h>
#include "rtseis/filterImplementations/firFilter.hpp"
#include "rtseis/utilities/characteristicFunction/carlSTALTA.hpp"
#include "rtseis/enums.hpp"

using namespace RTSeis::Utilities::CharacteristicFunction;

namespace 
{
template<RTSeis::ProcessingMode E, typename T>
void compute(const int n, const int mChunkSize,
             const T mRatio, const T mQuiet,
             const T *x, T *y,
             T *__attribute__((aligned(64))) xWork,
             T *__attribute__((aligned(64))) yWork,
             T *__attribute__((aligned(64))) staSignal,
             T *__attribute__((aligned(64))) ltaSignal,
             RTSeis::FilterImplementations::FIRFilter<E, T> &mSTAFilter,
             RTSeis::FilterImplementations::FIRFilter<E, T> &mLTAFilter,
             RTSeis::FilterImplementations::FIRFilter<E, T> &mSTARFilter,
             RTSeis::FilterImplementations::FIRFilter<E, T> &mLTARFilter)
{
    // Loop on chunks
    for (int i=0; i<n; i=i+mChunkSize)
    {
        auto nLocal = std::min(mChunkSize, n - i);
        auto xLoc = x + i;
        auto yLoc = y + i;
        // Compute the STA and LTA
        mSTAFilter.apply(nLocal, xLoc, &staSignal);
        mLTAFilter.apply(nLocal, xLoc, &ltaSignal);
        // Compute rectified short term average
        // as well as part of output: y[n] = |sta[n] - lta[n]|
        #pragma omp simd
        for (int j=0; j<nLocal; ++j)
        {
            xWork[j] = std::abs(xLoc[j] - ltaSignal[j]);
            yLoc[j] = std::abs(staSignal[j] - ltaSignal[j]);
        }
        // Compute the rectified short-term average
        mSTARFilter.apply(nLocal, xWork, &yWork);
        // Compute rectified long term average and part of output:
        // y[n] = rectifiedSTA[n] - |sta[n] - lta[n]|
        #pragma omp simd
        for (int j=0; j<nLocal; ++j)
        {
             xWork[j] = std::abs(xLoc[j] - ltaSignal[j]);
             yLoc[j] = yWork[j] - yLoc[j];
        }
        // Compute LTA
        mLTARFilter.apply(nLocal, xWork, &yWork);
        // Finish calculation:
        // y[n] = rectifiedSTA[n] - R*rectifiedLTA[n]
        //      - |sta[n] - lta[n] - Q
        #pragma omp simd
        for (int j=0; j<nLocal; ++j)
        {
            yLoc[j] = yLoc[j] - mRatio*yWork[j] - mQuiet;
        }
    }
}
}

template<class T, RTSeis::ProcessingMode E>
class CarlSTALTA<T, E>::CarlSTALTAImpl
{
public:
    /// C'tor
    CarlSTALTAImpl() = default;
    /// Copy c'tor
    CarlSTALTAImpl(const CarlSTALTAImpl &stalta)
    {
        *this = stalta;
    }
    /// Destructor
    ~CarlSTALTAImpl()
    {
        clear();
    }
    /// Copy
    CarlSTALTAImpl& operator=(const CarlSTALTAImpl &stalta)
    {
        mSTAFilter = stalta.mSTAFilter;
        mLTAFilter = stalta.mLTAFilter;
        mSTARFilter = stalta.mSTARFilter;
        mLTARFilter = stalta.mLTARFilter;
        mRatio = stalta.mRatio;
        mQuiet = stalta.mQuiet;
        mLTALength = stalta.mLTALength;
        mSTALength = stalta.mSTALength;
        mInitialized = stalta.mInitialized;;
        auto nbytes = static_cast<size_t> (mChunkSize)*sizeof(T);
        if (mDoublePrecision)
        {
            mXWork = ippsMalloc_64f(mChunkSize);
            mYWork = ippsMalloc_64f(mChunkSize);
            mSTA = ippsMalloc_64f(mChunkSize);
            mLTA = ippsMalloc_64f(mChunkSize);
        }
        else
        {
            mXWork = ippsMalloc_32f(mChunkSize);
            mYWork = ippsMalloc_32f(mChunkSize);
            mSTA = ippsMalloc_32f(mChunkSize);
            mLTA = ippsMalloc_32f(mChunkSize);
        }
        std::memcpy(mXWork, stalta.mXWork, nbytes);
        std::memcpy(mYWork, stalta.mYWork, nbytes); 
        std::memcpy(mSTA,   stalta.mSTA,   nbytes);
        std::memcpy(mLTA,   stalta.mLTA,   nbytes);
        return *this;
    }
    /// Resets the class
    void clear() noexcept
    {
        mSTAFilter.clear();
        mLTAFilter.clear();
        mSTARFilter.clear();
        mLTARFilter.clear();
        if (mXWork){ippsFree(mXWork);}
        if (mYWork){ippsFree(mYWork);}
        if (mSTA){ippsFree(mSTA);}
        if (mLTA){ippsFree(mLTA);}
        mXWork = nullptr;
        mYWork = nullptr;
        mSTA = nullptr;
        mLTA = nullptr;
        mRatio = 1;
        mQuiet = 0;
        mLTALength = 0;
        mSTALength = 0;
        mInitialized = false;
    } 
    /// Initialize
    void initialize(const int nSTA, const int nLTA,
                    const double ratio, const double quiet)
    {
        mSTALength = nSTA;
        mLTALength = nLTA;
        mRatio = static_cast<T> (ratio);
        mQuiet = static_cast<T> (quiet);
        std::vector<double> firCoeffs(mLTALength, 0.0);
        std::fill(firCoeffs.data(), firCoeffs.data() + mSTALength,
                  1./static_cast<double> (mSTALength));
        mSTAFilter.initialize(mSTALength, firCoeffs.data());
        mSTARFilter.initialize(mSTALength, firCoeffs.data());

        std::fill(firCoeffs.data(), firCoeffs.data() + mLTALength,
                  1./static_cast<double> (mLTALength));
        mLTAFilter.initialize(mLTALength, firCoeffs.data());
        mLTARFilter.initialize(mLTALength, firCoeffs.data());

        if (mDoublePrecision)
        {
            mXWork = ippsMalloc_64f(mChunkSize);
            mYWork = ippsMalloc_64f(mChunkSize);
            mSTA = ippsMalloc_64f(mChunkSize);
            mLTA = ippsMalloc_64f(mChunkSize);
        }
        else
        {
            mXWork = ippsMalloc_32f(mChunkSize);
            mYWork = ippsMalloc_32f(mChunkSize);
            mSTA = ippsMalloc_32f(mChunkSize);
            mLTA = ippsMalloc_32f(mChunkSize);
        }
        mInitialized = true;
    }
    /// Applies the filter to the signal (double precision implementation)
    void apply(const int n, const T x[], T y[])
    {
        auto xWork = static_cast<T *> (mXWork);
        auto yWork = static_cast<T *> (mYWork);
        auto staSignal = static_cast<T *> (mSTA);
        auto ltaSignal = static_cast<T *> (mLTA);
        compute(n, mChunkSize, mRatio, mQuiet, x, y,
                xWork, yWork, staSignal, ltaSignal,
                mSTAFilter, mLTAFilter, mSTARFilter, mLTARFilter);
        if (!mRealTime){resetInitialConditions();}
    }
    // Reset initial conditions
    void resetInitialConditions()
    {
        mSTAFilter.resetInitialConditions();
        mLTAFilter.resetInitialConditions();
        mSTARFilter.resetInitialConditions();
        mLTARFilter.resetInitialConditions(); 
    }
//private:
    /// Workspace for signals to filter.  This has dimension [mChunkSize].
    void *__attribute__((aligned(64))) mXWork = nullptr;
    /// Workspace for filtered signals.  This has dimension [mChunkSize].
    void *__attribute__((aligned(64))) mYWork = nullptr;
    /// Workspace for STA signal in block.  This has dimension [mChunkSize].
    void *__attribute__((aligned(64))) mSTA = nullptr;
    /// Workspace for LTA signal in block.  This has dimension [mChunkSize].
    void *__attribute__((aligned(64))) mLTA = nullptr;
    /// STA FIR filter implementation.
    RTSeis::FilterImplementations::FIRFilter<RTSeis::ProcessingMode::REAL_TIME, T> mSTAFilter;
    /// LTA FIR filter implementation.
    RTSeis::FilterImplementations::FIRFilter<RTSeis::ProcessingMode::REAL_TIME, T> mLTAFilter;
    /// Rectified STA filter implementation.
    RTSeis::FilterImplementations::FIRFilter<RTSeis::ProcessingMode::REAL_TIME, T> mSTARFilter;
    /// Rectified LTA filter implementation.
    RTSeis::FilterImplementations::FIRFilter<RTSeis::ProcessingMode::REAL_TIME, T> mLTARFilter;
    /// Ratio - scales the rectified LTA
    T mRatio = 1;
    /// Shift to increase sensitifity of characteristic function
    T mQuiet = 0;
    /// Number of samples in LTA
    int mLTALength = 0;
    /// Number of samples in STA
    int mSTALength = 0;
    /// Chunk size
    const int mChunkSize = 512; //2048;
    /// Flag indicating whether or not this is real-time
    const bool mRealTime
        = (E == RTSeis::ProcessingMode::REAL_TIME) ? true : false;
    const bool mDoublePrecision = (sizeof(T) == sizeof(double)) ? true : false;
    /// True indicates the class is initialized
    bool mInitialized = false;
};

/// C'tor
template<class T, RTSeis::ProcessingMode E>
CarlSTALTA<T, E>::CarlSTALTA() :
    pImpl(std::make_unique<CarlSTALTAImpl> ())
{
}

/// Copy c'tor
template<class T, RTSeis::ProcessingMode E>
CarlSTALTA<T, E>::CarlSTALTA(const CarlSTALTA &stalta)
{
    *this = stalta;
}

/// Move c'tor
template<class T, RTSeis::ProcessingMode E>
[[maybe_unused]]
CarlSTALTA<T, E>::CarlSTALTA(CarlSTALTA &&stalta) noexcept
{
    *this = std::move(stalta);
}

/// Copy assignment
template<class T, RTSeis::ProcessingMode E>
CarlSTALTA<T, E>& CarlSTALTA<T, E>::operator=(const CarlSTALTA<T, E> &stalta)
{
    if (&stalta == this){return *this;}
    pImpl = std::make_unique<CarlSTALTAImpl> (*stalta.pImpl);
    return *this;
}

/// Move assignment
template<class T, RTSeis::ProcessingMode E>
CarlSTALTA<T,E>& CarlSTALTA<T,E>::operator=(CarlSTALTA<T, E> &&stalta) noexcept
{
    if (&stalta == this){return *this;}
    pImpl = std::move(stalta.pImpl);
    return *this;
}

/// Destructor
template<class T, RTSeis::ProcessingMode E>
CarlSTALTA<T, E>::~CarlSTALTA() = default;

/// Reset the class
template<class T, RTSeis::ProcessingMode E>
void CarlSTALTA<T, E>::clear() noexcept
{
    pImpl->clear();
}

/// Initialize the class
template<class T, RTSeis::ProcessingMode E>
void CarlSTALTA<T, E>::initialize(
    const int nSTA, const int nLTA, const double ratio, const double quiet)
{
    if (nSTA < 2){throw std::invalid_argument("nSTA must be at least 2");}
    if (nLTA <= nSTA)
    {
        throw std::invalid_argument("nLTA = " + std::to_string(nLTA)
                                  + " must be greater than nSTA = "
                                  + std::to_string(nSTA)); 
    }
    pImpl->initialize(nSTA, nLTA, ratio, quiet);
}

/// Initialized?
template<class T, RTSeis::ProcessingMode E>
bool CarlSTALTA<T, E>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}


/// Apply the filter
template<class T, RTSeis::ProcessingMode E>
void CarlSTALTA<T, E>::apply(const int nx, const T x[], T *yIn[])
{
    if (nx < 1){return;}
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    auto y = *yIn;
    if (y == nullptr){throw std::invalid_argument("y is NULL");}
    pImpl->apply(nx, x, y);    
}

/// Gets the initial conditions length
template<class T, RTSeis::ProcessingMode E>
std::pair<int, int> CarlSTALTA<T, E>::getInitialConditionLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int nzSTA = pImpl->mSTAFilter.getInitialConditionLength();
    int nzLTA = pImpl->mLTAFilter.getInitialConditionLength();
    return std::pair(nzSTA, nzLTA);
}

/// Sets the initial conditions
template<class T, RTSeis::ProcessingMode E>
void CarlSTALTA<T, E>::setInitialConditions(const int nzSTA,
                                            const double zSTA[],
                                            const double zSTAR[],
                                            const int nzLTA,
                                            const double zLTA[],
                                            const double zLTAR[])
{
    auto [nzSTARef, nzLTARef] = getInitialConditionLength(); // Throws
    if (nzSTA != nzSTARef)
    {
        throw std::invalid_argument("STA initial condition length = " 
                                  + std::to_string(nzSTA) + " must equal "
                                  + std::to_string(nzSTARef));
     }
     if (nzLTA != nzLTARef)
     {
         throw std::invalid_argument("LTA initial condition length = "
                                   + std::to_string(nzLTA) + " must equal "
                                   + std::to_string(nzLTARef));
      }
      if (zSTA == nullptr){throw std::invalid_argument("zSTA is NULL");} 
      if (zSTAR == nullptr){throw std::invalid_argument("zSTAR is NULL");}
      if (zLTA == nullptr){throw std::invalid_argument("zLTA is NULL");}
      if (zLTAR == nullptr){throw std::invalid_argument("zLTAR is NULL");}
      pImpl->mSTAFilter.setInitialConditions(nzSTA,  zSTA);
      pImpl->mSTARFilter.setInitialConditions(nzSTA, zSTAR);
      pImpl->mLTAFilter.setInitialConditions(nzLTA,  zLTA);
      pImpl->mLTARFilter.setInitialConditions(nzLTA, zLTAR);
}

/// Reset the class
template<class T, RTSeis::ProcessingMode E>
void CarlSTALTA<T, E>::resetInitialConditions()
{
    pImpl->resetInitialConditions();
}

///--------------------------------------------------------------------------///
///                             Template Instantiation                       ///
///--------------------------------------------------------------------------///
template class RTSeis::Utilities::CharacteristicFunction::CarlSTALTA<double, RTSeis::ProcessingMode::REAL_TIME>; 
template class RTSeis::Utilities::CharacteristicFunction::CarlSTALTA<double, RTSeis::ProcessingMode::POST_PROCESSING>;
template class RTSeis::Utilities::CharacteristicFunction::CarlSTALTA<float, RTSeis::ProcessingMode::REAL_TIME>;
template class RTSeis::Utilities::CharacteristicFunction::CarlSTALTA<float, RTSeis::ProcessingMode::POST_PROCESSING>;
