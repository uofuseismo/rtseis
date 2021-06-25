#include <cstdio>
#include <cstdlib>
#include <array>
#include <string>
#include <ipps.h>
#include "rtseis/transforms/firEnvelope.hpp"
#include "rtseis/utilities/filterDesign/fir.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include "rtseis/filterImplementations/firFilter.hpp"
#include "private/throw.hpp"

using namespace RTSeis::Utilities;
using namespace RTSeis::Transforms;

namespace
{
void ippsMean(const double x[], const int n, double *pMean)
{
    ippsMean_64f(x, n, pMean);
}
void ippsMean(const float x[], const int n, float *pMean)
{
    ippsMean_32f(x, n, pMean, ippAlgHintAccurate);
}
void ippsSubC(const double x[], const double pMean, double xPad[],
              const int n)
{
    ippsSubC_64f(x, pMean, xPad, n);
}
void ippsSubC(const float x[], const float pMean, float xPad[],
              const int n)
{
    ippsSubC_32f(x, pMean, xPad, n);
}
void ippsAddC_I(const double pMean, double y[], const int n)
{
    ippsAddC_64f_I(pMean, y, n);
}
void ippsAddC_I(const float pMean, float y[], const int n)
{
    ippsAddC_32f_I(pMean, y, n);
}
void ippsMagnitude(const double pSrcRe[], const double pSrcIm[],
                   double pDst[], const int len)
{
    ippsMagnitude_64f(pSrcRe, pSrcIm, pDst, len);
}
void ippsMagnitude(const float pSrcRe[], const float pSrcIm[],
                   float pDst[], const int len)
{
    ippsMagnitude_32f(pSrcRe, pSrcIm, pDst, len);
}
}

template<RTSeis::ProcessingMode E, class T>
class FIREnvelope<E, T>::FIREnvelopeImpl
{
public:
    FIREnvelopeImpl() = default;
    ~FIREnvelopeImpl() = default;
    FIREnvelopeImpl(const FIREnvelopeImpl &firEnvelope)
    {
        *this = firEnvelope;
    }
    FIREnvelopeImpl& operator=(const FIREnvelopeImpl &firEnvelope)
    {
        if (&firEnvelope == this){return *this;}
        mRealFIRFilter = firEnvelope.mRealFIRFilter;
        mImagFIRFilter = firEnvelope.mImagFIRFilter;
        mMean = firEnvelope.mMean;
        mNumberOfTaps = firEnvelope.mNumberOfTaps;
        mZeroPhase = firEnvelope.mZeroPhase;
        mType3 = firEnvelope.mType3;
        mHaveInitialCondition = firEnvelope.mHaveInitialCondition;
        mInitialized = firEnvelope.mInitialized;
        return *this;
    }

    FilterImplementations::FIRFilter<E, T> mRealFIRFilter;
    FilterImplementations::FIRFilter<E, T> mImagFIRFilter;
    double mMean = 0;
    int mNumberOfTaps = 0;
    bool mZeroPhase = true;
    bool mType3 = false;
    bool mHaveInitialCondition = false;
    bool mInitialized = false;
    const ProcessingMode mMode = E;
};

/// Constructor
template<RTSeis::ProcessingMode E, class T>
FIREnvelope<E, T>::FIREnvelope() :
    pImpl(std::make_unique<FIREnvelopeImpl> ())
{
}

/// Copy constructor
template<RTSeis::ProcessingMode E, class T>
FIREnvelope<E, T>::FIREnvelope(const FIREnvelope &firEnvelope)
{
    *this = firEnvelope;
}

/// Move constructor
template<RTSeis::ProcessingMode E, class T>
FIREnvelope<E, T>::FIREnvelope(FIREnvelope &&firEnvelope) noexcept
{
    *this = std::move(firEnvelope);
}

/// Copy assignment
template<RTSeis::ProcessingMode E, class T>
FIREnvelope<E, T>& FIREnvelope<E, T>::operator=(const FIREnvelope &firEnvelope)
{
    if (&firEnvelope == this){return *this;}
    pImpl = std::make_unique<FIREnvelopeImpl> (*firEnvelope.pImpl);
    return *this; 
}

/// Move assignment
template<RTSeis::ProcessingMode E, class T>
FIREnvelope<E, T>&
FIREnvelope<E, T>::operator=(FIREnvelope &&firEnvelope) noexcept
{
    if (&firEnvelope == this){return *this;}
    pImpl = std::move(firEnvelope.pImpl);
    return *this;
}

/// Destructor
template<RTSeis::ProcessingMode E, class T>
FIREnvelope<E, T>::~FIREnvelope() = default;

/// Clear the filter
template<RTSeis::ProcessingMode E, class T>
void FIREnvelope<E, T>::clear() noexcept
{
    pImpl->mRealFIRFilter.clear();
    pImpl->mImagFIRFilter.clear();
    pImpl->mMean = 0;
    pImpl->mNumberOfTaps = 0;
    pImpl->mZeroPhase = true;
    pImpl->mType3 = false;
    pImpl->mHaveInitialCondition = false;
    pImpl->mInitialized = false;
}

/// Check if initialized
template<RTSeis::ProcessingMode E, class T>
bool FIREnvelope<E, T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Initialize
template<RTSeis::ProcessingMode E, class T>
void FIREnvelope<E, T>::initialize(const int ntaps)
{
    clear();
    if (ntaps < 1)
    {
        throw std::invalid_argument("ntaps = " + std::to_string(ntaps)
                                  + " must be positive");
    }
    pImpl->mType3 = false;
    if (ntaps%2 == 1){pImpl->mType3 = true;}
    pImpl->mNumberOfTaps = ntaps;
    std::pair<FilterRepresentations::FIR, FilterRepresentations::FIR> zfir;
    constexpr double beta = 8;
    constexpr auto direct = FilterImplementations::FIRImplementation::DIRECT;
    // Create an FIR hilbert transform
    try
    {
        auto zfir = FilterDesign::FIR::HilbertTransformer(ntaps - 1, beta);
        auto rfir = zfir.first.getFilterTaps();
        pImpl->mRealFIRFilter.initialize(rfir.size(), rfir.data(), direct);
        auto cfir = zfir.second.getFilterTaps();
        pImpl->mImagFIRFilter.initialize(cfir.size(), cfir.data(), direct);
    }
    catch (const std::exception &e) 
    {
        clear();
        auto errmsg = "Failed to initialize Hilbert transformer: "
                    + std::string(e.what());
        throw std::runtime_error(errmsg);
    }
    pImpl->mInitialized = true;
}

/// Get the initial condition length
template<RTSeis::ProcessingMode E, class T>
int FIREnvelope<E, T>::getInitialConditionLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mImagFIRFilter.getInitialConditionLength();
}

/// Sets the initial conditions
template<RTSeis::ProcessingMode E, class T>
void FIREnvelope<E, T>::setInitialConditions(const int nz, const double zi[])
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    int nzRef = getInitialConditionLength();
    if (nz != nzRef)
    {
       throw std::invalid_argument("nz = " + std::to_string(nz)
                                 + " must equal " + std::to_string(nzRef));
    }
    pImpl->mRealFIRFilter.setInitialConditions(nz, zi);
    pImpl->mImagFIRFilter.setInitialConditions(nz, zi);
    pImpl->mHaveInitialCondition = true;
}

/// Resets the initial conditions
template<RTSeis::ProcessingMode E, class T>
void FIREnvelope<E, T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pImpl->mRealFIRFilter.resetInitialConditions();
    pImpl->mImagFIRFilter.resetInitialConditions();
}

/// Perform transform
template<RTSeis::ProcessingMode E, class T>
void FIREnvelope<E, T>::transform(const int n, const T x[], T *yIn[])
{
    pImpl->mMean = 0;
    if (n < 1){return;} // Nothing to do
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    auto y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Post-processing removes the phase shift
    if (pImpl->mMode == RTSeis::ProcessingMode::POST_PROCESSING)
    {
        // Compute the mean
        T pMean;
        ippsMean(x, n, &pMean);
        pImpl->mMean = pMean;
        // Remove the mean and pad out the signal
        // N.B. The group delay is actually + 1 but C wants to shift relative to
        // a base address so we subtract the one.  Hence, n/2 instead of n/2+1.
        int groupDelay = pImpl->mNumberOfTaps/2;
        int npad = n + groupDelay;
        auto xPad = reinterpret_cast<T *> (ippsMalloc_8u(npad*sizeof(T))); // ippsMalloc_64f(npad);
        ippsSubC(x, pMean, xPad, n);
        std::fill(xPad + n, xPad + n  + groupDelay, 0); // Post-pad with zeros
        //ippsZero_64f(&xPad[n], groupDelay); // Post-pad with zeros
        // Now apply the filter and compute absolute value - Type III
        if (pImpl->mType3)
        {
            T *yPadr = xPad;
            auto yPadi = reinterpret_cast<T *>(ippsMalloc_8u(npad*sizeof(T)));// ippsMalloc_64f(npad);
            pImpl->mImagFIRFilter.apply(npad, xPad, &yPadi); 
            ippsMagnitude(yPadr, &yPadi[groupDelay], y, n);
            ippsFree(yPadi);
        }
        else
        {
            auto yPadr = reinterpret_cast<T *>(ippsMalloc_8u(npad*sizeof(T))); //ippsMalloc_64f(npad);
            pImpl->mRealFIRFilter.apply(npad, xPad, &yPadr);
            auto yPadi = reinterpret_cast<T *>(ippsMalloc_8u(npad*sizeof(T))); //ippsMalloc_64f(npad);
            pImpl->mImagFIRFilter.apply(npad, xPad, &yPadi);
            ippsMagnitude(&yPadr[groupDelay], &yPadi[groupDelay], y, n);
            ippsFree(yPadr);
            ippsFree(yPadi);
        }
        ippsFree(xPad);
        // Reconstitute the mean
        ippsAddC_I(pMean, y, n);
    }
    else
    {
        // TODO - add a sparse FIR filter for the type III case
        constexpr int chunkSize = 1024;
        std::array<T, chunkSize> yrTemp;
        std::array<T, chunkSize> yiTemp;
        for (auto ic=0; ic<n; ic=ic+chunkSize)
        {
            auto npfilt = std::min(n - ic, chunkSize);
            auto yrTempDataPtr = yrTemp.data();
            auto yiTempDataPtr = yiTemp.data();
            pImpl->mRealFIRFilter.apply(npfilt, &x[ic], &yrTempDataPtr); //yrTemp.data());
            pImpl->mImagFIRFilter.apply(npfilt, &x[ic], &yiTempDataPtr); //yiTemp.data());
            ippsMagnitude(yrTemp.data(), yiTemp.data(), &y[ic], npfilt);
        }
    }
}

/*
template<>
void FIREnvelope<float>::transform(const int n, const float x[], float *yIn[])
{
    pImpl->mMean = 0;
    if (n < 1){return;} // Nothing to do
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    float *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Post-processing removes the phase shift
    if (pImpl->mMode == RTSeis::ProcessingMode::POST_PROCESSING)
    {
        // Compute the mean
        float pMean;
        ippsMean(x, n, &pMean, ippAlgHintAccurate);
        pImpl->mMean = pMean;
        // Remove the mean and pad out the signal
        // N.B. The group delay is actually + 1 but C wants to shift relative to
        // a base address so we subtract the one.  Hence, n/2 instead of n/2+1.
        int groupDelay = pImpl->mNumberOfTaps/2;
        int npad = n + groupDelay;
        float *xPad = ippsMalloc_32f(npad);
        ippsSubC(x, pMean, xPad, n);
        //ippsZero_32f(&xPad[n], groupDelay);
        std::fill(xPad + n, xPad + n  + groupDelay, 0); // Post-pad with zeros
        // Now apply the filter and compute absolute value - Type III
        if (pImpl->mType3)
        {
            float *yPadr = xPad;
            float *yPadi = ippsMalloc_32f(npad);
            pImpl->mImagFIRFilter.apply(npad, xPad, &yPadi);
            ippsMagnitude(yPadr, &yPadi[groupDelay], y, n);
            ippsFree(yPadi);
        }
        else
        {
            float *yPadr = ippsMalloc_32f(npad);
            pImpl->mRealFIRFilter.apply(npad, xPad, &yPadr);
            float *yPadi = ippsMalloc_32f(npad);
            pImpl->mImagFIRFilter.apply(npad, xPad, &yPadi);
            ippsMagnitude(&yPadr[groupDelay], &yPadi[groupDelay], y, n);
            ippsFree(yPadr);
            ippsFree(yPadi);
        }
        ippsFree(xPad);
        // Reconstitute the mean
        ippsAddC_I(pMean, y, n);
    }
    else
    {
        // TODO - add a sparse FIR filter for the type III case
        constexpr int chunkSize = 1024;
        std::array<float, chunkSize> yrTemp;
        std::array<float, chunkSize> yiTemp;
        for (auto ic=0; ic<n; ic=ic+chunkSize)
        {
            auto npfilt = std::min(n - ic, chunkSize);
            auto yrTempDataPtr = yrTemp.data();
            auto yiTempDataPtr = yiTemp.data();
            pImpl->mRealFIRFilter.apply(npfilt, &x[ic], &yrTempDataPtr); //yrTemp.data());
            pImpl->mImagFIRFilter.apply(npfilt, &x[ic], &yiTempDataPtr); //yiTemp.data());
            ippsMagnitude(yrTemp.data(), yiTemp.data(), &y[ic], npfilt);
        }
    }
}
*/

///--------------------------------------------------------------------------///
///                          Template Instantiation                          ///
///--------------------------------------------------------------------------///
template class RTSeis::Transforms::FIREnvelope<RTSeis::ProcessingMode::POST, double>;
template class RTSeis::Transforms::FIREnvelope<RTSeis::ProcessingMode::REAL_TIME, double>;
template class RTSeis::Transforms::FIREnvelope<RTSeis::ProcessingMode::POST, float>;
template class RTSeis::Transforms::FIREnvelope<RTSeis::ProcessingMode::REAL_TIME, float>;
