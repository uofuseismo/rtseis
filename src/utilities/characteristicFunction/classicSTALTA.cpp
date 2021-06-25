#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <vector>
#include <ipps.h>
#include "private/throw.hpp"
#include "rtseis/enums.hpp"
#include "rtseis/utilities/characteristicFunction/classicSTALTA.hpp"
#include "rtseis/filterImplementations/firFilter.hpp"

#define CHUNK_SIZE 1024

namespace RealTime = RTSeis::Utilities::CharacteristicFunction::RealTime;
namespace PostProcessing 
    = RTSeis::Utilities::CharacteristicFunction::PostProcessing;

namespace
{

void ippsSqr(const int n, const double x[], double y[])
{
    ippsSqr_64f(x, y, n);
}

void ippsSqr(const int n, const float x[], float y[])
{
    ippsSqr_32f(x, y, n);
}

int ippsDiv(const int n, const double x1[], const double x2[], double y[])
{
    return ippsDiv_64f(x1, x2, y, n);
}

int ippsDiv(const int n, const float x1[], const float x2[], float y[])
{
    return ippsDiv_32f(x1, x2, y, n);
}

void ippsSet(const int n, const double val, double x[])
{
    ippsSet_64f(val, x, n);
}

void ippsThresholdToZero(const int n, const double x[], double y[], double tol)
{
    ippsThreshold_LTAbsVal_64f(x, y, n, tol, 0);
}

void ippsThresholdToZero(const int n, const float x[], float y[], float tol)
{
    ippsThreshold_LTAbsVal_32f(x, y, n, tol, 0);
}

template<class T>
void ippsCalloc(const int n, T **result)
{
    int len = n*static_cast<int> (sizeof(T));
    auto temp = ippsMalloc_8u(len);
    ippsZero_8u(temp, len);
    *result = reinterpret_cast<T *> (temp);
}

}

template<RTSeis::ProcessingMode E, class T>
class RTSeis::Utilities::CharacteristicFunction::ClassicSTALTAImpl
{
public:
    /// Default c'tor
    ClassicSTALTAImpl() = default;
    /// Copy c'tor
    ClassicSTALTAImpl(const ClassicSTALTAImpl &stalta)
    {
        *this = stalta;
    }
    /// Destrutor
    ~ClassicSTALTAImpl()
    {
        clear();
    }
    /// Copy assignment
    ClassicSTALTAImpl& operator=(const ClassicSTALTAImpl &stalta)
    {
        if (&stalta == this){return *this;}
        mSTAFilter = stalta.mSTAFilter;
        mLTAFilter = stalta.mLTAFilter;
        mChunkSize = stalta.mChunkSize;
        mSta = stalta.mSta;
        mLta = stalta.mLta;
        mInitialized = stalta.mInitialized;
        if (mChunkSize > 0)
        {
            if (stalta.mX2 && mChunkSize > 0)
            {
                ippsCalloc(mChunkSize, &mX2);
                std::copy(stalta.mX2, stalta.mX2+mChunkSize, mX2);
            }
            if (stalta.mYNum && mChunkSize > 0)
            {
                ippsCalloc(mChunkSize, &mYNum);
                std::copy(stalta.mX2, stalta.mX2+mChunkSize, mYNum);
            }
            if (stalta.mYDen && mChunkSize > 0)
            {
                ippsCalloc(mChunkSize, &mYDen);
                std::copy(stalta.mX2, stalta.mX2+mChunkSize, mYDen);
            }
        }
        return *this;
    }
    /// Initialize
    void initialize(const int nSta, const int nLta, const int chunkSize)
    {
        // Set space
        mSta = nSta;
        mLta = nLta;
        mChunkSize = chunkSize; 
        ippsCalloc(mChunkSize, &mX2);
        ippsCalloc(mChunkSize, &mYNum);
        ippsCalloc(mChunkSize, &mYDen);
        // Create the numerator averaging FIR filter
        std::vector<double> filterCoeffs(std::max(nSta, nLta));
        double div =  1/static_cast<double> (nSta);
        ippsSet(nSta, div, filterCoeffs.data());
        mSTAFilter.initialize(nSta, filterCoeffs.data(),
           RTSeis::FilterImplementations::FIRImplementation::DIRECT);
        // Set the numerator's initial conditions to 0 
        ippsSet(nSta-1, 0, filterCoeffs.data());
        mSTAFilter.setInitialConditions(nSta-1, filterCoeffs.data());
 
        // Create the denominator averaging FIR filter
        div = 1/static_cast<double> (nLta);
        ippsSet(nLta, div, filterCoeffs.data());
        mLTAFilter.initialize(nLta, filterCoeffs.data(),
           RTSeis::FilterImplementations::FIRImplementation::DIRECT);

        // Set the denominator's initial conditions to a large number.  This
        // will result in a calculation on startup like 0/big which is 0.
        div = std::numeric_limits<T>::max()/static_cast<T> ((4*nLta));
        ippsSet(nLta-1, div, filterCoeffs.data());
        mLTAFilter.setInitialConditions(nLta-1, filterCoeffs.data());

        mInitialized = true;        
    }
    /// Applies the filter
    void apply(const int nx, const T x[], T y[])
    {
        for (int i=0; i<nx; i=i+mChunkSize)
        {
            auto nloc = std::min(mChunkSize, nx - i);
            // Square input signal 
            ippsSqr(nloc, &x[i], mX2);
            // Compute numerator average of squared input signal
            mSTAFilter.apply(nloc, mX2, &mYNum);
            // Compute denominator average of squared input signal
            mLTAFilter.apply(nloc, mX2, &mYDen);
            // Divide numerator by denominator
            auto status = ippsDiv(nloc, mYDen, mYNum, &y[i]);
            if (status == ippStsDivByZero)
            {
                //RTSEIS_WARNMSG("%s", "Division by zero detected");
                // If |yden| < epsilon then y = 0
                ippsThresholdToZero(nloc, mYDen, &y[i], 
                                    std::numeric_limits<T>::epsilon());
            }
        }
        if (mMode == RTSeis::ProcessingMode::POST_PROCESSING)
        {
            resetInitialConditions();
        }
    }
    /// Clears the module
    void clear()
    {
        if (mX2){ippsFree(mX2);}
        if (mYNum){ippsFree(mYNum);}
        if (mYDen){ippsFree(mYDen);}
        if (mSTAFilter.isInitialized()){mSTAFilter.clear();}
        if (mLTAFilter.isInitialized()){mLTAFilter.clear();}
        mX2 = nullptr;
        mYNum = nullptr;
        mYDen = nullptr;
        mChunkSize = CHUNK_SIZE;
        mSta = 0;
        mLta = 0;
        mInitialized = false;
        // Note: Do not touch mMode 
    }
    /// Numerator and Denominator initial condition length
    std::pair<int, int> getInitialConditionLength() const
    {
        int nNum = mSTAFilter.getInitialConditionLength();
        int nDen = mLTAFilter.getInitialConditionLength();
        return std::pair(nNum, nDen);
    }
    /// Resets the initial conditions
    void resetInitialConditions()
    {
        mSTAFilter.resetInitialConditions();
        mLTAFilter.resetInitialConditions();
    }
    /// Sets the initial conditions
    void setInitialConditions(const int nzNum, const double zNum[],
                              const int nzDen, const double zDen[])
    {
        resetInitialConditions();
        mSTAFilter.setInitialConditions(nzNum, zNum);
        mLTAFilter.setInitialConditions(nzDen, zDen);
    }
///private:
    RTSeis::FilterImplementations::FIRFilter<
        RTSeis::ProcessingMode::REAL_TIME, T> mSTAFilter;
    RTSeis::FilterImplementations::FIRFilter<
        RTSeis::ProcessingMode::REAL_TIME, T> mLTAFilter;
    T *mX2 = nullptr;
    T *mYNum = nullptr;
    T *mYDen = nullptr;
    int mChunkSize = CHUNK_SIZE;
    int mSta = 0;
    int mLta = 0;
    const RTSeis::ProcessingMode mMode = E;
    bool mInitialized = false; 
};

//----------------------------------------------------------------------------//
//                                  Post Processing                           //
//----------------------------------------------------------------------------//
/// Constructor
template<class T>
PostProcessing::ClassicSTALTA<T>::ClassicSTALTA() :
    pImpl(std::make_unique<
             ClassicSTALTAImpl<RTSeis::ProcessingMode::POST, T>> ())
{
}

/// Copy c'tor
template<class T>
PostProcessing::ClassicSTALTA<T>::ClassicSTALTA(const ClassicSTALTA &stalta)
{
    *this = stalta;
}

/// Move c'tor
template<class T>
PostProcessing::ClassicSTALTA<T>::ClassicSTALTA(ClassicSTALTA &&stalta) noexcept
{
    *this = std::move(stalta);
}

/// Copy assignment
template<class T>
PostProcessing::ClassicSTALTA<T>&
PostProcessing::ClassicSTALTA<T>::operator=(const ClassicSTALTA &stalta)
{
    if (&stalta == this){return *this;}
    pImpl = std::make_unique<ClassicSTALTAImpl
                             <RTSeis::ProcessingMode::POST, T>> (*stalta.pImpl);
    return *this;
}

/// Move assignment
template<class T>
PostProcessing::ClassicSTALTA<T>&
PostProcessing::ClassicSTALTA<T>::operator=(ClassicSTALTA &&stalta) noexcept
{
    if (&stalta == this){return *this;}
    pImpl = std::move(stalta.pImpl);
    return *this;
}

/// Destructor
template<class T>
PostProcessing::ClassicSTALTA<T>::~ClassicSTALTA()
{
    clear();
}

/// Clear
template<class T>
void PostProcessing::ClassicSTALTA<T>::clear() noexcept
{
    if (pImpl){pImpl->clear();}
}

/// Initialize
template<class T>
void PostProcessing::ClassicSTALTA<T>::initialize(
    const int nSTA, const int nLTA)
{
    clear();
    if (nSTA < 2){RTSEIS_THROW_IA("nSTA = %d must be positive", nSTA);}
    if (nLTA < nSTA)
    {
        RTSEIS_THROW_IA("nLTA = %d must be greater than %d", nLTA, nLTA);
    }
    pImpl->initialize(nSTA, nLTA, CHUNK_SIZE);
}

/// Determines if class is initialized
template<class T>
bool PostProcessing::ClassicSTALTA<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Gets the initial conditions length
template<class T>
std::pair<int, int>
PostProcessing::ClassicSTALTA<T>::getInitialConditionLength() const
{
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    return pImpl->getInitialConditionLength();
}

/// Sets the initial conditions
template<class T>
void PostProcessing::ClassicSTALTA<T>::setInitialConditions(
    const int nzNum, const double zNum[],
    const int nzDen, const double zDen[])
{
   auto nzLen = getInitialConditionLength(); // Throws
   if (nzNum != nzLen.first)
   {
       RTSEIS_THROW_IA("Number of numerator IC's = %d must equal %d",
                       nzNum, nzLen.first);
   }
   if (nzDen != nzLen.second)
   {
       RTSEIS_THROW_IA("Number of denominator IC's = %d must equal %d",
                       nzDen, nzLen.second);
   }
   // zNum must be at least length 1
   if (zNum == nullptr || zDen == nullptr)
   {
       if (zNum == nullptr){RTSEIS_THROW_IA("%s", "zNum is NULL");}
       RTSEIS_THROW_IA("%s", "zDen is NULL");
   }
   pImpl->setInitialConditions(nzNum, zNum, nzDen, zDen);
}

/// Applies the STA/LTA
template<class T>
void PostProcessing::ClassicSTALTA<T>::apply(
    const int nx, const T x[], T *yIn[])
{
    if (nx < 1){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}    
    auto y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y IS NULL");
    }
    pImpl->apply(nx, x, y);
}

//----------------------------------------------------------------------------//
//                                   Real Time                                //
//----------------------------------------------------------------------------//
/// Constructor
template<class T>
RealTime::ClassicSTALTA<T>::ClassicSTALTA() :
    pImpl(std::make_unique<ClassicSTALTAImpl<
            RTSeis::ProcessingMode::REAL_TIME, T>> ())
{
}

/// Copy c'tor
template<class T>
RealTime::ClassicSTALTA<T>::ClassicSTALTA(const ClassicSTALTA &stalta)
{
    *this = stalta;
}

/// Move c'tor
template<class T>
RealTime::ClassicSTALTA<T>::ClassicSTALTA(ClassicSTALTA &&stalta) noexcept
{
    *this = std::move(stalta);
}

/// Copy assignment
template<class T>
RealTime::ClassicSTALTA<T>&
RealTime::ClassicSTALTA<T>::operator=(const ClassicSTALTA &stalta)
{
    if (&stalta == this){return *this;}
    pImpl = std::make_unique<ClassicSTALTAImpl<
                 RTSeis::ProcessingMode::REAL_TIME, T>> (*stalta.pImpl);
    return *this;
}

/// Move assignment
template<class T>
RealTime::ClassicSTALTA<T>&
RealTime::ClassicSTALTA<T>::operator=(ClassicSTALTA &&stalta) noexcept
{
    if (&stalta == this){return *this;}
    pImpl = std::move(stalta.pImpl);
    return *this;
}

/// Destructor
template<class T>
RealTime::ClassicSTALTA<T>::~ClassicSTALTA()
{
    clear();
}

/// Clear
template<class T>
void RealTime::ClassicSTALTA<T>::clear() noexcept
{
    if (pImpl){pImpl->clear();}
}

/// Initialize
template<class T>
void RealTime::ClassicSTALTA<T>::initialize(
    const int nSTA, const int nLTA)
{
    clear();
    if (nSTA < 2){RTSEIS_THROW_IA("nSTA = %d must be positive", nSTA);}
    if (nLTA < nSTA)
    {
        RTSEIS_THROW_IA("nLTA = %d must be greater than %d", nLTA, nLTA);
    }
    pImpl->initialize(nSTA, nLTA, CHUNK_SIZE);
}

/// Determines if class is initialized
template<class T>
bool RealTime::ClassicSTALTA<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Gets the initial conditions length
template<class T>
std::pair<int, int>
RealTime::ClassicSTALTA<T>::getInitialConditionLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->getInitialConditionLength();
}

/// Sets the initial conditions
template<class T>
void RealTime::ClassicSTALTA<T>::setInitialConditions(
    const int nzNum, const double zNum[],
    const int nzDen, const double zDen[])
{
   auto nzLen = getInitialConditionLength(); // Throws
   if (nzNum != nzLen.first)
   {
       RTSEIS_THROW_IA("Number of numerator IC's = %d must equal %d",
                       nzNum, nzLen.first);
   }
   if (nzDen != nzLen.second)
   {
       RTSEIS_THROW_IA("Number of denominator IC's = %d must equal %d",
                       nzDen, nzLen.second);
   }
   // zNum must be at least length 1
   if (zNum == nullptr || zDen == nullptr)
   {
       if (zNum == nullptr){RTSEIS_THROW_IA("%s", "zNum is NULL");}
       RTSEIS_THROW_IA("%s", "zDen is NULL");
   }
   pImpl->setInitialConditions(nzNum, zNum, nzDen, zDen);
}

/// Applies the STA/LTA
template<class T>
void RealTime::ClassicSTALTA<T>::apply(
    const int nx, const T x[], T *yIn[])
{
    if (nx < 1){return;}
    if (!isInitialized()){RTSEIS_THROW_RTE("%s", "Class not initialized");}
    auto y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){RTSEIS_THROW_IA("%s", "x is NULL");}
        RTSEIS_THROW_IA("%s", "y IS NULL");
    }
    pImpl->apply(nx, x, y);
}

/// Resets the initial conditions
template<class T>
void RealTime::ClassicSTALTA<T>::resetInitialConditions()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pImpl->resetInitialConditions();
}

//----------------------------------------------------------------------------//
//                          Template instantiation                            //
//----------------------------------------------------------------------------//
template class
RTSeis::Utilities::CharacteristicFunction::PostProcessing::ClassicSTALTA<double>;
template class
RTSeis::Utilities::CharacteristicFunction::PostProcessing::ClassicSTALTA<float>;
template class
RTSeis::Utilities::CharacteristicFunction::RealTime::ClassicSTALTA<double>;
template class
RTSeis::Utilities::CharacteristicFunction::RealTime::ClassicSTALTA<float>;
