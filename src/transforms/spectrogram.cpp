#include "rtseis/transforms/spectrogram.hpp"
#include "rtseis/transforms/slidingWindowRealDFT.hpp"
#include "rtseis/transforms/slidingWindowRealDFTParameters.hpp"
#include "rtseis/transforms/utilities.hpp"

using namespace RTSeis::Transforms;

template<class T>
class Spectrogram<T>::SpectrogramImpl
{
public:
    SlidingWindowRealDFT<T> mSlidingWindowRealDFT;
    SlidingWindowRealDFTParameters mParameters;
    std::vector<T> mAmplitude;
    std::vector<T> mPhase;
    double mSamplingRate = 1;
    bool mHavePhase = false; 
    bool mHaveAmplitude = false;
    bool mHaveTransform = false;
    bool mInitialized = false;
};

/// Reset the class
template<class T>
void Spectrogram<T>::clear() noexcept
{
    pImpl->mSlidingWindowRealDFT.clear();
    pImpl->mParameters.clear();
    pImpl->mAmplitude.clear();
    pImpl->mPhase.clear();
    pImpl->mSamplingRate = 1;
    pImpl->mHavePhase = false;
    pImpl->mHaveAmplitude = false;
    pImpl->mHaveTransform = false;
    pImpl->mInitialized = false;
}

/// C'tor
template<class T>
Spectrogram<T>::Spectrogram() :
    pImpl(std::make_unique<SpectrogramImpl> ())
{
}

/// Copy c'tor
template<class T>
Spectrogram<T>::Spectrogram(const Spectrogram &spectrogram)
{
    *this = spectrogram;
}

/// Move c'tor
template<class T>
Spectrogram<T>::Spectrogram(Spectrogram &&spectrogram) noexcept
{
    *this = std::move(spectrogram);
}

/// Copy assignment 
template<class T>
Spectrogram<T>& Spectrogram<T>::operator=(const Spectrogram &spectrogram)
{
    if (&spectrogram == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<SpectrogramImpl> (*spectrogram.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Spectrogram<T>& Spectrogram<T>::operator=(Spectrogram &&spectrogram) noexcept
{
    if (&spectrogram == this){return *this;}
    pImpl = std::move(spectrogram.pImpl);
    return *this;
}

/// Destructor
template<class T>
Spectrogram<T>::~Spectrogram() = default;

/// Initialize
template<class T>
void Spectrogram<T>::initialize(
    const SlidingWindowRealDFTParameters &parameters,
    const double samplingRate)
{
    clear();
    if (samplingRate <= 0)
    {
        throw std::invalid_argument("samplingRate = "
                                   + std::to_string(samplingRate)
                                   + " must be positive");
    }
    if (!parameters.isValid())
    {
        throw std::invalid_argument("parameters class is invalid");
    }
    // Initialize the sliding window DFT
    pImpl->mSamplingRate = samplingRate;
    try
    {
        pImpl->mParameters = parameters;
        pImpl->mSlidingWindowRealDFT.initialize(pImpl->mParameters);
    }
    catch (const std::exception &e)
    {
        clear();
        throw std::runtime_error("Failed to initialize DFT");
    }
    pImpl->mInitialized = true;
}

/// Actually perform the transform
template<class T>
void Spectrogram<T>::transform(const int nSamples, const T x[])
{
    pImpl->mHaveTransform = false;
    pImpl->mHavePhase = false;
    pImpl->mHaveAmplitude = false;
    // Check the class is initialized and that the inputs are as expected
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (nSamples != getNumberOfSamples())
    {
        throw std::invalid_argument("Number of samples = "
                                  + std::to_string(nSamples) + " must equal "
                                  + std::to_string(getNumberOfSamples()));
    }
    if (x == nullptr){throw std::invalid_argument("x is NULL");}
    pImpl->mSlidingWindowRealDFT.transform(nSamples, x);
    pImpl->mHaveTransform = true;
}

/// Get frequencies
template<class T>
std::vector<T> Spectrogram<T>::getFrequencies() const
{
    return pImpl->mSlidingWindowRealDFT.getFrequencies(pImpl->mSamplingRate);
}

/// Get frequencies
template<class T>
void Spectrogram<T>::getFrequencies(const int nFrequencies, T *freqsIn[]) const
{
    pImpl->mSlidingWindowRealDFT.getFrequencies(pImpl->mSamplingRate,
                                                nFrequencies, freqsIn);
}

/// Get the number of transform windows
template<class T>
int Spectrogram<T>::getNumberOfTransformWindows() const
{
    return pImpl->mSlidingWindowRealDFT.getNumberOfTransformWindows();
}

/// Initialized?
template<class T>
bool Spectrogram<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Have transform?
template<class T>
bool Spectrogram<T>::haveTransform() const noexcept
{
    return pImpl->mHaveTransform;
}

/// Have number of frequncies?
template<class T>
int Spectrogram<T>::getNumberOfFrequencies() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mSlidingWindowRealDFT.getNumberOfFrequencies();
}

/// Get number of samples
template<class T>
int Spectrogram<T>::getNumberOfSamples() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mSlidingWindowRealDFT.getNumberOfSamples();
}

/// Get the time domain windows
template<class T>
std::vector<T> Spectrogram<T>::getTimeWindows() const
{
    if (!isInitialized()){throw std::invalid_argument("Class not initialized");}
    return pImpl->mSlidingWindowRealDFT.getTimeWindows(pImpl->mSamplingRate);
}

/// Get the time domain windows
template<class T>
void Spectrogram<T>::getTimeWindows(const int nSamples, T *times[]) const
{
    if (!isInitialized()){throw std::invalid_argument("Class not initialized");}
    pImpl->mSlidingWindowRealDFT.getTimeWindows(pImpl->mSamplingRate,
                                                nSamples, times);
}

/// Gets the amplitude spectrum
template<class T>
std::vector<T> Spectrogram<T>::getAmplitude() const
{
    if (!haveTransform())
    {
        throw std::runtime_error("Transform not yet computed");
    }
    auto nFrequencies = getNumberOfFrequencies();
    auto nWindows = getNumberOfTransformWindows();
    std::vector<T> result(nFrequencies*nWindows);
    // Compute the amplitude spectrogram (or just get a pointer)
    auto ampPtr = getAmplitudePointer(); 
    std::copy(ampPtr, ampPtr + result.size(), result.data());
    return result;    
}

/// Gets the amplitude spectrum
template<class T>
const T *Spectrogram<T>::getAmplitudePointer() const
{
    if (!haveTransform())
    {
        throw std::runtime_error("Transform not yet computed");
    }
    if (!pImpl->mHaveAmplitude)
    {
        auto nFrequencies = getNumberOfFrequencies();
        auto nWindows = getNumberOfTransformWindows();
        pImpl->mAmplitude.resize(nWindows*nFrequencies, 0);
        for (int iw = 0; iw < nWindows; ++iw)
        {
            auto spectraPtr = pImpl->mSlidingWindowRealDFT.getTransform(iw);
            auto ampPtr = pImpl->mAmplitude.data() + iw*nFrequencies;
            DFTUtilities::magnitude(nFrequencies, spectraPtr, &ampPtr);
        }
        pImpl->mHaveAmplitude = true;
    }
    return pImpl->mAmplitude.data();
}

/// Gets the phase spectrum
template<class T>
std::vector<T> Spectrogram<T>::getPhase() const
{
    if (!haveTransform())
    {
        throw std::runtime_error("Transform not yet computed");
    }
    auto nFrequencies = getNumberOfFrequencies();
    auto nWindows = getNumberOfTransformWindows();
    std::vector<T> result(nFrequencies*nWindows);
    // Compute the amplitude spectrogram (or just get a pointer)
    auto phasePtr = getPhasePointer(); 
    std::copy(phasePtr, phasePtr + result.size(), result.data());
    return result;    
}

/// Gets the phase spectrum
template<class T>
const T *Spectrogram<T>::getPhasePointer() const
{
    if (!haveTransform())
    {
        throw std::runtime_error("Transform not yet computed");
    }
    if (!pImpl->mHavePhase)
    {
        constexpr bool wantDegrees = false;
        auto nFrequencies = getNumberOfFrequencies();
        auto nWindows = getNumberOfTransformWindows();
        pImpl->mPhase.resize(nWindows*nFrequencies, 0); 
        for (int iw = 0; iw < nWindows; ++iw)
        {
            auto spectraPtr = pImpl->mSlidingWindowRealDFT.getTransform(iw);
            auto phasePtr = pImpl->mPhase.data() + iw*nFrequencies;
            DFTUtilities::phase(nFrequencies, spectraPtr,
                                &phasePtr, wantDegrees);
        }
        pImpl->mHavePhase = true;
    }
    return pImpl->mPhase.data();
}


///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class RTSeis::Transforms::Spectrogram<double>;
template class RTSeis::Transforms::Spectrogram<float>;

