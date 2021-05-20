#include "rtseis/utilities/transforms/spectrogram.hpp"
#include "rtseis/utilities/transforms/slidingWindowRealDFT.hpp"
#include "rtseis/utilities/transforms/slidingWindowRealDFTParameters.hpp"
#include "rtseis/utilities/transforms/utilities.hpp"

using namespace RTSeis::Utilities::Transforms;

template<class T>
class Spectrogram<T>::SpectrogramImpl
{
public:
    class SlidingWindowRealDFT<T> mSlidingWindowRealDFT;
    class SlidingWindowRealDFTParameters mParameters;
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



///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class RTSeis::Utilities::Transforms::Spectrogram<double>;
template class RTSeis::Utilities::Transforms::Spectrogram<float>;

