#include <cstdio>
#include <cstdlib>
#include <complex>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/transforms/envelope.hpp"
#include "rtseis/transforms/hilbert.hpp"

using namespace RTSeis::Transforms;

template<>
class Envelope<double>::EnvelopeImpl
{
public: 
    /// Constructor
    EnvelopeImpl() = default;
    /// Destructor
    ~EnvelopeImpl()
    {
        clear();
    }
    /// Copy constructor
    EnvelopeImpl(const EnvelopeImpl &envelope)
    {
        *this = envelope;
    }
    /// Clears the memory
    void clear() noexcept
    {
       mHilbert.clear();
       if (mHilb64fc){ippsFree(mHilb64fc);}
       mHilb64fc = nullptr;
       mEnvelopeLength = 0;
       mMean = 0;
       mInitialized = false;
    }
    /// Initializes the class
    void initialize(const int n)
    {
        if (n > 1){mHilbert.initialize(n);}
        mHilb64fc = ippsMalloc_64fc(n);
        mEnvelopeLength = n;
        mInitialized = true; 
    }
    /// Copy assignment operator
    EnvelopeImpl& operator=(const EnvelopeImpl &envelope)
    {
        if (&envelope == this){return *this;}
        clear();
        mHilbert = envelope.mHilbert;
        mEnvelopeLength = envelope.mEnvelopeLength;
        mMean = envelope.mMean;
        mInitialized = envelope.mInitialized;
        if (mEnvelopeLength > 0)
        {
            mHilb64fc = ippsMalloc_64fc(mEnvelopeLength);
            ippsCopy_64fc(envelope.mHilb64fc, mHilb64fc, mEnvelopeLength);
        }
        return *this;
    }
    void transform(const int n, const double x[], double y[])
    {
        // Special case
        if (n == 1)
        {
            mMean = x[0];
            y[0] = x[0];
            return;
        }
        // Remove the mean
        ippsMean_64f(x, n, &mMean);
        ippsSubC_64f(x, mMean, y, n);
        // Compute the analytic signal of the data.
        auto yhilb = reinterpret_cast<std::complex<double> *> (mHilb64fc);
        mHilbert.transform(n, y, &yhilb);
        // Take the absolute value to obtain the envelope.
        // |x_r + i x_h| = x_r^2 + x_h^2 where x_r is the input signal
        // and x_h is the Hilbert transform of the data.
        ippsMagnitude_64fc(mHilb64fc, y, n);
        // Restore the mean
        ippsAddC_64f_I(mMean, y, n);
    }
//private::
    Transforms::Hilbert<double> mHilbert;
    Ipp64fc *mHilb64fc = nullptr;
    double mMean = 0;
    int mEnvelopeLength = 0;
    bool mInitialized = false;
};

template<>
class Envelope<float>::EnvelopeImpl
{
public:
    /// Constructor
    EnvelopeImpl() = default;
    /// Destructor
    ~EnvelopeImpl()
    {
        clear();
    }
    /// Copy constructor
    EnvelopeImpl(const EnvelopeImpl &envelope)
    {
        *this = envelope;
    }
    /// Clears the memory
    void clear() noexcept
    {
        mHilbert.clear();
        if (mHilb32fc){ippsFree(mHilb32fc);}
        mHilb32fc = nullptr;
        mEnvelopeLength = 0;
        mMean = 0;
        mInitialized = false;
    }
    /// Initializes the class
    void initialize(const int n)
    {
        if (n > 1){mHilbert.initialize(n);}
        mHilb32fc = ippsMalloc_32fc(n);
        mEnvelopeLength = n;
        mInitialized = true;
    }
    /// Copy assignment operator
    EnvelopeImpl& operator=(const EnvelopeImpl &envelope)
    {
        if (&envelope == this){return *this;}
        clear();
        mHilbert = envelope.mHilbert;
        mEnvelopeLength = envelope.mEnvelopeLength;
        mMean = envelope.mMean;
        mInitialized = envelope.mInitialized;
        if (mEnvelopeLength > 0)
        {
            mHilb32fc = ippsMalloc_32fc(mEnvelopeLength);
            ippsCopy_32fc(envelope.mHilb32fc, mHilb32fc, mEnvelopeLength);
        }
        return *this;
    }
    void transform(const int n, const float x[], float y[])
    {
        // Special case
        if (n == 1)
        {
            mMean = x[0];
            y[0] = x[0];
            return;
        }
        // Remove the mean
        float pMean32f;
        ippsMean_32f(x, n, &pMean32f, ippAlgHintAccurate);
        ippsSubC_32f(x, pMean32f, y, n);
        mMean = static_cast<double> (pMean32f);
        // Compute the analytic signal of the data.
        auto yhilb = reinterpret_cast<std::complex<float> *> (mHilb32fc);
        mHilbert.transform(n, y, &yhilb);
        // Take the absolute value to obtain the envelope.
        // |x_r + i x_h| = x_r^2 + x_h^2 where x_r is the input signal
        // and x_h is the Hilbert transform of the data.
        ippsMagnitude_32fc(mHilb32fc, y, n);
        // Restore the mean
        ippsAddC_32f_I(pMean32f, y, n);
    }
//private
    Transforms::Hilbert<float> mHilbert;
    Ipp32fc *mHilb32fc = nullptr;
    double mMean = 0;
    int mEnvelopeLength = 0;
    bool mInitialized = false;
};

//============================================================================//
/// Constructor
template<class T>
Envelope<T>::Envelope() :
    pImpl(std::make_unique<EnvelopeImpl> ())
{
}

/// Copy constructor
template<class T>
Envelope<T>::Envelope(const Envelope &envelope)
{
    *this = envelope;
} 

/// Move constructor
template<class T>
Envelope<T>::Envelope(Envelope &&envelope) noexcept
{
    *this = std::move(envelope);
}

/// Copy assignment operator
template<class T>
Envelope<T>& Envelope<T>::operator=(const Envelope &envelope)
{
    if (&envelope == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::make_unique<EnvelopeImpl> (*envelope.pImpl);
    return *this;
}

/// Move assignment oeprator
template<class T>
Envelope<T>& Envelope<T>::operator=(Envelope &&envelope) noexcept
{
    if (&envelope == this){return *this;}
    if (pImpl){pImpl.reset();}
    pImpl = std::move(envelope.pImpl);
    return *this;
}
/// Destructor
template<class T>
Envelope<T>::~Envelope() = default;

/// Clears the module
template<class T>
void Envelope<T>::clear() noexcept
{
    pImpl->clear();
}

/// Checks if the module is initialized
template<class T>
bool Envelope<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Gets the envelope length
template<class T>
int Envelope<T>::getTransformLength() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mEnvelopeLength;
}

/// Initializes the class
template<class T>
void Envelope<T>::initialize(const int n)
{
    clear();
    if (n < 1)
    {
        throw std::invalid_argument("n = " + std::to_string(n)
                                  + " must be positive");
    }
    pImpl->initialize(n);
}

/// Compute the upper and lower envelope
template<>
void Envelope<double>::transform(const int n, const double x[],
                                 double *yupperIn[], double *ylowerIn[])
{
    if (pImpl){pImpl->mMean = 0;}
    double *ylower = *ylowerIn;
    double *yupper = *yupperIn;
    if (ylower == nullptr){throw std::invalid_argument("ylower is NULL");}
    if (yupper == nullptr){throw std::invalid_argument("yupper is NULL");}
    // Compute the upper envelope
    transform(n, x, &yupper); // will throw
    // Special case
    if (n == 1)
    {
        ylower[0] = yupper[0];
        return;
    } 
    // Compute the lower envelope from the upper envelope.
    //  (1) yupper = |Hilbert| + mean
    //  (2) ylower = |Hilbert| - mean
    // Expressing |Hilbert| as a function of yupper in Eqn 2
    //   |Hilbert| = yupper - mean
    // And substituting into the second equation
    //   ylower = yupper - mean - mean 
    //          = yupper - 2*mean
    ippsSubCRev_64f(yupper, 2*pImpl->mMean, ylower, n);
}

template<>
void Envelope<float>::transform(const int n, const float x[],
                                float *yupperIn[], float *ylowerIn[])
{
    if (pImpl){pImpl->mMean = 0;}
    float *ylower = *ylowerIn;
    float *yupper = *yupperIn;
    if (ylower == nullptr){throw std::invalid_argument("ylower is NULL");}
    if (yupper == nullptr){throw std::invalid_argument("yupper is NULL");}
    // Compute the upper envelope
    transform(n, x, &yupper); // will throw
    // Special case
    if (n == 1)
    {
        ylower[0] = yupper[0];
        return;
    }
    // Compute the lower envelope from the upper envelope.
    //  (1) yupper = |Hilbert| + mean
    //  (2) ylower = |Hilbert| - mean
    // Expressing |Hilbert| as a function of yupper in Eqn 2
    //   |Hilbert| = yupper - mean
    // And substituting into the second equation
    //   ylower = yupper - mean - mean
    //          = yupper - 2*mean
    auto twoMean = static_cast<float> (2*pImpl->mMean);
    ippsSubCRev_32f(yupper, twoMean, ylower, n);
}

/// Compute the envelope
template<class T>
void Envelope<T>::transform(const int n, const T x[], T *yIn[])
{
    pImpl->mMean = 0;
    if (!isInitialized()){throw std::invalid_argument("Class not initialized");}
    if (n != pImpl->mEnvelopeLength)
    {
        throw std::invalid_argument("n = " + std::to_string(n)
                                  + " must equal "
                                  + std::to_string(pImpl->mEnvelopeLength));
    }
    T *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    pImpl->transform(n, x, y);
/*
    // Special case
    if (n == 1)
    {
        pImpl->mMean = x[0];
        y[0] = x[0];
        return;
    }
    // Remove the mean
    double pMean;
    ippsMean_64f(x, n, &pMean);
    ippsSubC_64f(x, pMean, y, n); 
    pImpl->mMean = pMean;
    // Compute the analytic signal of the data.
    auto yhilb = reinterpret_cast<std::complex<double> *> (pImpl->mHilb64fc);
    pImpl->mHilbert.transform(n, y, &yhilb);
    // Take the absolute value to obtain the envelope.
    // |x_r + i x_h| = x_r^2 + x_h^2 where x_r is the input signal 
    // and x_h is the Hilbert transform of the data.
    ippsMagnitude_64fc(pImpl->mHilb64fc, y, n); 
    // Restore the mean
    ippsAddC_64f_I(pMean, y, n);
*/
}

/// Compute the envelope
/*
template<>
void Envelope<float>::transform(const int n, const float x[], float *yIn[])
{
    pImpl->mMean = 0;
    if (!isInitialized())
    {
        throw std::runtime_error("Class not initialized");
    }
    if (n != pImpl->mEnvelopeLength)
    {
        throw std::invalid_argument("n = %d must equal %d", n, pImpl->mEnvelopeLength);
    }
    float *y = *yIn;
    if (x == nullptr || y == nullptr)
    {
        if (x == nullptr){throw std::invalid_argument("x is NULL");}
        throw std::invalid_argument("y is NULL");
    }
    // Special case
    if (n == 1)
    {
        pImpl->mMean = static_cast<double> (x[0]);
        y[0] = x[0];
        return;
    }
    // Remove the mean
    float pMean;
    ippsMean_32f(x, n, &pMean, ippAlgHintAccurate);
    ippsSubC_32f(x, pMean, y, n); 
    pImpl->mMean = static_cast<double> (pMean);
    // Hilbert transform the signal (yields the analytic signal)
    auto yhilb = reinterpret_cast<std::complex<float> *> (pImpl->mHilb32fc);
    pImpl->mHilbert.transform(n, y, &yhilb);
    // Take the absolute value to obtain the envelope
    ippsMagnitude_32fc(pImpl->mHilb32fc, y, n); 
    // Restore the mean
    ippsAddC_32f_I(pMean, y, n);
}
*/

///--------------------------------------------------------------------------///
///                            Template instantiation                        ///
///--------------------------------------------------------------------------///
template class RTSeis::Transforms::Envelope<double>;
template class RTSeis::Transforms::Envelope<float>;
