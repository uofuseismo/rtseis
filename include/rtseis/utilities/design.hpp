#ifndef RTSEIS_UTILS_DESIGN_HPP
#define RTSEIS_UTILS_DESIGN_HPP 1
#include "rtseis/config.h"
#include "rtseis/utilities/zpk.hpp"
#include "rtseis/utilities/ba.hpp"
#include "rtseis/utilities/sos.hpp"
#include <complex>
#include <vector>

namespace RTSeis
{
namespace Utilities
{
namespace FilterDesign
{

namespace IIR
{
    enum Prototype
    {
        /*!
         * @brief Defines the analog prototype for the IIR filter design.
         * @ingroup rtseis_utils_design_iir
         */
        BUTTERWORTH = 0, /*!< Butterworth filter design. */
        BESSEL      = 1, /*!< Bessel filter design. */
        CHEBYSHEV1  = 2, /*!< Chebyshev I filter design. */
        CHEBYSHEV2  = 3  /*!< Chebyshev II filter design. */
    };
    enum Bandtype
    {
        /*!
         * @brief Defines the bandtype for the IIR filter design.
         * @ingroup rtseis_utils_design_iir
         */
        LOWPASS = 0,  /*!< Lowpass filter. */  
        HIGHPASS = 1, /*!< Highpass filter. */
        BANDPASS = 2, /*!< Bandpass filter. */
        BANDSTOP = 3  /*!< Bandstop filter. */
    };
    enum Pairing
    {
        /*!
         * @brief Pairing strategy when converting a ZPK filter to an 
         *        SOS filter.
         * @ingroup rtseis_utils_design_iir
         */ 
        NEAREST = 0,  /*!< This attempts to minimize the peak gain. */
        KEEP_ODD = 1  /*!< This attempts to minimize the peak gain
                           subject to the constraint that odd-order
                           systems should retain one section as first order. */
    };
    /* Generalized analog protoytpe filter design. */
    int iirfilter(const int n, const double *W,
                  const double rp, const double rs,
                  const Bandtype btype,
                  const Prototype ftype,
                  BA &ba,
                  const bool lanalog = false);
    /* Generalized analog protoytpe filter design. */
    int iirfilter(const int n, const double *W,
                  const double rp, const double rs,
                  const Bandtype btype,
                  const Prototype ftype,
                  ZPK &zpk,
                  const bool lanalog = false);
    /* Generalized analog prototype filter design. */
    int iirfilter(const int n, const double *W,
                   const double rp, const double rs,
                   const Bandtype btype,
                   const Prototype ftype,
                   SOS &sos,
                   const bool lanalog,
                   const Pairing pairing = Pairing::NEAREST);
    /* Convert a ZPK structure to second order sections. */
    int zpk2sos(const ZPK zpk, SOS &sos,
                const Pairing pairing = Pairing::NEAREST);
    /* Convert a ZPK structure to a transfer function. */
    int zpk2tf(const ZPK zpk, BA &ba);
    /* Convert lowpass prototype filter to a lowpass filter. */
    int zpklp2lp(const ZPK zpkIn, const double w0, ZPK &zpkOut);
    /* Convert lowpass prototype filter to a highpass filter. */
    int zpklp2hp(const ZPK zpkIn, const double w0, ZPK &zpkOut);
    /* Convert lowpass prototype filter to a bandpass filter. */
    int zpklp2bp(const ZPK zpkIn, const double w0, const double bw, ZPK &zpkOut);
    /* Convert lowpass prototype filter to a bandstop filter. */
    int zpklp2bs(const ZPK zpkIn, const double w0, const double bw, ZPK &zpkOut);
    /* Bilinear transform. */
    int zpkbilinear(const ZPK zpk, const double fs, ZPK &zpkbl);
    namespace AnalogPrototype
    {
        /* Chebyshev I analog prototype. */
        int cheb1ap(const int n, const double rp, ZPK &zpk);
        /* Chebyshev II analog prototype. */
        int cheb2ap(const int n, const double rs, ZPK &zpk);
        /* Butterworth analog prototype. */
        int butter(const int n, ZPK &zpk);
        /* Bessel analog prototype. */
        int bessel(const int n, ZPK &zpk);
    };
}; /* End IIR */

namespace FIR
{
    enum class Window
    {
        /*!
         * @brief Defines the available windows for FIR filter design.
         * @ingroup rtseis_utils_design_fir
         */
        HAMMING  = 0,     /*!< Hamming window. */
        BARTLETT = 1,     /*!< Bartlett (triangle) window. */
        HANN = 2,         /*!< Hann window. */
        BLACKMAN_OPT = 3  /*!< Optimal Blackman window. */ 
    };
    /* Window-based FIR lowpass filter design. */
    int FIR1Lowpass(const int order, const double r,
                    std::vector<double> &taps,
                    const Window window = Window::HAMMING);
    /* Window-based FIR lowpass filter design. */
    int FIR1Highpass(const int order, const double r,
                     std::vector<double> &taps,
                     const Window window = Window::HAMMING);
    /* Window-based FIR bandpass filter design. */
    int FIR1Bandpass(const int order, const double r[2],
                     std::vector<double> &taps,
                     const Window window = Window::HAMMING);
    /* Window-based FIR bandstop filter design. */
    int FIR1Bandstop(const int order, const double r[2],
                     std::vector<double> &taps,
                     const Window window = Window::HAMMING);
}; /* End FIR */

namespace Response
{
    int freqs(const BA ba, const std::vector<double> w,
              std::vector<std::complex<double>> &h);
    int freqz(const BA ba, const std::vector<double> w,
              std::vector<std::complex<double>> &h);
}; /* End Response */

}; /* End FilterDesign */
}; /* End Utils */

}; /* End RTseis */
/*
#ifdef __cplusplus
class AnalogPrototype
{
    public:
        AnalogPrototype(void);
        ~AnalogPrototype(void);
        int cheb1ap(const int n, const double rp);
        int cheb2ap(const int n, const double rs);
        int butter(const int n);
        int bessel(const int n);
        ZPK getTransferFunction(void) const{return zpk_;};
    private:
        ZPK zpk_;
};

class Design : public AnalogPrototype
{
    public:
        Design(void);
        ~Design(void);
        
        int zpklp2lp(const ZPK zpkIn, const double w0, ZPK &zpkOut);
        int zpklp2hp(const ZPK zpkIn, const double w0, ZPK &zpkOut);
        int zpklp2bp(const ZPK zpkIn, const double w0, const double bw, ZPK &zpkOut);
        int zpklp2bs(const ZPK zpkIn, const double w0, const double bw, ZPK &zpkOut);
        int freqs(const BA ba, const std::vector<double> w,
                  std::vector<std::complex<double>> &h);
        int freqz(const BA ba, const std::vector<double> w,
                  std::vector<std::complex<double>> &h);
        
};
#endif
*/

#endif
