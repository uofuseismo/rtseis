#ifndef RTSEIS_UTILS_DESIGN_HPP
#define RTSEIS_UTILS_DESIGN_HPP 1
#include "rtseis/config.h"
#include "rtseis/utils/zpk.hpp"
#include "rtseis/utils/ba.hpp"
#include <complex>
#include <vector>

namespace RTSeis
{
namespace Utils
{
namespace FilterDesign
{

namespace IIR
{
    /* Convert lowpass prototype filter to a lowpass filter. */
    int zpklp2lp(const ZPK zpkIn, const double w0, ZPK &zpkOut);
    /* Convert lowpass prototype filter to a highpass filter. */
    int zpklp2hp(const ZPK zpkIn, const double w0, ZPK &zpkOut);
    /* Convert lowpass prototype filter to a bandpass filter. */
    int zpklp2bp(const ZPK zpkIn, const double w0, const double bw, ZPK &zpkOut);
    /* Convert lowpass prototype filter to a bandstop filter. */
    int zpklp2bs(const ZPK zpkIn, const double w0, const double bw, ZPK &zpkOut);
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
