#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "rtseis/filterDesign/enums.hpp"
#include "rtseis/filterDesign/filterDesigner.hpp"
#include "rtseis/filterDesign/fir.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/fir.hpp"
#include "rtseis/filterRepresentations/sos.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"

using namespace RTSeis::FilterDesign;

namespace
{
std::pair<double,double> 
iirPrototypeToRipple(const IIRPrototype ftype, const double r)
{
    std::pair<double,double> ripple(0,0);
    if (ftype == IIRPrototype::CHEBYSHEV1)
    {   
        ripple = std::make_pair(r, 0); 
    }   
    else if (ftype == IIRPrototype::CHEBYSHEV2)
    {   
        ripple = std::make_pair(0, r); 
    }   
    return ripple;
}
}

struct FIRDesignParameters
{
    FIRDesignParameters(const int orderIn,
                        const double r,
                        const FIRWindow windowIn,
                        const Bandtype btypeIn) :
        r1(r),
        r2(0),
        order(orderIn),
        window(windowIn),
        btype(btypeIn)
    {
    } 
    FIRDesignParameters(const int orderIn,
                        const std::pair<double,double> r,
                        const FIRWindow windowIn,
                        const Bandtype btypeIn) :
        r1(r.first),
        r2(r.second),
        order(orderIn),
        window(windowIn),
        btype(btypeIn)
    {
    } 
    FIRDesignParameters& operator=(const FIRDesignParameters &parms)
    {
        if (&parms == this){return *this;}
        r1 = parms.r1;
        r2 = parms.r2;
        order = parms.order;
        window = parms.window;
        btype  = parms.btype;
        return *this;
    }
    bool operator==(const FIRDesignParameters &parms) const
    {
        if (window != parms.window){return false;}
        if (btype  != parms.btype){return false;}
        if (r1 != parms.r1){return false;}
        if (btype == Bandtype::BANDPASS ||
            btype == Bandtype::BANDSTOP)
        {
            if(r2 != parms.r2){return false;}
        }
        return true; 
    }
    bool operator!=(const FIRDesignParameters &parms) const
    {
        return !(*this == parms);
    }
    void clear() noexcept
    {
        r1 = 0;
        r2 = 0;
        order = 0;
        window = FIRWindow::HAMMING;
        btype = Bandtype::LOWPASS;
    }

    /// First critical frequency
    double r1 = 0;
    /// Second critical frequency
    double r2 = 0;
    /// Filter order
    int order = 0;
    /// The window type
    FIRWindow window = FIRWindow::HAMMING;
    /// The filter band
    Bandtype btype = Bandtype::LOWPASS;
};
//----------------------------------------------------------------------------//
struct IIRDesignParameters
{
    /// IIR direct form constructor for lowpass/highpass
    IIRDesignParameters(const int orderIn,
                        const double r,
                        const double rippleIn,
                        const IIRPrototype prototypeIn,
                        const Bandtype btypeIn,
                        const IIRFilterDomain ldigitalIn) :
        r1(r),
        r2(0),
        ripple(rippleIn),
        order(orderIn),
        prototype(prototypeIn),
        btype(btypeIn),
        ldigital(ldigitalIn)
    {
    }
    /// IIR direct form constructor for bandpass/bandstop
    IIRDesignParameters(const int orderIn,
                        const std::pair<double,double> r,
                        const double rippleIn,
                        const IIRPrototype prototypeIn,
                        const Bandtype btypeIn,
                        const IIRFilterDomain ldigitalIn) :
        r1(r.first),
        r2(r.second),
        ripple(rippleIn),
        order(orderIn),
        prototype(prototypeIn),
        btype(btypeIn),
        ldigital(ldigitalIn)
    {
    }
    void clear() noexcept
    {
        r1 = 0;
        r2 = 0;
        ripple = 0;
        order = 0;
        prototype = IIRPrototype::BUTTERWORTH;
        btype = Bandtype::LOWPASS;
        ldigital = IIRFilterDomain::DIGITAL;
    }
    IIRDesignParameters& operator=(const IIRDesignParameters &parms)
    {   
        if (&parms == this){return *this;}
        r1 = parms.r1;
        r2 = parms.r2;
        ripple = parms.ripple;
        order = parms.order;
        prototype  = parms.prototype;
        btype  = parms.btype;
        ldigital = parms.ldigital;
        return *this;
    }
    bool operator==(const IIRDesignParameters &parms) const
    {   
        if (prototype != parms.prototype){return false;}
        if (btype  != parms.btype){return false;}
        if (r1 != parms.r1){return false;}
        if (btype == Bandtype::BANDPASS ||
            btype == Bandtype::BANDSTOP)
        {
            if (r2 != parms.r2){return false;}
        }
        if (prototype == IIRPrototype::CHEBYSHEV1 ||
            prototype == IIRPrototype::CHEBYSHEV2)
        {
            if (ripple != parms.ripple){return false;}
        }
        if (ldigital != parms.ldigital){return false;} 
        return true; 
    }
    bool operator!=(const IIRDesignParameters &parms) const
    {   
        return !(*this == parms);
    }
    /// First critical frequency
    double r1 = 0;
    /// Second critical frequency
    double r2 = 0;
    /// The ripple for Cheby1 and Cheby2 filters
    double ripple = 0;
    /// Filter order
    int order = 0;
    /// The window type
    IIRPrototype prototype = IIRPrototype::BUTTERWORTH;
    /// The filter band
    Bandtype btype = Bandtype::LOWPASS;
    /// Digital?
    IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL;
};
//----------------------------------------------------------------------------//
struct SOSDesignParameters
{
    /// SOS constructor for lowpass/highpass
    SOSDesignParameters(const int orderIn,
                        const double r,
                        const double rippleIn,
                        const IIRPrototype prototypeIn,
                        const Bandtype btypeIn,
                        const SOSPairing pairingIn,
                        const IIRFilterDomain ldigitalIn) :
        r1(r),
        r2(0),
        ripple(rippleIn),
        order(orderIn),
        prototype(prototypeIn),
        btype(btypeIn),
        pairing(pairingIn),
        ldigital(ldigitalIn)
    {
    }
    /// IIR sos constructor for bandpass/bandstop
    SOSDesignParameters(const int orderIn,
                        const std::pair<double,double> r,
                        const double rippleIn,
                        const IIRPrototype prototypeIn,
                        const Bandtype btypeIn,
                        const SOSPairing pairingIn,
                        const IIRFilterDomain ldigitalIn) :
        r1(r.first),
        r2(r.second),
        ripple(rippleIn),
        order(orderIn),
        prototype(prototypeIn),
        btype(btypeIn),
        pairing(pairingIn),
        ldigital(ldigitalIn)
    {   
    }
    void clear() noexcept
    {
        r1 = 0;
        r2 = 0;
        ripple = 0;
        order = 0;
        prototype = IIRPrototype::BUTTERWORTH;
        btype = Bandtype::LOWPASS;
        pairing = SOSPairing::NEAREST;
        ldigital = IIRFilterDomain::DIGITAL;
    }
    SOSDesignParameters& operator=(const SOSDesignParameters &parms)
    {   
        if (&parms == this){return *this;}
        r1 = parms.r1;
        r2 = parms.r2;
        ripple = parms.ripple;
        order = parms.order;
        prototype  = parms.prototype;
        btype  = parms.btype;
        pairing = parms.pairing;
        ldigital = parms.ldigital;
        return *this;
    }
    bool operator==(const SOSDesignParameters &parms) const
    {   
        if (prototype != parms.prototype){return false;}
        if (btype  != parms.btype){return false;}
        if (r1 != parms.r1){return false;}
        if (btype == Bandtype::BANDPASS ||
            btype == Bandtype::BANDSTOP)
        {
            if (r2 != parms.r2){return false;}
        }
        if (prototype == IIRPrototype::CHEBYSHEV1 ||
            prototype == IIRPrototype::CHEBYSHEV2)
        {
            if (ripple != parms.ripple){return false;}
        }
        if (ldigital != parms.ldigital){return false;} 
        if (pairing != parms.pairing){return false;}
        return true; 
    }
    bool operator!=(const SOSDesignParameters &parms) const
    {   
        return !(*this == parms);
    }
    /// First critical frequency
    double r1 = 0;
    /// Second critical frequency
    double r2 = 0;
    /// The ripple for Cheby1 and Cheby2 filters
    double ripple = 0;
    /// Filter order
    int order = 0;
    /// The window type
    IIRPrototype prototype = IIRPrototype::BUTTERWORTH;
    /// The filter band
    Bandtype btype = Bandtype::LOWPASS;
    /// Pole pairing
    SOSPairing pairing = SOSPairing::NEAREST;
    /// Digital?
    IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL;
};
//----------------------------------------------------------------------------//
class FilterDesigner::FilterDesignerImpl
{
public:
    FilterDesignerImpl()
    {
        zpkDesigns.reserve(64);
        zpkCache.reserve(64);
        baDesigns.reserve(64);
        baCache.reserve(64);
        sosDesigns.reserve(64);
        sosCache.reserve(64);
        firDesigns.reserve(64);
        firCache.reserve(64);
    }
    void clear() noexcept
    {
        zpkDesigns.clear();
        zpkCache.clear();
        baDesigns.clear();
        baCache.clear();
        sosDesigns.clear();
        sosCache.clear();
        firDesigns.clear();
        firCache.clear();
    }

    std::vector<IIRDesignParameters> zpkDesigns;
    std::vector<RTSeis::FilterRepresentations::ZPK> zpkCache;
    std::vector<IIRDesignParameters> baDesigns;
    std::vector<RTSeis::FilterRepresentations::BA> baCache;
    std::vector<SOSDesignParameters> sosDesigns;
    std::vector<RTSeis::FilterRepresentations::SOS> sosCache; 
    std::vector<FIRDesignParameters> firDesigns;
    std::vector<RTSeis::FilterRepresentations::FIR> firCache;
};

//=============================================================================//

/// C'tor
FilterDesigner::FilterDesigner() :
    pImpl(std::make_unique<FilterDesignerImpl>())
{
}

/// Copy c'tor
FilterDesigner::FilterDesigner(const FilterDesigner &design)
{
    *this = design;
}

FilterDesigner::FilterDesigner(FilterDesigner &&design) noexcept
{
    *this = std::move(design);
}

/// Destructor
FilterDesigner::~FilterDesigner() = default;

/// Copy assignment
FilterDesigner& FilterDesigner::operator=(const FilterDesigner &design)
{
    if (&design == this){return *this;}
    pImpl = std::make_unique<FilterDesignerImpl> (*design.pImpl);
    return *this;
}

/// Move assignment
FilterDesigner& FilterDesigner::operator=(FilterDesigner &&design) noexcept
{
    if (&design == this){return *this;}
    pImpl = std::move(design.pImpl);
    return *this;
}

/// Reset class
void FilterDesigner::clear() noexcept
{
    if (pImpl){pImpl->clear();}
}

//============================================================================//

void FilterDesigner::designLowpassIIRFilter(
    const int n, const double r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::ZPK &zpk,
    const IIRFilterDomain ldigital)
{
    zpk.clear();
    IIRDesignParameters parms(n, r, ripple, ftype, Bandtype::LOWPASS, ldigital);
    // Look for design
    auto it = std::find(pImpl->zpkDesigns.begin(),
                        pImpl->zpkDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->zpkDesigns.end())
    {
        size_t indx = std::distance(pImpl->zpkDesigns.begin(), it);
        zpk = pImpl->zpkCache[indx];
    }
    else
    {
        double W[1] = {r};
        std::pair<double, double> rp = iirPrototypeToRipple(ftype, ripple);
        zpk = IIR::designZPKIIRFilter(n, W, rp.first, rp.second,
                                      Bandtype::LOWPASS, ftype, ldigital);
        pImpl->zpkDesigns.push_back(parms);
        pImpl->zpkCache.push_back(zpk);
    }   
}

void FilterDesigner::designHighpassIIRFilter(
    const int n, const double r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::ZPK &zpk,
    const IIRFilterDomain ldigital)
{
    zpk.clear();
    IIRDesignParameters parms(n, r, ripple, ftype, Bandtype::HIGHPASS,ldigital);
    // Look for design
    auto it = std::find(pImpl->zpkDesigns.begin(),
                        pImpl->zpkDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->zpkDesigns.end())
    {
        size_t indx = std::distance(pImpl->zpkDesigns.begin(), it);
        zpk = pImpl->zpkCache[indx];
    }
    else
    {
        double W[1] = {r};
        std::pair<double, double> rp = iirPrototypeToRipple(ftype, ripple);
        zpk = IIR::designZPKIIRFilter(n, W, rp.first, rp.second,
                                      Bandtype::HIGHPASS, ftype, ldigital);
        pImpl->zpkDesigns.push_back(parms);
        pImpl->zpkCache.push_back(zpk);
    }
}

void FilterDesigner::designBandpassIIRFilter(
    const int n, const std::pair<double,double> &r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::ZPK &zpk,
    const IIRFilterDomain ldigital)
{
    zpk.clear();
    IIRDesignParameters parms(n, r, ripple, ftype, Bandtype::BANDPASS,ldigital);
    // Look for design
    auto it = std::find(pImpl->zpkDesigns.begin(),
                        pImpl->zpkDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->zpkDesigns.end())
    {
        size_t indx = std::distance(pImpl->zpkDesigns.begin(), it);
        zpk = pImpl->zpkCache[indx];
    }
    else
    {
        double W[2] = {r.first, r.second};
        std::pair<double, double> rp = iirPrototypeToRipple(ftype, ripple);
        zpk = IIR::designZPKIIRFilter(n, W, rp.first, rp.second,
                                      Bandtype::BANDPASS, ftype, ldigital);
        pImpl->zpkDesigns.push_back(parms);
        pImpl->zpkCache.push_back(zpk);
    }
}

void FilterDesigner::designBandstopIIRFilter(
    const int n, const std::pair<double,double> &r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::ZPK &zpk,
    const IIRFilterDomain ldigital)
{
    zpk.clear();
    IIRDesignParameters parms(n, r, ripple, ftype, Bandtype::BANDSTOP,ldigital);
    // Look for design
    auto it = std::find(pImpl->zpkDesigns.begin(),
                        pImpl->zpkDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->zpkDesigns.end())
    {
        size_t indx = std::distance(pImpl->zpkDesigns.begin(), it);
        zpk = pImpl->zpkCache[indx];
    }
    else
    {
        double W[2] = {r.first, r.second};
        std::pair<double, double> rp = iirPrototypeToRipple(ftype, ripple);
        zpk = IIR::designZPKIIRFilter(n, W, rp.first, rp.second,
                                      Bandtype::BANDSTOP, ftype, ldigital);
        pImpl->zpkDesigns.push_back(parms);
        pImpl->zpkCache.push_back(zpk);
    }
}

//============================================================================//

void FilterDesigner::designLowpassIIRFilter(
    const int n, const double r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::BA &ba,
    const IIRFilterDomain ldigital)
{
    ba.clear();
    IIRDesignParameters parms(n, r, ripple, ftype, Bandtype::LOWPASS, ldigital);
    // Look for design
    auto it = std::find(pImpl->baDesigns.begin(),
                        pImpl->baDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->baDesigns.end())
    {
        size_t indx = std::distance(pImpl->baDesigns.begin(), it);
        ba = pImpl->baCache[indx];
    }
    else
    {
        RTSeis::FilterRepresentations::ZPK zpk;
        designLowpassIIRFilter(n, r, ftype, ripple, zpk, ldigital);
        ba = IIR::zpk2tf(zpk);
        pImpl->baDesigns.push_back(parms);
        pImpl->baCache.push_back(ba);
    }
}

void FilterDesigner::designHighpassIIRFilter(
    const int n, const double r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::BA &ba,
    const IIRFilterDomain ldigital)
{
    ba.clear();
    IIRDesignParameters parms(n, r, ripple, ftype, Bandtype::HIGHPASS,ldigital);
    // Look for design
    auto it = std::find(pImpl->baDesigns.begin(),
                        pImpl->baDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->baDesigns.end())
    {
        size_t indx = std::distance(pImpl->baDesigns.begin(), it);
        ba = pImpl->baCache[indx];
    }
    else
    {
        RTSeis::FilterRepresentations::ZPK zpk;
        designHighpassIIRFilter(n, r, ftype, ripple, zpk, ldigital);
        ba = IIR::zpk2tf(zpk);
        pImpl->baDesigns.push_back(parms);
        pImpl->baCache.push_back(ba);
    }
}

void FilterDesigner::designBandpassIIRFilter(
    const int n, const std::pair<double,double> &r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::BA &ba,
    const IIRFilterDomain ldigital)
{
    ba.clear();
    IIRDesignParameters parms(n, r, ripple, ftype, Bandtype::BANDPASS,ldigital);
    // Look for design
    auto it = std::find(pImpl->baDesigns.begin(),
                        pImpl->baDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->baDesigns.end())
    {
        size_t indx = std::distance(pImpl->baDesigns.begin(), it);
        ba = pImpl->baCache[indx];
    }
    else
    {
        RTSeis::FilterRepresentations::ZPK zpk;
        designBandpassIIRFilter(n, r, ftype, ripple, zpk, ldigital);
        ba = IIR::zpk2tf(zpk);
        pImpl->baDesigns.push_back(parms);
        pImpl->baCache.push_back(ba);
    }
}

void FilterDesigner::designBandstopIIRFilter(
    const int n, const std::pair<double,double> &r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::BA &ba,
    const IIRFilterDomain ldigital)
{
    ba.clear();
    IIRDesignParameters parms(n, r, ripple, ftype, Bandtype::BANDSTOP,ldigital);
    // Look for design
    auto it = std::find(pImpl->baDesigns.begin(),
                        pImpl->baDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->baDesigns.end())
    {
        size_t indx = std::distance(pImpl->baDesigns.begin(), it);
        ba = pImpl->baCache[indx];
    }
    else
    {
        RTSeis::FilterRepresentations::ZPK zpk;
        designBandstopIIRFilter(n, r, ftype, ripple, zpk, ldigital);
        ba = IIR::zpk2tf(zpk);
        pImpl->baDesigns.push_back(parms);
        pImpl->baCache.push_back(ba);
    }
}

//============================================================================//

void FilterDesigner::designLowpassIIRFilter(
    const int n, const double r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::SOS &sos,
    const SOSPairing pairing,   
    const IIRFilterDomain ldigital)
{
    sos.clear();
    SOSDesignParameters parms(n, r, ripple, ftype, Bandtype::LOWPASS,
                              pairing, ldigital);
    // Look for design
    auto it = std::find(pImpl->sosDesigns.begin(),
                        pImpl->sosDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->sosDesigns.end())
    {
        size_t indx = std::distance(pImpl->sosDesigns.begin(), it);
        sos = pImpl->sosCache[indx];
    }
    else
    {
        RTSeis::FilterRepresentations::ZPK zpk;
        designLowpassIIRFilter(n, r, ftype, ripple, zpk, ldigital);
        sos = IIR::zpk2sos(zpk);
        pImpl->sosDesigns.push_back(parms);
        pImpl->sosCache.push_back(sos);
    }
    return;
}

void FilterDesigner::designHighpassIIRFilter(
    const int n, const double r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::SOS &sos,
    const SOSPairing pairing,
    const IIRFilterDomain ldigital)
{
    sos.clear();
    SOSDesignParameters parms(n, r, ripple, ftype, Bandtype::HIGHPASS,
                              pairing, ldigital);
    // Look for design
    auto it = std::find(pImpl->sosDesigns.begin(),
                        pImpl->sosDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->sosDesigns.end())
    {
        size_t indx = std::distance(pImpl->sosDesigns.begin(), it);
        sos = pImpl->sosCache[indx];
    }
    else
    {
        RTSeis::FilterRepresentations::ZPK zpk;
        designHighpassIIRFilter(n, r, ftype, ripple, zpk, ldigital);
        sos = IIR::zpk2sos(zpk);
        pImpl->sosDesigns.push_back(parms);
        pImpl->sosCache.push_back(sos);
    }
}

void FilterDesigner::designBandpassIIRFilter(
    const int n, const std::pair<double,double> &r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::SOS &sos,
    const SOSPairing pairing,
    const IIRFilterDomain ldigital)
{
    sos.clear();
    SOSDesignParameters parms(n, r, ripple, ftype, Bandtype::BANDPASS,
                              pairing, ldigital);
    // Look for design
    auto it = std::find(pImpl->sosDesigns.begin(),
                        pImpl->sosDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->sosDesigns.end())
    {
        size_t indx = std::distance(pImpl->sosDesigns.begin(), it);
        sos = pImpl->sosCache[indx];
    }
    else
    {
        RTSeis::FilterRepresentations::ZPK zpk;
        designBandpassIIRFilter(n, r, ftype, ripple, zpk, ldigital);
        sos = IIR::zpk2sos(zpk);
        pImpl->sosDesigns.push_back(parms);
        pImpl->sosCache.push_back(sos);
    }
}

void FilterDesigner::designBandstopIIRFilter(
    const int n, const std::pair<double,double> &r,
    const IIRPrototype ftype,
    const double ripple,
    RTSeis::FilterRepresentations::SOS &sos,
    const SOSPairing pairing,
    const IIRFilterDomain ldigital)
{
    sos.clear();
    SOSDesignParameters parms(n, r, ripple, ftype, Bandtype::BANDSTOP,
                              pairing, ldigital);
    // Look for design
    auto it = std::find(pImpl->sosDesigns.begin(),
                        pImpl->sosDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->sosDesigns.end())
    {
        size_t indx = std::distance(pImpl->sosDesigns.begin(), it);
        sos = pImpl->sosCache[indx];
    }
    else
    {
        RTSeis::FilterRepresentations::ZPK zpk;
        designBandstopIIRFilter(n, r, ftype, ripple, zpk, ldigital);
        sos = IIR::zpk2sos(zpk);
        pImpl->sosDesigns.push_back(parms);
        pImpl->sosCache.push_back(sos);
    }
}

//============================================================================//

void FilterDesigner::designLowpassFIRFilter(
    const int order,
    const double r,
    const FIRWindow window,
    RTSeis::FilterRepresentations::FIR &fir) const
{
    fir.clear();
    FIRDesignParameters parms(order, r, window, Bandtype::HIGHPASS);
    // Look for this design
    auto it = std::find(pImpl->firDesigns.begin(),
                        pImpl->firDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->firDesigns.end())
    {
        size_t indx = std::distance(pImpl->firDesigns.begin(), it);
        fir = pImpl->firCache[indx];
    }
    else
    {
        fir = FIR::FIR1Lowpass(order, r, window); // Throws error
        pImpl->firDesigns.push_back(parms); 
        pImpl->firCache.push_back(fir);
    }
}

void FilterDesigner::designHighpassFIRFilter(
    const int order,
    const double r,
    const FIRWindow window,
    RTSeis::FilterRepresentations::FIR &fir) const
{
    fir.clear();
    FIRDesignParameters parms(order, r, window, Bandtype::HIGHPASS);
    // Look for this design
    auto it = std::find(pImpl->firDesigns.begin(),
                        pImpl->firDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->firDesigns.end())
    {
        size_t indx = std::distance(pImpl->firDesigns.begin(), it);
        fir = pImpl->firCache[indx];
    }
    else
    {
        fir = FIR::FIR1Highpass(order, r, window); // Throws error
        pImpl->firDesigns.push_back(parms); 
        pImpl->firCache.push_back(fir);
    }
}

void FilterDesigner::designBandpassFIRFilter(
    const int order,
    const std::pair<double,double> &r,
    const FIRWindow window,
    RTSeis::FilterRepresentations::FIR &fir) const
{
    fir.clear();
    FIRDesignParameters parms(order, r, window, Bandtype::BANDPASS);
    // Look for this design
    auto it = std::find(pImpl->firDesigns.begin(),
                        pImpl->firDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->firDesigns.end())
    {
        size_t indx = std::distance(pImpl->firDesigns.begin(), it);
        fir = pImpl->firCache[indx];
    }
    else
    {
        fir = FIR::FIR1Bandpass(order, r, window); // Throws error
        pImpl->firDesigns.push_back(parms); 
        pImpl->firCache.push_back(fir);
    }
}

void FilterDesigner::designBandstopFIRFilter(
    const int order,
    const std::pair<double,double> &r,
    const FIRWindow window,
    RTSeis::FilterRepresentations::FIR &fir) const
{
    fir.clear();
    FIRDesignParameters parms(order, r, window, Bandtype::BANDSTOP);
    // Look for this design
    auto it = std::find(pImpl->firDesigns.begin(),
                        pImpl->firDesigns.end(), parms);
    // Found! Copy the design
    if (it != pImpl->firDesigns.end())
    {
        size_t indx = std::distance(pImpl->firDesigns.begin(), it);
        fir = pImpl->firCache[indx];
    }
    else
    {
        fir = FIR::FIR1Bandstop(order, r, window); // Throws error
        pImpl->firDesigns.push_back(parms);
        pImpl->firCache.push_back(fir);
    }
}

