#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "rtseis/utilities/design/enums.hpp"
#include "rtseis/utilities/design/filterDesigner.hpp"
#include "rtseis/utilities/design/fir.hpp"
#include "rtseis/utilities/design/iir.hpp"
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterRepresentations/sos.hpp"
#include "rtseis/utilities/filterRepresentations/zpk.hpp"

using namespace RTSeis::Utilities;
using namespace RTSeis::Utilities::FilterDesign;


struct FIRDesignParameters
{
    FIRDesignParameters(void){return;}
    ~FIRDesignParameters(void){return;}
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
        return;
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
        return;
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
    void clear(void)
    {
        r1 = 0;
        r2 = 0;
        order = 0;
        window = FIRWindow::HAMMING;
        btype = Bandtype::LOWPASS;
        return;
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
    IIRDesignParameters(void){return;}
    ~IIRDesignParameters(void){return;}
    void clear(void)
    {
        r1 = 0;
        r2 = 0;
        order = 0;
        prototype = IIRPrototype::BUTTERWORTH;
        btype = Bandtype::LOWPASS;
        pairing = SOSPairing::NEAREST;
        lsos = false;
        ldigital = IIRFilterDomain::DIGITAL;
    }
    IIRDesignParameters& operator=(const IIRDesignParameters &parms)
    {   
        if (&parms == this){return *this;}
        r1 = parms.r1;
        r2 = parms.r2;
        order = parms.order;
        prototype  = parms.prototype;
        btype  = parms.btype;
        pairing = parms.pairing;
        lsos = parms.lsos;
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
            if(r2 != parms.r2){return false;}
        }
        if (ldigital != parms.ldigital){return false;} 
        if (( lsos && !parms.lsos) ||
            (!lsos &&  parms.lsos)){return false;}
        if (lsos)
        {
            if (pairing != parms.pairing){return false;}
        }
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
    /// Filter order
    int order = 0;
    /// The window type
    IIRPrototype prototype = IIRPrototype::BUTTERWORTH;
    /// The filter band
    Bandtype btype = Bandtype::LOWPASS;
    /// Pole pairing
    SOSPairing pairing = SOSPairing::NEAREST;
    /// SOS?
    bool lsos = false;
    /// Digital?
    IIRFilterDomain ldigital = IIRFilterDomain::DIGITAL;
};
//----------------------------------------------------------------------------//
class FilterDesigner::FilterDesignerImpl
{
    public:
        FilterDesignerImpl(void)
        {
            zpkDesigns.reserve(64);
            zpkCache.reserve(64);
            baDesigns.reserve(64);
            baCache.reserve(64);
            sosDesigns.reserve(64);
            sosCache.reserve(64);
            firDesigns.reserve(64);
            firCache.reserve(64);
            return;
        }
        ~FilterDesignerImpl(void)
        {
            clear();
            return;
        }
        void clear(void)
        {
            zpkDesigns.clear();
            zpkCache.clear();
            baDesigns.clear();
            baCache.clear();
            sosDesigns.clear();
            sosCache.clear();
            firDesigns.clear();
            firCache.clear();
            return;
        }

        std::vector<IIRDesignParameters> zpkDesigns;
        std::vector<FilterRepresentations::ZPK> zpkCache;
        std::vector<IIRDesignParameters> baDesigns;
        std::vector<FilterRepresentations::BA> baCache;
        std::vector<IIRDesignParameters> sosDesigns;
        std::vector<FilterRepresentations::SOS> sosCache; 
        std::vector<FIRDesignParameters> firDesigns;
        std::vector<FilterRepresentations::FIR> firCache;
};

//=============================================================================//

FilterDesigner::FilterDesigner(void) :
    pImpl(new FilterDesignerImpl())
{
    return;
}

FilterDesigner::FilterDesigner(const FilterDesigner &design)
{
    *this = design;
    return;
}

FilterDesigner::~FilterDesigner(void)
{
    if (pImpl){pImpl->clear();}
}

FilterDesigner& FilterDesigner::operator=(const FilterDesigner &design)
{
    if (&design == this){return *this;}
    if (pImpl){pImpl->clear();}
    pImpl = std::unique_ptr<FilterDesignerImpl> (new FilterDesignerImpl());
    pImpl->zpkDesigns = design.pImpl->zpkDesigns;
    pImpl->zpkCache   = design.pImpl->zpkCache;
    pImpl->baDesigns  = design.pImpl->baDesigns;
    pImpl->baCache    = design.pImpl->baCache;
    pImpl->sosDesigns = design.pImpl->sosDesigns;
    pImpl->sosCache   = design.pImpl->sosCache;
    pImpl->firDesigns = design.pImpl->firDesigns;
    pImpl->firCache   = design.pImpl->firCache;
    return *this;
}

void FilterDesigner::clear(void)
{
    if (pImpl){pImpl->clear();}
}

//============================================================================//

//void FilterDesigner::designLowpassIIRFilter(
//    const int order,
//============================================================================//

void FilterDesigner::designLowpassFIRFilter(
    const int order,
    const double r,
    const FIRWindow window,
    FilterRepresentations::FIR &fir) const
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
        FIR::FIR1Lowpass(order, r, fir, window); // Throws error
        pImpl->firDesigns.push_back(parms); 
        pImpl->firCache.push_back(fir);
    }
    return;
}

void FilterDesigner::designHighpassFIRFilter(
    const int order,
    const double r,
    const FIRWindow window,
    FilterRepresentations::FIR &fir) const
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
        FIR::FIR1Highpass(order, r, fir, window); // Throws error
        pImpl->firDesigns.push_back(parms); 
        pImpl->firCache.push_back(fir);
    }
    return;
}

void FilterDesigner::designBandpassFIRFilter(
    const int order,
    const std::pair<double,double> r,
    const FIRWindow window,
    FilterRepresentations::FIR &fir) const
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
        FIR::FIR1Bandpass(order, r, fir, window); // Throws error
        pImpl->firDesigns.push_back(parms); 
        pImpl->firCache.push_back(fir);
    }
    return;
}

void FilterDesigner::designBandstopFIRFilter(
    const int order,
    const std::pair<double,double> r,
    const FIRWindow window,
    FilterRepresentations::FIR &fir) const
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
        FIR::FIR1Bandstop(order, r, fir, window); // Throws error
        pImpl->firDesigns.push_back(parms);
        pImpl->firCache.push_back(fir);
    }
    return;
}
