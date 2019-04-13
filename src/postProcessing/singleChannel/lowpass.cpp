#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
//#include "rtseis/postProcessing/singleChannel/lowpass.hpp"

namespace FilterDesign = RTSeis::Utilities::FilterDesign;
using namespace RTSeis::PostProcessing::SingleChannel;

class LowpassFilterDesignParameters : public std::exception
{
     
    private:
        class ParmsImpl;
        std::unique_ptr<ParmsImpl> pParms;
};

class LowpassFilter : public std::exception
{
    public:
        LowpassFilter(void);
        
    private:
        class Impl
        std::unique_ptr<Impl> pImpl;
};

class LowpassFilterDesignParameters::ParmsImpl : public std::exception
{
    public:
        ParsImpl(void)
        {
            return;
        }
        ParsImpl& operator=(const ParmsImpl &impl)
        {
            ba = impls.ba;
        } 
        void clear(void)
        {
            f0 = 0;
            dt = 1;
            prototype = FilterDesign::Prototype::BUTTERWORTH;
            return;
        }
           
        /// Critical (corner) frequency (Hz)
        double f0 = 0;
        /// Sampling period (seconds) 
        double dt = 1;
        /// Set the nqyuist frequency
        double fnyq = 1/(2*dt);
        /// Analog prototype
        const FilterDesign::Prototype prototype
             = FilterDesign::Prototype::BUTTERWORTH;
        /// Window-based FIR design

        /// FIR filter design
        bool isFIR = false;
        /// 
        bool isInitialized = false;
};

class LowpassFilter::Impl : public std::exception
{
    public:
};

/*
class Lowpass::LowpassImpl
{
    public:
        
};
*/

void LowpassIIRFilterParameters::LowpassFilterParameters(void)
{
    return;
}

void LowpassIIRFilterParameters::LowpassFilterParameters( 
    const LowpassFilterParameters &lowpass)
{
    *this = lowpass;
    return;
}

bool LowpassIIRFilterDesign::isFIR(void) const
{
    return  
}

//============================================================================//

LowpassFilter::LowpassFilter(void) :
    pImpl(new ParmsImpl())
{
    return;
}

LowpassFilter::LowpassFilter(const LowpassIIRFilterParameters &parms) :
    pImpl(new ParmsImpl())
{
    initialize(parms);
}

void LowpassFilter::initialize(const LowpassIIRFilterParameters &parms)
{
    clear();
    if (!parms.isValid())
    {
        throw std::invalid_argument("Filter parameters are invalid");
    }
    pImpl->parms = parms;
    return;
}
/*

void LowpassIIRFilterDesign::setCriticalFrequency(const double f0, const double dt)
{

}
 
void LowpassFilterDesign::setPrototype(
    const RTSeis::Utilities::FilterDesign::Prototype prototype, = RTSeis::Utilities::FilterDesign::Prototype::BUTTERWORTH)
    const RTSeis::Utilities::FilterRepresentations::Representation representation RTSeis::Utilities::FilterRepresentations::Representation::SOS,
    const FilterAppliciation application=FilterApplication::ZERO_PHASE)
{

}

Lowpass::Lowpass( )
{
    
}
*/
