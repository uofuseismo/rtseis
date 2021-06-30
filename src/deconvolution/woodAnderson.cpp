#include <vector>
#include <cmath>
#include "rtseis/deconvolution/instruments/woodAnderson.hpp"
#include "rtseis/filterRepresentations/ba.hpp"

#define T0 0.8
#define H  0.8
#define GAIN 2800
// rad = sqrt(1 - h**2) = sqrt(1 - 0.8^2) = sqrt(1 - 0.64) = sqrt(0.36)
#define RAD 0.6 
// om0 = twopoi/t0
// note 2/0.8 = 2.5
#define OM0 2.5*M_PI
#define TWOR 2*OM0*H
#define MAG (OM0*OM0)*(H*H + RAD*RAD)

using namespace RTSeis::Deconvolution::Instruments;
/*
 
class WoodAnderson::WoodAndersonImpl
{
public:
    /// Transfer function representation of Wood-Anderson instrument.
    /// This uses h = 0.8, t0 = 0.8, and k = 2800.
    /// The zeros are: {0, 0}.
    /// The poles are: { -2*pi*h/t0 - 2*pi/t0*sqrt(1 - h^2)i,
    ///                  -2*pi*h/t0 + 2*pi/t0*sqrt(1 - h^2)i }
    /// Expanding gain*(z - z_0)*(z - z_0) and (p- p_0)*(p - \bar{p}_0) we get
    std::vector<double> b{GAIN*1, 0, 0};
    std::vector<double> a{1, TWOR, MAG};
    RTSeis::FilterRepresentations::BA mBA{b, a};
    /// This filter design is taken directly from the mechanical instrument 
    /// so it must be analog.
    const bool mAnalog = true;
};

/// C'tor
WoodAnderson::WoodAnderson() :
    pImpl(std::make_unique<WoodAndersonImpl> ())
{
}

/// Destructor
WoodAnderson::~WoodAnderson() = default;
*/

/// 
RTSeis::FilterRepresentations::BA 
WoodAnderson::getTransferFunction() noexcept
{
    const std::vector<double> b{GAIN*1, 0, 0};
    const std::vector<double> a{1, TWOR, MAG};
    RTSeis::FilterRepresentations::BA ba{b, a};
    return ba;
}

