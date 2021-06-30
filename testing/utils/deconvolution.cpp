#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "rtseis/filterDesign/response.hpp"
#include "rtseis/filterDesign/iir.hpp"
#include "rtseis/filterRepresentations/ba.hpp"
#include "rtseis/filterRepresentations/zpk.hpp"
#include "rtseis/deconvolution/instruments/woodAnderson.hpp"
#include <gtest/gtest.h>
namespace
{

using namespace RTSeis::Deconvolution::Instruments;
using namespace RTSeis::FilterRepresentations;
using namespace RTSeis::FilterDesign;

TEST(Deconvolution, WoodAnderson)
{
    constexpr std::complex<double> zero(0, 0);
    WoodAnderson wa;
    auto ba = wa.getTransferFunction();
    auto zpk = IIR::tf2zpk(ba);
    auto k = zpk.getGain();
    auto p = zpk.getPoles();
    auto z = zpk.getZeros();

    EXPECT_NEAR(k, 2800, 1.e-10);
    EXPECT_EQ(static_cast<int> (p.size()), 2);
    EXPECT_EQ(static_cast<int> (z.size()), 2);
    // Check zeros
    for (const auto &zi : z)
    {
        EXPECT_NEAR(std::abs(zi - zero), 0, 1.e-14);
    }
    // Check poles (these can come in any order)
    EXPECT_EQ(p[0], std::conj(p[1])); // Complex conjugates
    auto pr = std::real(p[0]);
    auto pi = std::abs(std::imag(p[0]));
    const double t0 = 0.8;
    const double h = 0.8;
    const double om0 = M_PI*2/t0;
    const double rad = std::sqrt(1 - h*h);
    EXPECT_NEAR(pr, -om0*h,   1.e-13); // -6.283185307179586
    EXPECT_NEAR(pi,  om0*rad, 1.e-13); //  4.712388980384689
}

}
