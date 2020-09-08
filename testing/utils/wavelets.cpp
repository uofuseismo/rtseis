#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <exception>
#include <vector>
#include <ipps.h>
//#include "rtseis/utilities/transforms/wavelets/derivativeOfGaussian.hpp"
#include "rtseis/utilities/transforms/wavelets/morlet.hpp"
//#include "rtseis/utilities/transforms/wavelets/ricker.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis::Utilities::Transforms;

/*
TEST(UtilitiesTransformsWavelets, dog)
{
    int order = 3;
    double dt = 0.25;
    int nk = 20;
    double scale = 13.454342644059432;
    int npad = 1024;
    double df = (2*M_PI)/(dt*npad);
    std::vector<std::complex<double>> refWavelet1({
        std::complex<double> (0, 0), 
        std::complex<double> (-0.0000000000000000,0.34394694837912410),
        std::complex<double> (-0.0000000000000000,2.3363889616734714),
        std::complex<double> (-0.0000000000000000,6.0037843787718908),
        std::complex<double> (-0.0000000000000000,9.7160541230762583),
        std::complex<double> (-0.0000000000000000,11.617456603642266),
        std::complex<double> (-0.0000000000000000,11.020171594204173),
        std::complex<double> (-0.0000000000000000,8.6139976943725429),
        std::complex<double> (-0.0000000000000000,5.6754340694418390),
        std::complex<double> (-0.0000000000000000,3.1982882790899012),
        std::complex<double> (-0.0000000000000000,1.5570175216301099),
        std::complex<double> (-0.0000000000000000,0.65950401289063132),
        std::complex<double> (-0.0000000000000000,0.24432700854476191),
        std::complex<double> (-0.0000000000000000,7.94856288399509447E-002),
        std::complex<double> (-0.0000000000000000,2.27780188489948542E-002),
        std::complex<double> (-0.0000000000000000,5.76396542448160461E-003),
        std::complex<double> (-0.0000000000000000,1.29052300914506471E-003),
        std::complex<double> (-0.0000000000000000,2.56066784552561749E-004),
        std::complex<double> (-0.0000000000000000,4.50886408777352387E-005),
        std::complex<double> (-0.0000000000000000,7.05335987033177645E-006)
    } );
    // Initialize the wavelet
    Wavelets::DerivativeOfGaussian dog;
    EXPECT_NO_THROW(dog.setOrder(order));
    EXPECT_NO_THROW(dog.setSamplingPeriod(dt));
    EXPECT_EQ(dog.getOrder(), order);
    EXPECT_NEAR(dog.getSamplingPeriod(), dt, 1.e-14);
    EXPECT_NEAR(dog.computeConeOfInfluenceScalar(), 2.3748208234474517, 1.e-14);

    std::vector<double> kwave(nk, 0);
    kwave[0] = 0;
    for (int i=1; i<nk; ++i) 
    {
        kwave[i] = i*df;
    }
    std::vector<std::complex<double>> w(nk);
    auto wPtr = w.data();
    dog.evaluate(nk, scale, kwave.data(), &wPtr);
 
    EXPECT_EQ(w.size(), refWavelet1.size());
    double emax = 0;
    for (int i=0; i<static_cast<int> (w.size()); ++i)
    {
        emax = std::max(emax, std::abs(refWavelet1[i] - w[i]));
        // std::cout << i << " " << kwave[i] << " " << w[i] <<  std::endl;
    }
    EXPECT_NEAR(emax, 0, 1.e-14);

    Wavelets::DerivativeOfGaussian dogCopy(dog);
    EXPECT_EQ(dogCopy.getOrder(), order);
    EXPECT_NEAR(dogCopy.getSamplingPeriod(), dt, 1.e-14);
    dogCopy.evaluate(nk, scale, kwave.data(), &wPtr); 
    emax = 0;
    for (int i=0; i<static_cast<int> (w.size()); ++i)
    {
        emax = std::max(emax, std::abs(refWavelet1[i] - w[i]));
    }
    EXPECT_NEAR(emax, 0, 1.e-14);
    EXPECT_NEAR(dogCopy.computeConeOfInfluenceScalar(),
                2.3748208234474517, 1.e-14);
}
*/

/*
TEST(UtilitiesTransformsWavelets, ricker)
{
    std::vector<std::complex<double>> refWavelet1({ 
        std::complex<double> (0.0000000000000000,0.0000000000000000),
        std::complex<double> (1.6468690890961695,0.0000000000000000),
        std::complex<double> (5.5934887331000285,0.0000000000000000),
        std::complex<double> (9.5823371962771340,0.0000000000000000),
        std::complex<double> (11.630477664277823,0.0000000000000000),
        std::complex<double> (11.125221645150589,0.0000000000000000),
        std::complex<double> (8.7943698103045662,0.0000000000000000),
        std::complex<double> (5.8921572491094105,0.0000000000000000),
        std::complex<double> (3.3968526905871901,0.0000000000000000),
        std::complex<double> (1.7015421611687627,0.0000000000000000),
        std::complex<double> (0.74552312199243487,0.0000000000000000),
        std::complex<double> (0.28707301784948097,0.0000000000000000),
        std::complex<double> (9.74895013848125769E-002,0.0000000000000000),
        std::complex<double> (2.92760769291654788E-002,0.0000000000000000),
        std::complex<double> (7.79032417022442695E-003,0.0000000000000000),
        std::complex<double> (1.83991485371643175E-003,0.0000000000000000),
        std::complex<double> (3.86200993813166969E-004,0.0000000000000000),
        std::complex<double> (7.21226946072720622E-005,0.0000000000000000),
        std::complex<double> (1.19939506305377826E-005,0.0000000000000000),
        std::complex<double> (1.77750184094896598E-006,0.0000000000000000)
    });
    double dt = 0.25;
    int nk = 20;
    double scale = 13.454342644059432;
    int npad = 1024;
    double df = (2*M_PI)/(dt*npad);
    // Initialize the wavelet
    Wavelets::Ricker ricker;
    EXPECT_NO_THROW(ricker.setSamplingPeriod(dt));
    EXPECT_NEAR(ricker.getSamplingPeriod(), dt, 1.e-14);
    EXPECT_NEAR(ricker.computeConeOfInfluenceScalar(),
                2.8099258924162904, 1.e-14);

    std::vector<double> kwave(nk, 0);
    kwave[0] = 0;
    for (int i=1; i<nk; ++i)
    {
        kwave[i] = i*df;
    }
    std::vector<std::complex<double>> w(nk);
    auto wPtr = w.data();
    ricker.evaluate(nk, scale, kwave.data(), &wPtr);

    EXPECT_EQ(w.size(), refWavelet1.size());
    double emax = 0;
    for (int i=0; i<static_cast<int> (w.size()); ++i)
    {
        emax = std::max(emax, std::abs(refWavelet1[i] - w[i]));
    }
    EXPECT_NEAR(emax, 0, 1.e-14);

    Wavelets::Ricker rickerCopy(ricker);
    EXPECT_NEAR(rickerCopy.getSamplingPeriod(), dt, 1.e-14);
    rickerCopy.evaluate(nk, scale, kwave.data(), &wPtr);
    emax = 0;
    for (int i=0; i<static_cast<int> (w.size()); ++i)
    {
        emax = std::max(emax, std::abs(refWavelet1[i] - w[i]));
    }
    EXPECT_NEAR(emax, 0, 1.e-14);
    EXPECT_NEAR(rickerCopy.computeConeOfInfluenceScalar(),
                2.8099258924162904, 1.e-14);

}
*/

TEST(UtilitiesTransformsWavelets, morlet)
{
    // from scipy.signal import morlet2
    // morlet2(25, s = 4, w=8)
    const std::vector<std::complex<double>> wRef({
       std::complex<double> (0.0017697280686101674,+0.0037781866095412506),
       std::complex<double> (-0.008560310412059065,+7.577292337759703e-05),
       std::complex<double> (0.006733793359028701, -0.015064579509298948),
       std::complex<double> (0.019729992282475312, +0.02243919076717591),
       std::complex<double> (-0.04867485704121378, +0.014633231402493077),
       std::complex<double> (0.011105953507896739, -0.08045826426135003),
       std::complex<double> (0.10288890284209785,  +0.06542293042124467),
       std::complex<double> (-0.14427429382906345, +0.093541800528087),
       std::complex<double> (-0.03314350159009988, -0.22536624742831696),
       std::complex<double> (0.27219834187680325,  +0.07921140277081239),
       std::complex<double> (-0.21663905522635993, +0.25082930872918496),
       std::complex<double> (-0.1514807445857645,  -0.3309914654364319),
       std::complex<double> (0.37556277223247125,  +0),
       std::complex<double> (-0.1514807445857645,  +0.3309914654364319),
       std::complex<double> (-0.21663905522635993, -0.25082930872918496),
       std::complex<double> (0.27219834187680325,  -0.07921140277081239),
       std::complex<double> (-0.03314350159009988, +0.22536624742831696),
       std::complex<double> (-0.14427429382906345, -0.093541800528087),
       std::complex<double> (0.10288890284209785,  -0.06542293042124467),
       std::complex<double> (0.011105953507896739, +0.08045826426135003),
       std::complex<double> (-0.04867485704121378, -0.014633231402493077),
       std::complex<double> (0.019729992282475312, -0.02243919076717591),
       std::complex<double> (0.006733793359028701, +0.015064579509298948),
       std::complex<double> (-0.008560310412059065,-7.577292337759703e-05),
       std::complex<double> (0.0017697280686101674,-0.0037781866095412506)
    });
    // morlet2(8, s = 4, w=8)
    const std::vector<std::complex<double>> wRef8({
       std::complex<double> (0.19308308170504393,-0.16826186205005536),
       std::complex<double> (0.08763161987861195,+0.2962400060213443),
       std::complex<double> (-0.3465597394124581,-0.04940089282196845),
       std::complex<double> (0.201338315852819,-0.31356584837818885),
       std::complex<double> (0.201338315852819,+0.31356584837818885),
       std::complex<double> (-0.3465597394124581,+0.04940089282196845),
       std::complex<double> (0.08763161987861195,-0.2962400060213443),
       std::complex<double> (0.19308308170504393,+0.16826186205005536)
    });
    Wavelets::Morlet morlet;
    double omega0 = 8;
    const double scale = 4;
    EXPECT_NO_THROW(morlet.setParameter(omega0));
    EXPECT_NEAR(morlet.getParameter(), omega0, 1.e-14);
    //EXPECT_NEAR(morlet.computeConeOfInfluenceScalar(),
    //            0.55108811163349181, 1.e-14);
    auto n = static_cast<int> (wRef.size());
    std::vector<std::complex<double>> daughter(n);
    auto dPtr = daughter.data();
    morlet.evaluate(n, scale, &dPtr);
    double error = 0;
    for (int i=0; i<n; ++i)
    {
        // Normalize to conform with scipy
        auto xnorm = 1./std::sqrt(scale);
        error = std::max(error, std::abs(xnorm*daughter[i] - wRef[i]));
    }
    EXPECT_NEAR(error, 0, 1.e-14);
    
    Wavelets::Morlet mcopy(morlet);
    EXPECT_NEAR(mcopy.getParameter(), omega0, 1.e-14); 
    n = static_cast<int> (wRef8.size());
    dPtr = daughter.data();
    mcopy.evaluate(n, scale, &dPtr);
    error = 0;
    for (int i=0; i<n; ++i)
    {
        // Normalize to conform with scipy
        auto xnorm = 1./std::sqrt(scale);
        error = std::max(error, std::abs(xnorm*daughter[i] - wRef8[i]));
    }
    EXPECT_NEAR(error, 0, 1.e-14);
}


}
