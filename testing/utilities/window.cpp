#include <vector>
#include <limits>
#include "rtseis/vector.hpp"
#include "rtseis/window.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace RTSeis;

using MyTypes = ::testing::Types<double, float>;

template<class T>
class WindowTest : public testing::Test
{
public:
    std::vector<T> bartlett20{0,
                              0.105263157894737, 0.210526315789474, 0.315789473684211,
                              0.421052631578947, 0.526315789473684, 0.631578947368421,
                              0.736842105263158, 0.842105263157895, 0.947368421052632,
                              0.947368421052632, 0.842105263157895, 0.736842105263158,
                              0.631578947368421, 0.526315789473684, 0.421052631578947,
                              0.315789473684211, 0.210526315789474, 0.105263157894737,
                              0};
    std::vector<T> bartlett19{0,
                              0.111111111111111, 0.222222222222222, 0.333333333333333,
                              0.444444444444444, 0.555555555555556, 0.666666666666667,
                              0.777777777777778, 0.888888888888889, 1.000000000000000,
                              0.888888888888889, 0.777777777777778, 0.666666666666667,
                              0.555555555555556, 0.444444444444444, 0.333333333333333,
                              0.222222222222222, 0.111111111111111,
                              0};
    std::vector<T> blackman20{0,
                              0.010222619901394, 0.045068584273067, 0.114390286966549,
                              0.226899356333081, 0.382380768463948, 0.566665186596425,
                              0.752034438175084, 0.903492728253039, 0.988846031037412,
                              0.988846031037412, 0.903492728253039, 0.752034438175084,
                              0.566665186596425, 0.382380768463948, 0.226899356333081,
                              0.114390286966549, 0.045068584273067, 0.010222619901394,
                              0};
    std::vector<T> blackman19{0, 
                              0.011437245056564, 0.050869632653865, 0.130000000000000,
                              0.258000501503662, 0.431648679170593, 0.630000000000000,
                              0.816914075772843, 0.951129865842472, 1.000000000000000,
                              0.951129865842472, 0.816914075772843, 0.630000000000000,
                              0.431648679170593, 0.258000501503662, 0.130000000000000,
                              0.050869632653865, 0.011437245056564,
                              0};

    std::vector<T> hamming20{0.080000000000000,
                             0.104924068817708, 0.176995365677659, 0.288403847263684,
                             0.427076675915232, 0.577986498917273, 0.724779895340366,
                             0.851549522947841, 0.944557925554985, 0.993726199565252,
                             0.993726199565252, 0.944557925554985, 0.851549522947841,
                             0.724779895340366, 0.577986498917273, 0.427076675915232,
                             0.288403847263684, 0.176995365677659, 0.104924068817708,
                             0.080000000000000 };
    std::vector<T> hamming19{0.080000000000000,
                             0.107741394438482, 0.187619556165270, 0.310000000000000,
                             0.460121838273212, 0.619878161726788, 0.770000000000000,
                             0.892380443834730, 0.972258605561518, 1.000000000000000,
                             0.972258605561518, 0.892380443834730, 0.770000000000000,
                             0.619878161726788, 0.460121838273212, 0.310000000000000,
                             0.187619556165270, 0.107741394438482,
                             0.080000000000000};
    std::vector<T> hanning20{0,
                             0.027091379149683, 0.105429745301803, 0.226525920938787,
                             0.377257256429600, 0.541289672736166, 0.700847712326485,
                             0.838640785812870, 0.939736875603244, 0.993180651701361,
                             0.993180651701361, 0.939736875603244, 0.838640785812870,
                             0.700847712326485, 0.541289672736166, 0.377257256429600,
                             0.226525920938787, 0.105429745301803, 0.027091379149683,
                             0};
    std::vector<T> hanning19{0,
                             0.030153689607046, 0.116977778440511, 0.250000000000000,
                             0.413175911166535, 0.586824088833465, 0.750000000000000,
                             0.883022221559489, 0.969846310392954, 1.000000000000000,
                             0.969846310392954, 0.883022221559489, 0.750000000000000,
                             0.586824088833465, 0.413175911166535, 0.250000000000000,
                             0.116977778440511, 0.030153689607046,
                             0};
    // w20 = sin(linspace(0,N-1,N)*pi/(N-1))
    std::vector<T> sine20{0.000000000000000,
                          0.164594590280734, 0.324699469204683, 0.475947393037074,
                          0.614212712689668, 0.735723910673132, 0.837166478262529,
                          0.915773326655057, 0.969400265939330, 0.996584493006670,
                          0.996584493006670, 0.969400265939330, 0.915773326655057,
                          0.837166478262528, 0.735723910673132, 0.614212712689668,
                          0.475947393037074, 0.324699469204683, 0.164594590280734,
                          0.000000000000000};
    // w19 = sin(linspace(0,N-1,N)*pi/(N-1))
    std::vector<T> sine19{0.000000000000000,
                          0.173648177666930, 0.342020143325669, 0.500000000000000,
                          0.642787609686539, 0.766044443118978, 0.866025403784439,
                          0.939692620785908, 0.984807753012208, 1.000000000000000,
                          0.984807753012208, 0.939692620785908, 0.866025403784439,
                          0.766044443118978, 0.642787609686539, 0.500000000000000,
                          0.342020143325669, 0.173648177666930,
                          0.000000000000000};
    // beta = 2.5
    std::vector<T> kaiser20{0.303966229415369, 0.406333120777527,
                            0.511011861080447, 0.614155742009957,
                            0.711825746588635, 0.800182020135053,
                            0.875673636668776, 0.935216242888425,
                            0.976347754218060, 0.997353444300571,
                            0.997353444300571, 0.976347754218060,
                            0.935216242888425, 0.875673636668776,
                            0.800182020135053, 0.711825746588635,
                            0.614155742009957, 0.511011861080447,
                            0.406333120777527, 0.303966229415369};
    // beta = 5.5
    std::vector<T> kaiser19{0.023422141030647, 0.078225306346850,
                            0.166678217276960, 0.288508741610676,
                            0.436709600995736, 0.597737293640817,
                            0.753274215094738, 0.883279467813599,
                            0.969699049461490, 1.000000000000000,
                            0.969699049461490, 0.883279467813599,
                            0.753274215094738, 0.597737293640817,
                            0.436709600995736, 0.288508741610676,
                            0.166678217276960, 0.078225306346850,
                            0.023422141030647};

protected:
    WindowTest()
    {   
        x19.resize(19, 1);
        x20.resize(20, 1);
    }   
    ~WindowTest() = default;
public:
    Window<T> window;
    Vector<T> x19;
    Vector<T> x20;
    Vector<T> y;
    const T epsilon{std::numeric_limits<T>::epsilon()*10};
    typename Window<T>::Type mBartlett{Window<T>::Type::Bartlett};
    typename Window<T>::Type mBlackman{Window<T>::Type::Blackman};
    typename Window<T>::Type mSine{Window<T>::Type::Sine};
    typename Window<T>::Type mHanning{Window<T>::Type::Hanning};
    typename Window<T>::Type mHamming{Window<T>::Type::Hamming};
    typename Window<T>::Type mKaiser{Window<T>::Type::Kaiser};
};

TYPED_TEST_SUITE(WindowTest, MyTypes);

TYPED_TEST(WindowTest, SineDesignWindow)
{
    const auto type = this->mSine;
    auto epsilon = this->epsilon;
    auto x19 = this->x19;
    auto x20 = this->x20;
    auto yRef19 = this->sine19;
    auto yRef20 = this->sine20;
    this->window.initialize(yRef19.size(), this->mSine);
    EXPECT_TRUE(this->window.isInitialized());
    EXPECT_EQ(this->window.getType(), type);
    auto window19 = this->window;
    const auto &y = window19.getWindowReference();
    EXPECT_EQ(y.size(), yRef19.size());
    for (int i = 0; i < y.size(); ++i)
    {
        EXPECT_NEAR(std::abs(y.at(i) - yRef19.at(i)), 0, epsilon);
    }
 
    EXPECT_NO_THROW(this->window.initialize(yRef20.size(), type));
    auto y20 = this->window.getWindowReference(); 
    EXPECT_EQ(y20.size(), yRef20.size());
    for (int i = 0; i < y20.size(); ++i)
    {
        EXPECT_NEAR(std::abs(y20.at(i) - yRef20.at(i)), 0, epsilon);
    }
}

TYPED_TEST(WindowTest, HanningDesignWindow)
{
    const auto type = this->mHanning;
    auto epsilon = this->epsilon;
    auto x19 = this->x19;
    auto x20 = this->x20;
    auto yRef19 = this->hanning19;
    auto yRef20 = this->hanning20;
    this->window.initialize(yRef19.size(), type);
    EXPECT_TRUE(this->window.isInitialized());
    EXPECT_EQ(this->window.getType(), type);
    auto y = this->window.getWindow();
    EXPECT_EQ(y.size(), yRef19.size());
    for (int i = 0; i < y.size(); ++i)
    {   
        EXPECT_NEAR(std::abs(y.at(i) - yRef19.at(i)), 0, epsilon);
    }   
 
    EXPECT_NO_THROW(this->window.initialize(yRef20.size(), type));
    y = this->window.getWindowReference(); 
    EXPECT_EQ(y.size(), yRef20.size());
    for (int i = 0; i < y.size(); ++i)
    {   
        EXPECT_NEAR(std::abs(y.at(i) - yRef20.at(i)), 0, epsilon);
    }   
}

TYPED_TEST(WindowTest, HammingDesignWindow)
{
    const auto type = this->mHamming;
    auto epsilon = this->epsilon;
    auto x19 = this->x19;
    auto x20 = this->x20;
    auto yRef19 = this->hamming19;
    auto yRef20 = this->hamming20;
    this->window.initialize(yRef19.size(), type);
    EXPECT_TRUE(this->window.isInitialized());
    EXPECT_EQ(this->window.getType(), type);
    auto y = this->window.getWindow();
    EXPECT_EQ(y.size(), yRef19.size());
    for (int i = 0; i < y.size(); ++i)
    {   
        EXPECT_NEAR(std::abs(y.at(i) - yRef19.at(i)), 0, epsilon);
    }   
 
    EXPECT_NO_THROW(this->window.initialize(yRef20.size(), type));
    y = this->window.getWindowReference(); 
    EXPECT_EQ(y.size(), yRef20.size());
    for (int i = 0; i < y.size(); ++i)
    {   
        EXPECT_NEAR(std::abs(y.at(i) - yRef20.at(i)), 0, epsilon);
    }
}

TYPED_TEST(WindowTest, BlackmanDesignWindow)
{
    const auto type = this->mBlackman;
    auto epsilon = this->epsilon;
    auto x19 = this->x19;
    auto x20 = this->x20;
    auto yRef19 = this->blackman19;
    auto yRef20 = this->blackman20;
    this->window.initialize(yRef19.size(), type);
    EXPECT_TRUE(this->window.isInitialized());
    EXPECT_EQ(this->window.getType(), type);
    auto y = this->window.getWindow();
    EXPECT_EQ(y.size(), yRef19.size());
    for (int i = 0; i < y.size(); ++i)
    {
        EXPECT_NEAR(std::abs(y.at(i) - yRef19.at(i)), 0, epsilon);
    }

    EXPECT_NO_THROW(this->window.initialize(yRef20.size(), type));
    y = this->window.getWindowReference();
    EXPECT_EQ(y.size(), yRef20.size());
    for (int i = 0; i < y.size(); ++i)
    {
        EXPECT_NEAR(std::abs(y.at(i) - yRef20.at(i)), 0, epsilon);
    }
}

TYPED_TEST(WindowTest, BartlettDesignWindow)
{
    const auto type = this->mBartlett;
    auto epsilon = this->epsilon;
    auto x19 = this->x19;
    auto x20 = this->x20;
    auto yRef19 = this->bartlett19;
    auto yRef20 = this->bartlett20;
    this->window.initialize(yRef19.size(), type);
    EXPECT_TRUE(this->window.isInitialized());
    EXPECT_EQ(this->window.getType(), type);
    auto y = this->window.getWindow();
    EXPECT_EQ(y.size(), yRef19.size());
    for (int i = 0; i < y.size(); ++i)
    {   
        EXPECT_NEAR(std::abs(y.at(i) - yRef19.at(i)), 0, epsilon);
    }   
 
    EXPECT_NO_THROW(this->window.initialize(yRef20.size(), type));
    y = this->window.getWindowReference(); 
    EXPECT_EQ(y.size(), yRef20.size());
    for (int i = 0; i < y.size(); ++i)
    {   
        EXPECT_NEAR(std::abs(y.at(i) - yRef20.at(i)), 0, epsilon);
    }   
}

TYPED_TEST(WindowTest, KaiserDesignWindow)
{
    const auto type = this->mKaiser;
    auto epsilon = this->epsilon;
    auto x19 = this->x19;
    auto x20 = this->x20;
    auto yRef19 = this->kaiser19;
    auto yRef20 = this->kaiser20;
    this->window.initialize(yRef19.size(), type, 5.5);
    EXPECT_TRUE(this->window.isInitialized());
    EXPECT_EQ(this->window.getType(), type);
    auto y = this->window.getWindow();
    EXPECT_EQ(y.size(), yRef19.size());
    for (int i = 0; i < y.size(); ++i)
    {   
        EXPECT_NEAR(std::abs(y.at(i) - yRef19.at(i)), 0, epsilon);
    }   
 
    EXPECT_NO_THROW(this->window.initialize(yRef20.size(), type, 2.5));
    y = this->window.getWindowReference(); 
    EXPECT_EQ(y.size(), yRef20.size());
    for (int i = 0; i < y.size(); ++i)
    {   
        EXPECT_NEAR(std::abs(y.at(i) - yRef20.at(i)), 0, epsilon);
    }   
}

}
