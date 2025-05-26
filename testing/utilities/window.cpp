#include <vector>
#include <limits>
#include "rtseis/vector.hpp"
#include "rtseis/window.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace RTSeis;

template<class T>
class WindowReference
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
};

TEMPLATE_TEST_CASE("CoreTest::Window", "[TypeName][template]", float, double)
{
    WindowReference<TestType> windowReference;

    Window<TestType> window;

    SECTION("BartlettSize19")
    {   
        auto yRef = windowReference.bartlett19;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Bartlett));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Bartlett); 
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size()); 
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }

        SECTION("Copy")
        {
            auto windowCopy = window;
            auto yCopy = windowCopy.getWindowReference();
            REQUIRE(yCopy.size() == yRef.size());
            for (int i = 0; i < static_cast<int> (yCopy.size()); ++i)
            {
                CHECK(yCopy.at(i) == Catch::Approx(yRef.at(i)));
            }
        }
    }

    SECTION("BartlettSize20")
    {
        auto yRef = windowReference.bartlett20;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Bartlett));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Bartlett);
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("BlackmanSize19")
    {
        auto yRef = windowReference.blackman19;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Blackman));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Blackman); 
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size()); 
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }   
    }

    SECTION("BlackmanSize20")
    {
        auto yRef = windowReference.blackman20;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Blackman));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Blackman); 
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size()); 
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {   
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }   
    }

    SECTION("HammingSize19")
    {   
        auto yRef = windowReference.hamming19;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Hamming));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Hamming); 
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size()); 
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }   

    SECTION("HammingSize20")
    {   
        auto yRef = windowReference.hamming20;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Hamming));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Hamming); 
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size()); 
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("HanningSize19")
    {   
        auto yRef = windowReference.hanning19;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Hanning));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Hanning);
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size()); 
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }   

    SECTION("HanningSize20")
    {   
        auto yRef = windowReference.hanning20;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Hanning));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Hanning);
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("KaiserSize19")
    {   
        auto yRef = windowReference.kaiser19;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Kaiser, 5.5));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Kaiser);
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size()); 
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }   

    SECTION("KaiserSize20")
    {   
        auto yRef = windowReference.kaiser20;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Kaiser, 2.5));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Kaiser);
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("SineSize19")
    {
        auto yRef = windowReference.sine19;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Sine));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Sine);
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size()); 
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {   
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("SineSize20")
    {
        auto yRef = windowReference.sine20;
        REQUIRE_NOTHROW(window.initialize(yRef.size(), WindowType::Sine));
        REQUIRE(window.isInitialized());
        REQUIRE(window.getType() == WindowType::Sine);
        auto y = window.getWindowReference();
        REQUIRE(y.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }
}

