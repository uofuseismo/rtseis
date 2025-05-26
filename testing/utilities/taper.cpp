#include <fstream>
#include <filesystem>
#include <vector>
#include <limits>
#include "rtseis/vector.hpp"
#include "rtseis/taper.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace RTSeis;

std::filesystem::path dataDirectory{"data"};
std::filesystem::path taperSolutions100FileName{dataDirectory/"taper100.all.txt"};
std::filesystem::path taperSolutions101FileName{dataDirectory/"taper101.all.txt"};

template<class T>
class DataLoader
{
public:
    void load(const bool is100 = true)
    {
        std::string line;
        std::ifstream taperFile;
        auto fileName = taperSolutions101FileName;
        if (is100)
        {
            fileName = taperSolutions100FileName;
        }
        if (!std::filesystem::exists(fileName))
        {
            throw std::runtime_error(std::string {fileName} + " does not exist");
        }
        taperFile.open(fileName);
        mHammingReference.clear();
        mHanningReference.clear();
        mSineReference.clear();
        mHammingReference.reserve(101);
        mHanningReference.reserve(101);
        mSineReference.reserve(101);
        int nLines{0};
        while (std::getline(taperFile, line))
        {
            double yHamming, yHanning, ySine;
            std::sscanf(line.c_str(), "%lf, %lf, %lf\n",
                        &yHamming, &yHanning, &ySine);
            mHammingReference.push_back(yHamming);
            mHanningReference.push_back(yHanning);
            mSineReference.push_back(ySine);
            nLines = nLines + 1;
        }
        taperFile.close();
        if (is100)
        {
            if (mHammingReference.size() != 100)
            {
                throw std::runtime_error("Should be 100 lines in file");
            }
        }
        else
        {
            if (mHammingReference.size() != 101)
            {
                throw std::runtime_error("Should be 101 lines in file");
            }
        }
    }
public:
    //Taper<T> taper;
    //Vector<T> ones;
    std::vector<T> mHammingReference;
    std::vector<T> mHanningReference;
    std::vector<T> mSineReference;
    double mHammingPercentage100{40};
    double mHanningPercentage100{20};
    double mSinePercentage100{60};
    double mHammingPercentage101{10};
    double mHanningPercentage101{20};
    double mSinePercentage101{30};
    //const T eps{std::numeric_limits<T>::epsilon()};
    //typename Taper<T>::Window mHamming{Taper<T>::Window::Hamming};
    //typename Taper<T>::Window mHanning{Taper<T>::Window::Hanning};
    //typename Taper<T>::Window mSine{Taper<T>::Window::Sine};
    //typename Taper<T>::Window mBoxcar{Taper<T>::Window::Boxcar};
};

TEMPLATE_TEST_CASE("CoreTest::Taper", "[TypeName][template]", float, double)
{
    DataLoader<TestType> dataLoader100;
    REQUIRE_NOTHROW(dataLoader100.load(true));
    DataLoader<TestType> dataLoader101;
    REQUIRE_NOTHROW(dataLoader101.load(false));

    Taper<TestType> taper;
    SECTION("Hamming100")
    {
        constexpr int nSamples{100}; 
        constexpr TestType one{1};
        constexpr double hammingPercentage100{40};
        auto window = Taper<TestType>::Window::Hamming;
        Vector<TestType> x;
        x.resize(nSamples, one);

        auto yRef = dataLoader100.mHammingReference;
        
        REQUIRE_NOTHROW(taper.initialize(window, hammingPercentage100));
        REQUIRE(taper.isInitialized());
        REQUIRE_NOTHROW(taper.setInput(x));
        REQUIRE_NOTHROW(taper.apply());
        auto y = taper.getOutput(); 
        REQUIRE(yRef.size() == y.size()); 
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("Hanning100")
    {
        constexpr int nSamples{100};
        constexpr TestType one{1};
        constexpr double hanningPercentage100{20};
        auto window = Taper<TestType>::Window::Hanning;
        Vector<TestType> x;
        x.resize(nSamples, one);

        auto yRef = dataLoader100.mHanningReference;

        REQUIRE_NOTHROW(taper.initialize(window, hanningPercentage100));
        REQUIRE(taper.isInitialized());
        REQUIRE_NOTHROW(taper.setInput(x));
        REQUIRE_NOTHROW(taper.apply());
        auto y = taper.getOutput(); 
        REQUIRE(yRef.size() == y.size()); 
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("Sine100")
    {
        constexpr int nSamples{100};
        constexpr TestType one{1};
        constexpr double sinePercentage100{60};
        auto window = Taper<TestType>::Window::Sine;
        Vector<TestType> x;
        x.resize(nSamples, one);

        auto yRef = dataLoader100.mSineReference;

        REQUIRE_NOTHROW(taper.initialize(window, sinePercentage100));
        REQUIRE(taper.isInitialized());
        REQUIRE_NOTHROW(taper.setInput(x));
        REQUIRE_NOTHROW(taper.apply());
        auto y = taper.getOutput(); 
        REQUIRE(yRef.size() == y.size()); 
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }   

    SECTION("Boxcar100")
    {
        constexpr int nSamples{100};
        constexpr TestType one{1};
        constexpr TestType zero{0};
        constexpr double boxcarPercentage100{10};
        auto window = Taper<TestType>::Window::Boxcar;
        Vector<TestType> x{nSamples, one};
        //x.resize(nSamples, one);

        REQUIRE_NOTHROW(taper.initialize(window, boxcarPercentage100));
        REQUIRE(taper.isInitialized());
        REQUIRE_NOTHROW(taper.setInput(x));
        REQUIRE_NOTHROW(taper.apply());
        auto y = taper.getOutput();
        REQUIRE(nSamples == static_cast<int> (y.size()));

        for (int i = 0; i < nSamples; ++i)
        {
            if (i < 5 || i >= 95)
            {
                CHECK(y.at(i) == Catch::Approx(zero));
            }
            else
            {
                CHECK(y.at(i) == Catch::Approx(one));
            }
        }
    }

    SECTION("Hamming101")
    {
        constexpr int nSamples{101}; 
        constexpr TestType one{1};
        constexpr double hammingPercentage101{10};
        auto window = Taper<TestType>::Window::Hamming;
        Vector<TestType> x;
        x.resize(nSamples, one);

        auto yRef = dataLoader101.mHammingReference;
    
        REQUIRE_NOTHROW(taper.initialize(window, hammingPercentage101));
        REQUIRE(taper.isInitialized());
        REQUIRE_NOTHROW(taper.setInput(x));
        REQUIRE_NOTHROW(taper.apply());
        auto y = taper.getOutputReference(); 
        REQUIRE(yRef.size() == y.size()); 
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("Hanning101")
    {
        constexpr int nSamples{101};
        constexpr TestType one{1};
        constexpr double hanningPercentage101{20};
        auto window = Taper<TestType>::Window::Hanning;
        Vector<TestType> x;
        x.resize(nSamples, one);

        auto yRef = dataLoader101.mHanningReference;

        REQUIRE_NOTHROW(taper.initialize(window, hanningPercentage101));
        REQUIRE(taper.isInitialized());
        REQUIRE_NOTHROW(taper.setInput(x));
        REQUIRE_NOTHROW(taper.apply());
        auto y = taper.getOutput(); 
        REQUIRE(yRef.size() == y.size()); 
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {   
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }

    SECTION("Sine101")
    {
        constexpr int nSamples{101};
        constexpr TestType one{1};
        constexpr double sinePercentage101{30};
        auto window = Taper<TestType>::Window::Sine;
        Vector<TestType> x;
        x.resize(nSamples, one);

        auto yRef = dataLoader101.mSineReference;

        REQUIRE_NOTHROW(taper.initialize(window, sinePercentage101));
        REQUIRE(taper.isInitialized());
        REQUIRE_NOTHROW(taper.setInput(x));
        REQUIRE_NOTHROW(taper.apply());
        auto y = taper.getOutput(); 
        REQUIRE(yRef.size() == y.size()); 
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
        }
    }
}
/*
TYPED_TEST(TaperTest, HammingTaper100)
{
    constexpr bool is100{true};
    this->load(is100);
    EXPECT_NO_THROW(this->taper.initialize(this->mHamming, this->mHammingPercentage100));
    EXPECT_TRUE(this->taper.isInitialized());
    EXPECT_EQ(this->taper.getWindow(), this->mHamming);
    EXPECT_NEAR(this->taper.getPercentage(), this->mHammingPercentage100, 1.e-10);
    EXPECT_NO_THROW(this->taper.setInput(this->ones));
    EXPECT_NO_THROW(this->taper.apply());
    auto y = this->taper.getOutput();
    auto yReference = this->mHammingReference;
    EXPECT_EQ(yReference.size(), y.size());
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {
         EXPECT_NEAR(y[i], yReference[i], 1.e-6);
    } 
}

TYPED_TEST(TaperTest, HanningTaper100)
{
    constexpr bool is100{false};
    this->load(is100);
    EXPECT_NO_THROW(this->taper.initialize(this->mHanning, this->mHanningPercentage100));
    EXPECT_TRUE(this->taper.isInitialized());
    EXPECT_EQ(this->taper.getWindow(), this->mHanning);
    EXPECT_NEAR(this->taper.getPercentage(), this->mHanningPercentage100, 1.e-10);
    EXPECT_NO_THROW(this->taper.setInput(this->ones));
    EXPECT_NO_THROW(this->taper.apply());
    auto y = this->taper.getOutput();
    auto yReference = this->mHanningReference;
    EXPECT_EQ(yReference.size(), y.size());
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {   
         EXPECT_NEAR(y[i], yReference[i], 1.e-6);
    }   
}

TYPED_TEST(TaperTest, SineTaper100)
{
    constexpr bool is100{false};
    this->load(is100);
    EXPECT_NO_THROW(this->taper.initialize(this->mSine, this->mSinePercentage100));
    EXPECT_TRUE(this->taper.isInitialized());
    EXPECT_EQ(this->taper.getWindow(), this->mSine);
    EXPECT_NEAR(this->taper.getPercentage(), this->mSinePercentage100, 1.e-10);
    EXPECT_NO_THROW(this->taper.setInput(this->ones));
    EXPECT_NO_THROW(this->taper.apply());
    auto y = this->taper.getOutput();
    auto yReference = this->mSineReference;
    EXPECT_EQ(yReference.size(), y.size());
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {   
         EXPECT_NEAR(y[i], yReference[i], 1.e-6);
    }   
}

TYPED_TEST(TaperTest, HammingTaper101)
{
    constexpr bool is100{false};
    this->load(is100);
    EXPECT_NO_THROW(this->taper.initialize(this->mHamming, this->mHammingPercentage101));
    EXPECT_TRUE(this->taper.isInitialized());
    EXPECT_EQ(this->taper.getWindow(), this->mHamming);
    EXPECT_NEAR(this->taper.getPercentage(), this->mHammingPercentage101, 1.e-10);
    EXPECT_NO_THROW(this->taper.setInput(this->ones));
    EXPECT_NO_THROW(this->taper.apply());
    auto y = this->taper.getOutput();
    auto yReference = this->mHammingReference;
    EXPECT_EQ(yReference.size(), y.size());
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {   
         EXPECT_NEAR(y[i], yReference[i], 1.e-6);
    }   
}

TYPED_TEST(TaperTest, HanningTaper101)
{
    constexpr bool is100{false};
    this->load(is100);
    EXPECT_NO_THROW(this->taper.initialize(this->mHanning, this->mHanningPercentage101));
    EXPECT_TRUE(this->taper.isInitialized());
    EXPECT_EQ(this->taper.getWindow(), this->mHanning);
    EXPECT_NEAR(this->taper.getPercentage(), this->mHanningPercentage101, 1.e-10);
    EXPECT_NO_THROW(this->taper.setInput(this->ones));
    EXPECT_NO_THROW(this->taper.apply());
    auto y = this->taper.getOutput();
    auto yReference = this->mHanningReference;
    EXPECT_EQ(yReference.size(), y.size());
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {   
         EXPECT_NEAR(y[i], yReference[i], 1.e-6);
    }   
}

TYPED_TEST(TaperTest, SineTaper101)
{
    constexpr bool is100{false};
    this->load(is100);
    EXPECT_NO_THROW(this->taper.initialize(this->mSine, this->mSinePercentage101));
    EXPECT_TRUE(this->taper.isInitialized());
    EXPECT_EQ(this->taper.getWindow(), this->mSine);
    EXPECT_NEAR(this->taper.getPercentage(), this->mSinePercentage101, 1.e-10);
    EXPECT_NO_THROW(this->taper.setInput(this->ones));
    EXPECT_NO_THROW(this->taper.apply());
    auto y = this->taper.getOutput();
    auto yReference = this->mSineReference;
    EXPECT_EQ(yReference.size(), y.size());
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {   
         EXPECT_NEAR(y[i], yReference[i], 1.e-6);
    }   
}

TYPED_TEST(TaperTest, Boxcar100)
{
    constexpr bool is100{false};
    constexpr double percentage{10};
    this->ones.resize(100, 1);
    EXPECT_NO_THROW(this->taper.initialize(this->mBoxcar, percentage));
    EXPECT_TRUE(this->taper.isInitialized());
    EXPECT_EQ(this->taper.getWindow(), this->mBoxcar);
    EXPECT_NEAR(this->taper.getPercentage(), percentage, 1.e-10);
    EXPECT_NO_THROW(this->taper.setInput(this->ones));
    EXPECT_NO_THROW(this->taper.apply());
    auto y = this->taper.getOutput();
    EXPECT_EQ(this->ones.size(), y.size());
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {
        if (i < 5 || i >= 95)
        {
            EXPECT_NEAR(y[i], 0, this->eps);
        }
        else
        {
            EXPECT_NEAR(y[i], 1, this->eps);
        } 
    }   
}
*/
