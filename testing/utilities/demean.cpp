#include <iostream>
#include <vector>
#include "rtseis/vector.hpp"
#include "rtseis/demean.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace RTSeis;

TEMPLATE_TEST_CASE("CoreTest::Demean", "[TypeName][template]", float, double)
{
    std::vector<TestType> inputSignal{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<TestType> outputSignal{-4.5, -3.5, -2.5, -1.5, -0.5,
                                       0.5,  1.5,  2.5,  3.5,  4.5};
    constexpr TestType mean{5.5};
 
    Demean<TestType> demean;   
    Vector<TestType> x{inputSignal};
    Vector<TestType> yRef{outputSignal};
    REQUIRE(demean.isInitialized());
    REQUIRE_NOTHROW(demean.setInput(x));
    REQUIRE_NOTHROW(demean.apply());
    REQUIRE(mean == Catch::Approx(demean.getMean()));
    auto y = demean.getOutput(); 
    REQUIRE(yRef.size() == y.size());
    for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
    {
        CHECK(y.at(i) == Catch::Approx(yRef.at(i)));
    }

    SECTION("copy")
    {
        auto dCopy = demean;
        auto yCopy = dCopy.getOutputReference();
        REQUIRE(yCopy.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(yCopy.at(i) == Catch::Approx(yRef.at(i)));
        }   

    }

    SECTION("move")
    {
        auto dMove = std::move(demean);
        auto yMove = dMove.getOutputReference(); 
        REQUIRE(yMove.size() == yRef.size());
        for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
        {
            CHECK(yMove[i] == Catch::Approx(yRef[i]));
        }
    }
/*
    auto epsilon = this->epsilon;
    auto x = this->x;
    }
    EXPECT_NEAR(std::abs(this->mean - dCopy.getMean()), 0, epsilon);

    auto dMove = std::move(dCopy);
    const auto &yMove = dMove.getOutputReference();
    EXPECT_EQ(yMove.size(), yRef.size());
    EXPECT_NEAR(std::abs(this->mean - dMove.getMean()), 0, epsilon);
*/
}
