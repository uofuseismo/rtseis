#include <vector>
#include "rtseis/vector.hpp"
#include "rtseis/demean.hpp"
#include <gtest/gtest.h>

namespace
{
using namespace RTSeis;

using MyTypes = ::testing::Types<double, float>;

template<class T>
class DemeanTest : public testing::Test
{
public:
    std::vector<T> inputSignal{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<T> outputSignal{-4.5, -3.5, -2.5, -1.5, -0.5,
                                 0.5,  1.5,  2.5,  3.5,  4.5};
    const double mean{5.5};
protected:
    DemeanTest() :
        x{inputSignal},
        y{outputSignal}
    {
    }
    ~DemeanTest() = default;
public:
    Demean<T> demean;
    Vector<T> x;
    Vector<T> y;
    const T epsilon{std::numeric_limits<T>::epsilon()*10};
};

TYPED_TEST_SUITE(DemeanTest, MyTypes);

TYPED_TEST(DemeanTest, Demean)
{
    auto epsilon = this->epsilon;
    auto x = this->x;
    auto yRef = this->y;
    EXPECT_TRUE(this->demean.isInitialized());
    EXPECT_NO_THROW(this->demean.setInput(x)); 
    EXPECT_NO_THROW(this->demean.apply());
    auto y = this->demean.getOutput();
    EXPECT_EQ(yRef.size(), y.size());
    auto yRefPtr = yRef.data();
    auto yPtr = y.data();
    for (int i = 0; i < static_cast<int> (yRef.size()); ++i)
    {
        EXPECT_NEAR(std::abs(yRefPtr[i] - yPtr[i]), 0, epsilon);
    }
    EXPECT_NEAR(std::abs(this->mean - this->demean.getMean()), 0, epsilon);

    auto dCopy = this->demean;
    const auto &yCopy = dCopy.getOutputReference();
    EXPECT_EQ(yCopy.size(), yRef.size());
    for (int i = 0; i < static_cast<int> (yCopy.size()); ++i)
    {
        EXPECT_NEAR(std::abs(yRefPtr[i] - yCopy.data()[i]), 0, epsilon);
    }
    EXPECT_NEAR(std::abs(this->mean - dCopy.getMean()), 0, epsilon);

    auto dMove = std::move(dCopy);
    const auto &yMove = dMove.getOutputReference();
    EXPECT_EQ(yMove.size(), yRef.size());
    EXPECT_NEAR(std::abs(this->mean - dMove.getMean()), 0, epsilon);
}

}
