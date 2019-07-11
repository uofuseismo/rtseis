#include <mkl.h>
#include <gtest/gtest.h>

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    auto success = RUN_ALL_TESTS();
    mkl_free_buffers();
    return success;
}
