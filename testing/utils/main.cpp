#ifdef WITH_MKL
#include <mkl.h>
#endif
#include <gtest/gtest.h>

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    auto success = RUN_ALL_TESTS();
#ifdef WITH_MKL
    mkl_free_buffers();
#endif
    return success;
}
