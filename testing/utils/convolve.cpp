#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <ipps.h>
#include "rtseis/utilities/math/convolve.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis::Utilities::Math;

/*
int rtseis_test_utils_convolve(void)
{
    int ierr;
    ierr = convolve_convolve_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed convolution test");
        return -1;
    }
    ierr = convolve_correlate_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed correlation test");
        return -1;
    }
    ierr = convolve_autocorrelate_test();
    if (ierr != 0)
    {
        RTSEIS_ERRMSG("%s", "Failed autocorrelation test");
        return -1;
    }
    return EXIT_SUCCESS;
}
*/
//============================================================================//
//int convolve_convolve_test(void)
TEST(UtilitiesConvolve, Convolve)
{
    for (auto imp=0; imp<2; imp++)
    {
        Convolve::Implementation implementation;
        if (imp == 0)
        {
            implementation = Convolve::Implementation::DIRECT;
        }
        else
        {
            implementation = Convolve::Implementation::FFT;
        }
        std::vector<double> a1({1, 3, 2}); 
        std::vector<double> b1({-1, 1, 0.5});
        std::vector<double> r1({-1, -2,   1.5,   3.5,   1});
        std::vector<double> c;
        // Test 1
        EXPECT_NO_THROW(c = Convolve::convolve(a1, b1,
                                               Convolve::Mode::FULL,
                                               implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 5);
        double emax;
        ippsNormDiff_Inf_64f(c.data(), r1.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 2 - interchange 1
        EXPECT_NO_THROW(c = Convolve::convolve(b1, a1,
                                   Convolve::Mode::FULL, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 5);
        ippsNormDiff_Inf_64f(c.data(), r1.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 3
        EXPECT_NO_THROW(c = Convolve::convolve(a1, b1,
                                   Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 1);
        EXPECT_NEAR(r1[2], c[0], 1.e-10);
        // Test 4 - interchange 3
        EXPECT_NO_THROW(
            c = Convolve::convolve(b1, a1,
                                   Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 1);
        EXPECT_NEAR(r1[2], c[0], 1.e-10);
        // Test 5
        EXPECT_NO_THROW(
            c = Convolve::convolve(a1, b1,
                                   Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 3);
        ippsNormDiff_Inf_64f(c.data(), r1.data()+1, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 6 - interchange 5
        EXPECT_NO_THROW(
            c = Convolve::convolve(b1, a1,
                                   Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 3);
        ippsNormDiff_Inf_64f(c.data(), r1.data()+1, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        //--------------------------------------------------------------------//
        //                            Unequal Lengths                         //
        //--------------------------------------------------------------------//
        std::vector<double> a2({1, 3, 2, -0.2, 0.1, 4.0});
        std::vector<double> b2({-2, 2, 2.5});
        std::vector<double> r2({-2, -4, 4.5, 11.9, 4.4, -8.3, 8.25, 10});
        // Test 7
        EXPECT_NO_THROW(
            c = Convolve::convolve(a2, b2,
                                   Convolve::Mode::FULL, implementation));
        EXPECT_EQ(c.size(), r2.size());
        ippsNormDiff_Inf_64f(c.data(), r2.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10); 
        // Test 8 - interchange 7
        EXPECT_NO_THROW(
            c = Convolve::convolve(b2, a2,
                                   Convolve::Mode::FULL, implementation));
        EXPECT_EQ(c.size(), r2.size());
        ippsNormDiff_Inf_64f(c.data(), r2.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 9
        EXPECT_NO_THROW(
            c = Convolve::convolve(a2, b2,
                                   Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 4);
        ippsNormDiff_Inf_64f(c.data(), r2.data()+2, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 10 - interchange 9
        EXPECT_NO_THROW(
            c = Convolve::convolve(b2, a2,
                                   Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 4);
        ippsNormDiff_Inf_64f(c.data(), r2.data()+2, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 11 - different sizes
        EXPECT_NO_THROW(
            c = Convolve::convolve(a2, b2,
                                   Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 6);
        ippsNormDiff_Inf_64f(c.data(), r2.data()+1, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 12 - interchange 11
        EXPECT_NO_THROW(
            c = Convolve::convolve(b2, a2,
                                   Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 6);
        ippsNormDiff_Inf_64f(c.data(), r2.data()+1, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        //--------------------------------------------------------------------//
        //                   Do some extra tests on SAME and VALID            //
        //--------------------------------------------------------------------//
        std::vector<double> aSame4({1, 1, 1, 1});
        std::vector<double> bSame2({1, 1});
        std::vector<double> cSameRef1({1, 2, 2, 2});
        // Test 13 - Same, Even but unequal sizes
        EXPECT_NO_THROW(
            c = Convolve::convolve(aSame4, bSame2,
                                   Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 4);
        ippsNormDiff_Inf_64f(c.data(), cSameRef1.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 14 - Same, Even but unequal sizes reversed
        EXPECT_NO_THROW(
            c = Convolve::convolve(bSame2, aSame4,
                                   Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 4);
        ippsNormDiff_Inf_64f(c.data(), cSameRef1.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);

        // Test 15 - Same, Even and equal size
        std::vector<double> bSame4({1, 1, 1, 1});
        std::vector<double> cSameRef2({2, 3, 4, 3});
        EXPECT_NO_THROW(
            c = Convolve::convolve(aSame4, bSame4,
                                   Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 4);
        ippsNormDiff_Inf_64f(c.data(), cSameRef2.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 16 - Valid - Even and unequal sizes
        EXPECT_NO_THROW(
            c = Convolve::convolve(bSame2, aSame4,
                                   Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 3);
        for (int i=0; i<3; ++i){EXPECT_NEAR(c[i], 2.0, 1.e-10);}
        // Test 17 - Valid - Equal sizes 
        EXPECT_NO_THROW(
            c = Convolve::convolve(bSame4, aSame4,
                                   Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 1);
        EXPECT_NEAR(c[0], 4, 1.e-10);
    }
}
//============================================================================//
//int convolve_correlate_test(void)
TEST(UtilitiesConvolve, correlate)
{
    for (auto imp=0; imp<2; imp++)
    {   
        Convolve::Implementation implementation;
        if (imp == 0)
        {
            implementation = Convolve::Implementation::DIRECT;
        }
        else
        {
            implementation = Convolve::Implementation::FFT;
        }
        double emax;
        std::vector<double> a1({1, 3, 2});
        std::vector<double> b1({-1, 1, 0.5});
        std::vector<double> r1({0.5,  2.5,  3 , -1 , -2});
        std::vector<double> c;
        // Test 1
        EXPECT_NO_THROW(
            c = Convolve::correlate(a1, b1,
                                    Convolve::Mode::FULL, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 5);
        ippsNormDiff_Inf_64f(c.data(), r1.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 2 - interchange 1
        EXPECT_NO_THROW(
            c = Convolve::correlate(b1, a1,
                                    Convolve::Mode::FULL, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 5);
        emax = 0;
        for (size_t i=0; i<c.size(); i++)
        {
            emax = std::max(emax, std::abs(r1[4-i] - c[i]));
            if (std::abs(r1[4-i] - c[i]) > 1.e-10)
            {
                fprintf(stderr, "Failed corr %lf %lf on iter %d",
                        r1[4-i], c[i], imp);
            }
        }
        EXPECT_LE(emax, 1.e-10);
        // Test 3
        EXPECT_NO_THROW(
            c = Convolve::correlate(a1, b1,
                                    Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 1);
        EXPECT_NEAR(r1[2], c[0], 1.e-10);
        // Test 4 - interchange 3
        EXPECT_NO_THROW(
            c = Convolve::correlate(b1, a1,
                                    Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 1);
        EXPECT_NEAR(r1[2], c[0], 1.e-10);
        // Test 5
        EXPECT_NO_THROW(
            c = Convolve::correlate(a1, b1,
                                    Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 3);
        ippsNormDiff_Inf_64f(c.data(), r1.data()+1, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 6 - interchange 5
        EXPECT_NO_THROW(
            c = Convolve::correlate(b1, a1,
                                    Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 3);
        emax = 0;
        for (size_t i=0; i<c.size(); i++)
        {
            emax = std::max(emax, std::abs(r1[3-i] - c[i]));
            if (std::abs(r1[3-i] - c[i]) > 1.e-10)
            {
                fprintf(stderr, "Failed conv %lf %lf on iter %d",
                        r1[3-i], c[i], imp);
            }
        }
        EXPECT_LE(emax, 1.e-10);
        //--------------------------------------------------------------------//
        //                            Unequal Lengths                         //
        //--------------------------------------------------------------------//
        std::vector<double> a2({1, 3, 2, -0.2, 0.1, 4.0});
        std::vector<double> b2({-2, 2, 2.5});
        std::vector<double> r2({2.5, 9.5, 9, -2.5, -4.15, 10.6, 7.8, -8});
        // Test 7
        EXPECT_NO_THROW(
            c = Convolve::correlate(a2, b2,
                                    Convolve::Mode::FULL, implementation));
        EXPECT_EQ(c.size(), r2.size());
        ippsNormDiff_Inf_64f(c.data(), r2.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 8 - interchange 7
        EXPECT_NO_THROW(
            c = Convolve::correlate(b2, a2,
                                    Convolve::Mode::FULL, implementation));
        EXPECT_EQ(c.size(), r2.size());
        emax = 0;
        for (size_t i=0; i<c.size(); i++)
        {
            emax = std::max(emax, std::abs(r2[c.size()-1-i] - c[i]));
            if (std::abs(r2[c.size()-1-i] - c[i]) > 1.e-10)
            {
                fprintf(stderr, "Failed corr %lf %lf on iter %d",
                        r2[c.size()-1-i], c[i], imp);
            }
        }
        EXPECT_LE(emax, 1.e-10);
        // Test 9
        EXPECT_NO_THROW(
            c = Convolve::correlate(a2, b2,
                                    Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 4);
        ippsNormDiff_Inf_64f(c.data(), r2.data()+2, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 10 - interchange 9 
        EXPECT_NO_THROW(
            c = Convolve::correlate(b2, a2,
                                    Convolve::Mode::VALID, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 4);
        emax = 0;
        for (size_t i=0; i<c.size(); i++)
        {
            emax = std::max(emax, std::abs(r2[5-i] - c[i]));
            if (std::abs(r2[5-i] - c[i]) > 1.e-10)
            {
                fprintf(stderr, "Failed corr %lf %lf on iter %d",
                        r2[5-i], c[i], imp);
            }
        }
        EXPECT_LE(emax, 1.e-10);
        // Test 11
        EXPECT_NO_THROW(
            c = Convolve::correlate(a2, b2,
                                    Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 6);
        ippsNormDiff_Inf_64f(c.data(), r2.data()+1, c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 12 - interchange 11
        EXPECT_NO_THROW(
            c = Convolve::correlate(b2, a2,
                                    Convolve::Mode::SAME, implementation));
        EXPECT_EQ(static_cast<int> (c.size()), 6);
        emax = 0;
        for (size_t i=0; i<c.size(); i++)
        {
            emax = std::max(emax, std::abs(r2[6-i] - c[i]));
            if (std::abs(r2[6-i] - c[i]) > 1.e-10)
            {
                fprintf(stderr, "Failed corr %lf %lf on iter %d",
                        r2[6-i], c[i], imp);
            }
        }
        EXPECT_LE(emax, 1.e-10);
    }
}

//int convolve_autocorrelate_test(void)
TEST(UtilitiesConvolve, autoCorrelate)
{
    for (auto imp=0; imp<2; imp++)
    {
        Convolve::Implementation implementation;
        if (imp == 0)
        {
            implementation = Convolve::Implementation::DIRECT;
        }
        else
        {
            implementation = Convolve::Implementation::FFT;
        }
        std::vector<double> a1({1, 3, 2, 42, 31, -32, 34, -42, 3}); 
        std::vector<double> b1({1, 3, 2, 42, 31, -32, 34, -42, 3});
        std::vector<double> c, cref;
        // Test 1
        EXPECT_NO_THROW(
            cref = Convolve::correlate(a1, b1,
                                       Convolve::Mode::FULL, implementation));
        EXPECT_NO_THROW(
            c = Convolve::autocorrelate(a1,
                                        Convolve::Mode::FULL, implementation));
        EXPECT_EQ(cref.size(), c.size());
        double emax;
        ippsNormDiff_Inf_64f(c.data(), cref.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 2
        EXPECT_NO_THROW(
            cref = Convolve::correlate(a1, b1,
                                       Convolve::Mode::VALID, implementation));
        EXPECT_NO_THROW(
            c = Convolve::autocorrelate(a1,
                                        Convolve::Mode::VALID, implementation));
        EXPECT_EQ(cref.size(), c.size());
        ippsNormDiff_Inf_64f(c.data(), cref.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 3
        EXPECT_NO_THROW(
            cref = Convolve::correlate(a1, b1,
                                       Convolve::Mode::SAME, implementation));
        EXPECT_NO_THROW(
            c = Convolve::autocorrelate(a1,
                                        Convolve::Mode::SAME, implementation));
        EXPECT_EQ(cref.size(), c.size());
        ippsNormDiff_Inf_64f(c.data(), cref.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        //--------------------------------------------------------------------//
        //                         Different Length                           //
        //--------------------------------------------------------------------// 
        std::vector<double> a2({1, 3, 2, 42, 31, -32, 34, -42}); 
        std::vector<double> b2({1, 3, 2, 42, 31, -32, 34, -42});
        // Test 4
        EXPECT_NO_THROW(
            cref = Convolve::correlate(a2, b2,
                                       Convolve::Mode::FULL, implementation));
        EXPECT_NO_THROW(
            c = Convolve::autocorrelate(a2,
                                        Convolve::Mode::FULL, implementation));
        EXPECT_EQ(cref.size(), c.size());
        ippsNormDiff_Inf_64f(c.data(), cref.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 5
        EXPECT_NO_THROW(
            cref = Convolve::correlate(a2, b2,
                                       Convolve::Mode::VALID, implementation));
        EXPECT_NO_THROW(
            c = Convolve::autocorrelate(a2,
                                        Convolve::Mode::VALID, implementation));
        EXPECT_EQ(cref.size(), c.size());
        ippsNormDiff_Inf_64f(c.data(), cref.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
        // Test 6
        EXPECT_NO_THROW(
            cref = Convolve::correlate(a2, b2,
                                       Convolve::Mode::SAME, implementation));
        EXPECT_NO_THROW(
            c = Convolve::autocorrelate(a2,
                                        Convolve::Mode::SAME, implementation));
        EXPECT_EQ(cref.size(), c.size());
        ippsNormDiff_Inf_64f(c.data(), cref.data(), c.size(), &emax);
        EXPECT_LE(emax, 1.e-10);
    }
} 

}
