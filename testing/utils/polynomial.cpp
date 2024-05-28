#include <cstdio>
#include <cstdlib>
#include <complex>
#include <string>
#include <vector>
#ifdef WITH_IPP_2024
#include <ipp.h>
#else
#include <ipps.h>
#endif
#include "rtseis/utilities/math/polynomial.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace RTSeis::Utilities::Math;

/*
int polynomial_poly_test(void);
int polynomial_roots_test(void);
int polynomial_polyval_test(void);

int rtseis_test_utils_polynomial(void)
{
    int ierr;
    ierr = polynomial_poly_test();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed poly test");
        return EXIT_FAILURE;
    }

    ierr = polynomial_roots_test();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed roots test");
        return EXIT_FAILURE;
    }

    ierr = polynomial_polyval_test();
    if (ierr != EXIT_SUCCESS)
    {
        RTSEIS_ERRMSG("%s", "Failed polyval test");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
*/

//int polynomial_polyval_test(void)
TEST(UtilitiesPolynomial, polyval)
{
    std::vector<double> x({5, 7, 9, -2});
    std::vector<double> p0({3});
    std::vector<double> p1({3, 2});
    std::vector<double> p2({3, 2, 1});
    std::vector<double> p3({3, 2, 1, -2});
    std::vector<double> p4({3, 2, 1, -2, -3});
    std::vector<double> yref0({3, 3, 3, 3});
    std::vector<double> yref1({17, 23,  29, -4});
    std::vector<double> yref2({86, 162, 262, 9});
    std::vector<double> yref3({428, 1132, 2356, -20});
    std::vector<double> yref4({2137, 7921, 21201,  37});

    std::vector<std::complex<double>> zx;
    zx.push_back(std::complex<double> (5,0) );
    zx.push_back(std::complex<double> (0,7) );
    zx.push_back(std::complex<double> (2,9) );
    zx.push_back(std::complex<double> (4,-2) );
    std::vector<std::complex<double>> zp0;
    zp0.push_back( std::complex<double> (3,0));
    std::vector<std::complex<double>> zp1;
    zp1.push_back( std::complex<double> (3,0));
    zp1.push_back( std::complex<double> (0,2));
    std::vector<std::complex<double>> zp2;
    zp2.push_back( std::complex<double> (3,0));
    zp2.push_back( std::complex<double> (0,2));
    zp2.push_back( std::complex<double> (1,0));
    std::vector<std::complex<double>> zp3;
    zp3.push_back( std::complex<double> (3,0));
    zp3.push_back( std::complex<double> (0,2));
    zp3.push_back( std::complex<double> (1,0));
    zp3.push_back( std::complex<double> (1,-2));
    std::vector<std::complex<double>> zp4;
    zp4.push_back( std::complex<double> (3,0));
    zp4.push_back( std::complex<double> (0,2));
    zp4.push_back( std::complex<double> (1,0));
    zp4.push_back( std::complex<double> (1,-2));
    zp4.push_back( std::complex<double> (-3,0));
    std::vector<std::complex<double>> zyref0;
    zyref0.push_back(std::complex<double> (3,0));
    zyref0.push_back(std::complex<double> (3,0));
    zyref0.push_back(std::complex<double> (3,0));
    zyref0.push_back(std::complex<double> (3,0));
    std::vector<std::complex<double>> zyref1;
    zyref1.push_back(std::complex<double> (15,2));
    zyref1.push_back(std::complex<double> (0,23));
    zyref1.push_back(std::complex<double> (6,29));
    zyref1.push_back(std::complex<double> (12,-4));
    std::vector<std::complex<double>> zyref2;
    zyref2.push_back(std::complex<double> (76,10));
    zyref2.push_back(std::complex<double> (-160,0));
    zyref2.push_back(std::complex<double> (-248,112));
    zyref2.push_back(std::complex<double> (41,-40));
    std::vector<std::complex<double>> zyref3;
    zyref3.push_back(std::complex<double> (381,48));
    zyref3.push_back(std::complex<double> (1,-1122));
    zyref3.push_back(std::complex<double> (-1503,-2010));
    zyref3.push_back(std::complex<double> (85,-244));
    std::vector<std::complex<double>> zyref4;
    zyref4.push_back(std::complex<double> (1902,240));
    zyref4.push_back(std::complex<double> (7851,7));
    zyref4.push_back(std::complex<double> (15081,-17547));
    zyref4.push_back(std::complex<double> (-151,-1146));

    std::vector<double> y;
    EXPECT_NO_THROW(y = Polynomial::polyval(p0, x));
    EXPECT_EQ(x.size(), y.size());
    double error;
    ippsNormDiff_Inf_64f(y.data(), yref0.data(), y.size(), &error);
    ASSERT_LE(error, 1.e-14);

    EXPECT_NO_THROW(y = Polynomial::polyval(p1, x));
    EXPECT_EQ(x.size(), y.size());
    ippsNormDiff_Inf_64f(y.data(), yref1.data(), y.size(), &error);
    ASSERT_LE(error, 1.e-14);
 
    EXPECT_NO_THROW(y = Polynomial::polyval(p2, x));
    EXPECT_EQ(x.size(), y.size());
    ippsNormDiff_Inf_64f(y.data(), yref2.data(), y.size(), &error);
    ASSERT_LE(error, 1.e-14);

    EXPECT_NO_THROW(y = Polynomial::polyval(p3, x));
    EXPECT_EQ(x.size(), y.size());
    ippsNormDiff_Inf_64f(y.data(), yref3.data(), y.size(), &error);
    ASSERT_LE(error, 1.e-14);

    EXPECT_NO_THROW(y = Polynomial::polyval(p4, x));
    EXPECT_EQ(x.size(), y.size());
    ippsNormDiff_Inf_64f(y.data(), yref4.data(), y.size(), &error);
    ASSERT_LE(error, 1.e-14);

    // Complex
    std::vector<std::complex<double>> zy;
    EXPECT_NO_THROW(zy = Polynomial::polyval(zp0, zx));
    EXPECT_EQ(zx.size(), zy.size());
    error = 0;
    for (auto i=0; i<static_cast<int> (zy.size()); i++)
    {
        error = std::max(error, std::abs(zy[i] - zyref0[i]));
        if (std::abs(zy[i] - zyref0[i]) > 1.e-14)
        {
            fprintf(stderr, "Failed polyval0\n");
        }
    }
    ASSERT_LE(error, 1.e-14);

    EXPECT_NO_THROW(zy = Polynomial::polyval(zp1, zx));
    EXPECT_EQ(zx.size(), zy.size());
    error = 0;
    for (auto i=0; i<static_cast<int> (zy.size()); i++)
    {
        error = std::max(error, std::abs(zy[i] - zyref1[i]));
        if (std::abs(zy[i] - zyref1[i]) > 1.e-14)
        {
            fprintf(stderr, "Failed polyval1\n");
        }
    }
    ASSERT_LE(error, 1.e-14);

    EXPECT_NO_THROW(zy = Polynomial::polyval(zp2, zx));
    EXPECT_EQ(zx.size(), zy.size());
    error = 0;
    for (auto i=0; i<static_cast<int> (zy.size()); i++)
    {
        error = std::max(error, std::abs(zy[i] - zyref2[i]));
        if (std::abs(zy[i] - zyref2[i]) > 1.e-14)
        {
            fprintf(stderr, "Failed polyval2\n");
        }
    }
    ASSERT_LE(error, 1.e-14);

    EXPECT_NO_THROW(zy = Polynomial::polyval(zp3, zx));
    EXPECT_EQ(zx.size(), zy.size());
    error = 0;
    for (auto i=0; i<static_cast<int> (zy.size()); i++)
    {
        error = std::max(error, std::abs(zy[i] - zyref3[i]));
        if (std::abs(zy[i] - zyref3[i]) > 1.e-14)
        {
            fprintf(stderr, "Failed polyval3\n");
        }
    }
    ASSERT_LE(error, 1.e-14);

    EXPECT_NO_THROW(zy = Polynomial::polyval(zp4, zx));
    EXPECT_EQ(zx.size(), zy.size());
    error = 0;
    for (auto i=0; i<static_cast<int> (zy.size()); i++)
    {
        error = std::max(error, std::abs(zy[i] - zyref4[i]));
        if (std::abs(zy[i] - zyref4[i]) > 1.e-14)
        {
            fprintf(stderr, "Failed polyval4\n");
        }
    }
    ASSERT_LE(error, 1.e-14);
}

//int polynomial_roots_test(void)
TEST(UtilitiesPolynomial, roots)
{
    std::vector<std::complex<double>> ref1, ref2;
    ref1.resize(3);
    ref1[0] = std::complex<double> (12.122893784632391, 0);
    ref1[1] = std::complex<double> (-5.734509942225072, 0);
    ref1[2] = std::complex<double> (-0.388383842407320, 0);

    ref2.resize(3);
    ref2[0] = std::complex<double> (3.181666466658253, 8.011804223473870);
    ref2[1] = std::complex<double> (3.181666466658253,-8.011804223473870);
    ref2[2] = std::complex<double> (-0.363332933316507, 0);

    std::vector<double> p;
    p.resize(4);
    p[0] = 1; p[1] =-6; p[2] =-72; p[3] =-27; 
    std::vector<std::complex<double>> roots;
    EXPECT_NO_THROW(roots = Polynomial::roots(p));
    EXPECT_EQ(static_cast<int> (roots.size()), 3);
    double error = 0;
    for (auto i=0; i<static_cast<int> (roots.size()); i++)
    {
        error = std::max(error, std::abs(roots[i] - ref1[i]));
        if (std::abs(roots[i] - ref1[i]) > 1.e-14)
        {
            fprintf(stderr, "Failed roots 1 test\n");
        }
    }
    ASSERT_LE(error, 1.e-14);

    p[0] = 1;
    p[1] =-6;
    p[2] = 72;
    p[3] = 27;
    EXPECT_NO_THROW(roots = Polynomial::roots(p));
    EXPECT_EQ(static_cast<int> (roots.size()), 3);
    error = 0;
    for (auto i=0; i<static_cast<int> (roots.size()); i++)
    {
        error = std::max(error, std::abs(roots[i] - ref2[i]));
        if (std::abs(roots[i] - ref2[i]) > 1.e-14)
        {
            fprintf(stderr, "Failed roots 2 test\n");
        }
    }
    ASSERT_LE(error, 1.e-14);
}

//int polynomial_poly_test(void)
TEST(UtilitiesPolynomial, poly)
{
    std::vector<std::complex<double>> croots;
    croots.resize(4);
    croots[0] = std::complex<double> (1, 0);
    croots[1] = std::complex<double> (0, 1);
    croots[2] = std::complex<double> (4.3, 0);
    croots[3] = std::complex<double> (0, -2);

    std::vector<double> roots;
    roots.resize(4);
    roots[0] = 1;
    roots[1] = 1;
    roots[2] = 4.3;
    roots[3] =-2;

    std::vector<double> ref;
    ref.resize(5);
    ref[0] = 1;
    ref[1] =-4.3;
    ref[2] =-3;
    ref[3] = 14.9;
    ref[4] =-8.6;

    std::vector<std::complex<double>> cref;
    cref.resize(5);
    cref[0] = std::complex<double> (    1,   0);
    cref[1] = std::complex<double> ( -5.3, 1.0); 
    cref[2] = std::complex<double> (  6.3,-5.3); 
    cref[3] = std::complex<double> (-10.6, 4.3); 
    cref[4] = std::complex<double> (8.6, 0);

    std::vector<double> poly;
    EXPECT_NO_THROW(poly = Polynomial::poly(roots));
    EXPECT_EQ(poly.size(), ref.size());
    double error = 0;
    ippsNormDiff_Inf_64f(poly.data(), ref.data(), poly.size(), &error);
    ASSERT_LE(error, 1.e-14);

    std::vector<std::complex<double>> cpoly;
    EXPECT_NO_THROW(cpoly = Polynomial::poly(croots));
    EXPECT_EQ(cpoly.size(), cref.size());
    error = 0;
    for (size_t i=0; i<poly.size(); i++)
    {
        error = std::max(error, std::abs(cpoly[i] - cref[i]));
        if (std::abs(cpoly[i] - cref[i]) > 1.e-14)
        {
            fprintf(stderr, "Complex poly failed (%lf,%lf) (%lf,%lf)\n",
                    std::real(cpoly[i]), std::imag(cpoly[i]),
                    std::real(cref[i]),  std::imag(cref[i]));
        }
    }
    ASSERT_LE(error, 1.e-14);
}

}
