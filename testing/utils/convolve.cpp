#include <stdio.h>
#include <stdlib.h>
#include <string>
#define RTSEIS_LOGGING 1
#include "utils.hpp"
#include "rtseis/utils/convolve.hpp"
#include "rtseis/log.h"

using namespace RTSeis::Utils::Math;

int convolve_convolve_test(void);
int convolve_correlate_test(void);

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
    return EXIT_SUCCESS;
}
//============================================================================//
int convolve_convolve_test(void)
{
    int ierr;
    for (int imp=0; imp<2; imp++)
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
        ierr = Convolve::convolve(a1, b1, c,
                                  Convolve::Mode::FULL, implementation);
        if (ierr != 0 || c.size() != 5)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r1[i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r1[i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 2 - interchange 1
        ierr = Convolve::convolve(b1, a1, c,
                                  Convolve::Mode::FULL, implementation);
        if (ierr != 0 || c.size() != 5)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r1[i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r1[i], c[i], imp);
                return EXIT_FAILURE;
            }   
        }  
        // Test 3
        ierr = Convolve::convolve(a1, b1, c,
                                  Convolve::Mode::VALID, implementation);
        if (ierr != 0 || c.size() != 1)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        if (std::abs(r1[2] - c[0]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                          r1[2], c[0], imp);
            return EXIT_FAILURE;
        }
        // Test 4 - interchange 3
        ierr = Convolve::convolve(b1, a1, c,
                                  Convolve::Mode::VALID, implementation);
        if (ierr != 0 || c.size() != 1)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        if (std::abs(r1[2] - c[0]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                          r1[2], c[0], imp);
            return EXIT_FAILURE;
        }
        // Test 5
        ierr = Convolve::convolve(a1, b1, c,
                                  Convolve::Mode::SAME, implementation);
        if (ierr != 0 || c.size() != 3)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r1[1+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r1[1+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 6 - interchange 5
        ierr = Convolve::convolve(b1, a1, c,
                                  Convolve::Mode::SAME, implementation);
        if (ierr != 0 || c.size() != 3)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r1[1+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r1[1+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        //--------------------------------------------------------------------//
        //                            Unequal Lengths                         //
        //--------------------------------------------------------------------//
        std::vector<double> a2({1, 3, 2, -0.2, 0.1, 4.0});
        std::vector<double> b2({-2, 2, 2.5});
        std::vector<double> r2({-2, -4, 4.5, 11.9, 4.4, -8.3, 8.25, 10});
        // Test 7
        ierr = Convolve::convolve(a2, b2, c,
                                  Convolve::Mode::FULL, implementation);
        if (ierr != 0 || c.size() != r2.size())
        {
            RTSEIS_ERRMSG("Failed to call convolve unequal on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r2[i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r2[i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 8 - interchange 7
        ierr = Convolve::convolve(b2, a2, c,
                                  Convolve::Mode::FULL, implementation);
        if (ierr != 0 || c.size() != r2.size())
        {
            RTSEIS_ERRMSG("Failed to call convolve unequal on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r2[i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r2[i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 9
        ierr = Convolve::convolve(a2, b2, c,
                                  Convolve::Mode::VALID, implementation);
        if (ierr != 0 || c.size() != 4)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::fabs(r2[2+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r2[2+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 10 - interchange 9
        ierr = Convolve::convolve(b2, a2, c,
                                  Convolve::Mode::VALID, implementation);
        if (ierr != 0 || c.size() != 4)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        if (ierr != 0 || c.size() != 4)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::fabs(r2[2+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r2[2+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 11 - different sizes
        ierr = Convolve::convolve(a2, b2, c,
                                  Convolve::Mode::SAME, implementation);
        if (ierr != 0 || c.size() != 6)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::fabs(r2[1+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r2[1+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 12 - interchange 11
        ierr = Convolve::convolve(b2, a2, c,
                                  Convolve::Mode::SAME, implementation);
        if (ierr != 0 || c.size() != 6)
        {
            RTSEIS_ERRMSG("Failed to call convolve on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::fabs(r2[1+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r2[1+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}
//============================================================================//
int convolve_correlate_test(void)
{
    int ierr;
    for (int imp=0; imp<2; imp++)
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
        std::vector<double> r1({0.5,  2.5,  3 , -1 , -2});
        std::vector<double> c;
        // Test 1
        ierr = Convolve::correlate(a1, b1, c,
                                   Convolve::Mode::FULL, implementation);
        if (ierr != 0 || c.size() != 5)
        {
            RTSEIS_ERRMSG("Failed to call correlate on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r1[i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r1[i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 2 - interchange 1
        ierr = Convolve::correlate(b1, a1, c,
                                   Convolve::Mode::FULL, implementation);
        if (ierr != 0 || c.size() != 5)
        {
            RTSEIS_ERRMSG("Failed to call correlate on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r1[4-i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r1[4-i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 3
        ierr = Convolve::correlate(a1, b1, c,
                                   Convolve::Mode::VALID, implementation);
        if (ierr != 0 || c.size() != 1)
        {
            RTSEIS_ERRMSG("Failed to call correlate on round %d", imp);
            return EXIT_FAILURE;
        }
        if (std::abs(r1[2] - c[0]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                          r1[2], c[0], imp);
            return EXIT_FAILURE;
        }
        // Test 4 - interchange 3
        ierr = Convolve::correlate(b1, a1, c,
                                   Convolve::Mode::VALID, implementation);
        if (ierr != 0 || c.size() != 1)
        {
            RTSEIS_ERRMSG("Failed to call correlate on round %d", imp);
            return EXIT_FAILURE;
        }
        if (std::abs(r1[2] - c[0]) > 1.e-10)
        {
            RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                          r1[2], c[0], imp);
            return EXIT_FAILURE;
        }
        // Test 5
        ierr = Convolve::correlate(a1, b1, c,
                                   Convolve::Mode::SAME, implementation);
        if (ierr != 0 || c.size() != 3)
        {
            RTSEIS_ERRMSG("Failed to call correlate on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r1[1+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r1[1+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 6 - interchange 5
        ierr = Convolve::correlate(b1, a1, c,
                                   Convolve::Mode::SAME, implementation);
        if (ierr != 0 || c.size() != 3)
        {
            RTSEIS_ERRMSG("Failed to call correlate on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r1[3-i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed conv %lf %lf on iter %d",
                              r1[3-i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        //--------------------------------------------------------------------//
        //                            Unequal Lengths                         //
        //--------------------------------------------------------------------//
        std::vector<double> a2({1, 3, 2, -0.2, 0.1, 4.0});
        std::vector<double> b2({-2, 2, 2.5});
        std::vector<double> r2({2.5, 9.5, 9, -2.5, -4.15, 10.6, 7.8, -8});
        // Test 7
        ierr = Convolve::correlate(a2, b2, c,
                                   Convolve::Mode::FULL, implementation);
        if (ierr != 0 || c.size() != r2.size())
        {
            RTSEIS_ERRMSG("Failed to call correlate unequal on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r2[i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r2[i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 8 - interchange 7
        ierr = Convolve::correlate(b2, a2, c,
                                   Convolve::Mode::FULL, implementation);
        if (ierr != 0 || c.size() != r2.size())
        {
            RTSEIS_ERRMSG("Failed to call correlate unequal on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r2[c.size()-1-i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r2[c.size()-1-i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 9
        ierr = Convolve::correlate(a2, b2, c,
                                   Convolve::Mode::VALID, implementation);
        if (ierr != 0 || c.size() != 4)
        {
            RTSEIS_ERRMSG("Failed to call correlate unequal on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r2[2+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r2[2+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 10 - interchange 9 
        ierr = Convolve::correlate(b2, a2, c,
                                   Convolve::Mode::VALID, implementation);
        if (ierr != 0 || c.size() != 4)
        {
            RTSEIS_ERRMSG("Failed to call correlate unequal on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r2[5-i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r2[5-i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 11
        ierr = Convolve::correlate(a2, b2, c,
                                   Convolve::Mode::SAME, implementation);
        if (ierr != 0 || c.size() != 6)
        {
            RTSEIS_ERRMSG("Failed to call correlate unequal on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r2[1+i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r2[1+i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
        // Test 12 - interchange 11
        ierr = Convolve::correlate(b2, a2, c,
                                   Convolve::Mode::SAME, implementation);
        if (ierr != 0 || c.size() != 6)
        {
            RTSEIS_ERRMSG("Failed to call correlate unequal on round %d", imp);
            return EXIT_FAILURE;
        }
        for (int i=0; i<c.size(); i++)
        {
            if (std::abs(r2[6-i] - c[i]) > 1.e-10)
            {
                RTSEIS_ERRMSG("Failed corr %lf %lf on iter %d",
                              r2[6-i], c[i], imp);
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}
