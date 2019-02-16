#ifndef TEST_MODULES_HPP 
#define TEST_MODULES_HPP 1
#include <string>

int readTextFile(int *npts, double *xPtr[], const std::string fileName);

int rtseis_test_modules_demean(void);
int rtseis_test_modules_detrend(void);
int rtseis_test_modules_oneBitNormalization(void);
int rtseis_test_modules_classicSTALTA(const int npts, const double x[],
                                      const std::string fileName);

#endif
