#ifndef TEST_UTILS_HPP 
#define TEST_UTILS_HPP 1

#ifdef __cplusplus
extern "C"
{
#endif

int rtseis_test_utils_polynomial(void);
int rtseis_test_utils_convolve(void);
int rtseis_test_utils_transforms(void);
int rtseis_test_utils_design_iir_ap(void);
int rtseis_test_utils_design_zpk2sos(void);
int rtseis_test_utils_design_iir(void);
int rtseis_test_utils_design_fir_fir1(void);
int rtseis_test_utils_design_freqs(void);
int rtseis_test_utils_design_groupDelay(void);
int rtseis_test_utils_normalization(void);
int rtseis_test_utils_windowFunctions(void);

int rtseis_test_utils_filters(void);

#ifdef __cplusplus
}
#endif

#endif
