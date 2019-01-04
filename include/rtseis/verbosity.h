#ifndef RTSEIS_VERBOSITY_H__
#define RTSEIS_VERBOSITY_H__ 1
#include <stdbool.h>
#include "rtseis/config.h"
#include "rtseis/enums.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Sets the verbosity level. */
int rtseis_utils_verbosity_setLoggingLevel(const enum rtseisVerbosityLevel_enum verbosity);
/* Gets the verbosity level. */
enum rtseisVerbosityLevel_enum rtseis_utils_verbosity_getLoggingLevel(void);

bool rtseis_utils_verbosity_printError(void);
bool rtseis_utils_verbosity_printErrorAndWarning(void);
bool rtseis_utils_verbosity_printErrorAndWarningAndInfo(void);
bool rtseis_utils_verbosity_printAll(void);

#ifdef __cplusplus
}
#endif

#endif
