#include <stdio.h>
#include <stdlib.h>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
#include "rtseis/verbosity.h"

/*!
 * @defgroup rtseis_ispl_verbosity
 * @brief Manages the logging level.  This is from ISTI's ISPL.  The function
 *        names have been modified to conform with rtseis.
 * @ingroup rtseis_ispl
 * @author Ben Baker
 * @copyright ISTI distributed under the Apache 2 license.
 */
/*!
 * @brief Internal variable defining the logging level.
 * @ingroup rtseis_utils_verbosity
 */
static enum rtseisVerbosityLevel_enum verb = RTSEIS_SHOW_ERRORS;

/*!
 * @brief Sets the verbosity level.
 * @param[in] verbosity   RTSEIS_SHOW_NONE (0) will not output any message. 
 * @param[in] verbosity   RTSEIS_SHOW_ERRORS (1) is the default and will show
 *                        error messages.
 * @param[in] verbosity   RTSEIS_SHOW_ERRORS_AND_WARNINGS (2) will show errors
 *                        and warnings.
 * @param[in] verbosity   RTSEIS_SHOW_ERRORS_WARNINGS_AND_INFO (3) will show
 *                        errors, warnings, and info.
 * @param[in] verbosity   RTSEIS_SHOW_ALL (4) will show errors, warnings, and
 *                        intermediate debugging information.
 * @result 0 indicates success. 
 * @ingroup rtseis_ispl_verbosity
 */
int rtseis_utils_verbosity_setLoggingLevel(
    const enum rtseisVerbosityLevel_enum verbosity)
{
    if ((int) verbosity < 0 || (int) verbosity > 4)
    {
        RTSEIS_ERRMSG("Invalid verbosity level %d", (int) verbosity);
        return -1;
    }
    verb = verbosity;
    return 0;
}
//============================================================================//
/*!
 * @brief Gets the verbosity level. 
 * @result The verbosity level.
 * @ingroup rtseis_ispl_verbosity
 */
enum rtseisVerbosityLevel_enum rtseis_utils_verbosity_getLoggingLevel(void)
{
    return verb;
}
/*!
 * @brief Logical function to determine if errors are printed.
 * @result If true then the logging level prints warnings.
 * @ingroup rtseis_ispl_verbosity
 */
bool rtseis_utils_verbosity_printError(void)
{
    enum rtseisVerbosityLevel_enum verbosity;
    verbosity = rtseis_utils_verbosity_getLoggingLevel();
    if ((int) verbosity < 1){return false;}
    return true;
}
/*!
 * @brief Logical function to determine if errors and warnings are printed.
 * @result If true then the logging level prints warnings and errors.
 * @brief rtseis_ispl_verbosity
 */
bool rtseis_utils_verbosity_printErrorAndWarning(void)
{
    enum rtseisVerbosityLevel_enum verbosity;
    verbosity = rtseis_utils_verbosity_getLoggingLevel();
    if ((int) verbosity < 2){return false;}
    return true;
}
/*!
 * @brief Logical function to determine if errors, warnings, and general
 *        information are printed.
 * @result If true then the logging level prints errors, warnings, and
 *         general information.
 * @ingroup rtseis_ispl_verbosity
 */
bool rtseis_utils_verbosity_printErrorAndWarningAndInfo(void)
{
    enum rtseisVerbosityLevel_enum verbosity;
    verbosity = rtseis_utils_verbosity_getLoggingLevel();
    if ((int) verbosity < 3){return false;}
    return true;
}
/*!
 * @brief Logical function to determine if errors, warning, general
 *        information, and debugging information are printed.
 * @result If true then the logging level prints errors, warning,
 *         general information, and debugging information.
 * @ingroup rtseis_ispl_verbosity
 */
bool rtseis_utils_verbosity_printAll(void)
{
    enum rtseisVerbosityLevel_enum verbosity;
    verbosity = rtseis_utils_verbosity_getLoggingLevel();
    if ((int) verbosity < 4){return false;}
    return true;
}
