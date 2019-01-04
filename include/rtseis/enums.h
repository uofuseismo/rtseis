#ifndef RTSEIS_ENUMS_H__
#define RTSEIS_ENUMS_H__ 1

enum rtseisPrecision_enum
{
    /*!< Unknown precision. */
    RTSEIS_UNKNOWN = 0,
    /*!< Floating precision. */
    RTSEIS_FLOAT = 1,
    /*!< Double precision. */
    RTSEIS_DOUBLE = 2
};

enum rtseisVerbosityLevel_enum
{
    /*!< Don't write anything. */
    RTSEIS_SHOW_NONE = 0,
    /*!< Display errors only. */ 
    RTSEIS_SHOW_ERRORS = 1,
    /*!< Display errors and warnings. */
    RTSEIS_SHOW_ERRORS_AND_WARNINGS  = 2,
    /*!< Display errors, warnings, and info messages. */
    RTSEIS_SHOW_ERRORS_WARNINGS_AND_INFO = 3,
    /*!< Display errors, warnings, info, and debugging messages. */
    RTSEIS_SHOW_ALL = 4
};

#endif
