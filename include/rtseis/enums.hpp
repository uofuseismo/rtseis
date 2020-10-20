#ifndef RTSEIS_ENUMS_HPP
#define RTSEIS_ENUMS_HPP 1

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

#ifdef __cplusplus
namespace RTSeis
{
/*!
 * @brief Defines the log-level.
 */
enum class LogLevel
{
    NONE = 0,   /*!< No log information will be printed. */
    ERROR = 1,  /*!< Only errors will be printed. */
    WARN = 2,   /*!< Errors and warnings will be printed. */
    INFO = 3,   /*!< Errors, warnings, and information will be printed. */
    DEBUG = 4   /*!< Errors, warnings, information messages, and debugging
                     messages will be printed. */
};
/*!
 * @brief Defines the class's processing mode.
 */
enum class ProcessingMode
{
    POST = 0,      /*!< Indicates the module is to be used for post-processing. */
    REAL_TIME = 1, /*!< Indicates the module is to be used for real-time processing. */
    POST_PROCESSING = 0 /*!< Indicates the module is to be used for post-processing.   This is deprecated.  */
};
/*!
 * @brief Defines a classes precision.
 * @note This is deprecated and will be deleted.
 */
enum class Precision
{
    /*!< Indicates the computations underlying the module should be
        performed in floating point precision. */
    FLOAT = 1,
    /*!< Indicates the computations underlying the module should be
        performed in double precision. */
    DOUBLE = 2
};

    class Verbosity
    {
        public:
            /*!< Don't write anything. */
            static constexpr enum rtseisVerbosityLevel_enum NONE     = RTSEIS_SHOW_NONE;
            /*!< Display errors only. */
            static constexpr enum rtseisVerbosityLevel_enum ERRORS   = RTSEIS_SHOW_ERRORS;
            /*!< Display errors and warnings only. */
            static constexpr enum rtseisVerbosityLevel_enum WARNINGS = RTSEIS_SHOW_ERRORS_AND_WARNINGS;
            /*!< Display errors, warnings, and info messages. */
            static constexpr enum rtseisVerbosityLevel_enum INFO     = RTSEIS_SHOW_ERRORS_WARNINGS_AND_INFO;
            /*!< Display errors, warnings, info, and debug messages. */
            static constexpr enum rtseisVerbosityLevel_enum DEBUG    = RTSEIS_SHOW_ALL;
            /*!< Default constructor. */
            Verbosity(const enum rtseisVerbosityLevel_enum verbosity) :
                verbosity_(verbosity){};
            /*!< Sets the verbosity level. */
            void setVerbosity(const enum rtseisVerbosityLevel_enum verbosity){verbosity_ = verbosity;}
            /*!< Gets the verbosity level. */
            enum rtseisVerbosityLevel_enum getVerbosity(void) const{return verbosity_;}
        private:
            enum rtseisVerbosityLevel_enum verbosity_ = ERRORS;
    };
};
#endif

#endif
