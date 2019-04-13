#ifndef RTSEIS_THROW_HPP
#define RTSIES_THROW_HPP 1
#include <cstring>
#include <string>
#include <exception>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

inline std::string makeInvalidArgumentError(
  std::string msg, const char *file, const char *function, std::size_t line)
{
    std::string red("\x1b[31m");
    std::string reset("\x1b[0m");
    std::string errMsg;
    errMsg = red + "[ERROR]: (" + file + ":" + function + std::to_string(line) + "): "
           + msg + reset + "\n";
    return errMsg;
}

#ifndef RTSEIS_THROW_IA
#define RTSEIS_THROW_IA(fmt, ...) \
{ \
    char msgData[512]; \
    memset(msgData, 0, 512*sizeof(char)); \
    sprintf(msgData, "%s[ERROR]: (%s:%s:line=%d) ", ANSI_COLOR_RED, __FILE__, __func__, __LINE__ ); \
    do \
    {  \
        snprintf(&msgData[strlen(msgData)], 512, fmt, __VA_ARGS__); \
    } while(0); \
    sprintf(msgData, "%s%s\n", msgData, ANSI_COLOR_RESET); \
    std::string throwMsg(msgData); \
    throw std::invalid_argument(throwMsg); \
}
#endif

#ifndef RTSEIS_THROW_RTE
#define RTSEIS_THROW_RTE(fmt, ...) \
{ \
    char msgData[512]; \
    memset(msgData, 0, 512*sizeof(char)); \
    sprintf(msgData, "%s[ERROR]: (%s:%s:line=%d) ", ANSI_COLOR_RED, __FILE__, __func__, __LINE__ ); \
    do \
    {  \
        snprintf(&msgData[strlen(msgData)], 512, fmt, __VA_ARGS__); \
    } while(0); \
    sprintf(msgData, "%s%s\n", msgData, ANSI_COLOR_RESET); \
    std::string throwMsg(msgData); \
    throw std::runtime_error(throwMsg); \
}
#endif
/*

#ifndef RTSEIS_THROW_IA
#define RTSEIS_THROW_IA(...) \
{ \
    std::string errMsg = makeInvalidArgumentError(__VA_ARGS__, __FILE__, __func__, __LINE__ ); \
    throw std::invalid_argument(errMsg); \
}
#endif
*/

#endif
