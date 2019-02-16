#ifndef RTSEIS_LOG_H__
#define RTSEIS_LOG_H__ 1
#ifdef RTSEIS_LOGGING

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "rtseis/config.h"
#include "rtseis/verbosity.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

/*
#ifndef RTSEIS_PACK_ERRMSG
#define RTSEIS_PACK_ERRMSG(fmt, ...) \
({ \
    char errMsg[RTSEIS_MAXMSG_LEN]; \
    memset(errMsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char)); \
    sprintf(errMsg, "%s[ERROR]: (%s:%s:line=%d) ", ANSI_COLOR_RED, __FILE__, __func__, __LINE__ );\
    do \
    {  \
       snprintf(&errMsg[strlen(errMsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
    } while(0); \
    sprintf(errMsg, "%s%s\n", errMsg, ANSI_COLOR_RESET); \
    errMsg; \
})
#endif
*/

#ifndef RTSEIS_ERRMSG 
#define RTSEIS_ERRMSG(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printError()) {\
       char errMsg[RTSEIS_MAXMSG_LEN]; \
       memset(errMsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char)); \
       sprintf(errMsg, "%s[ERROR]: (%s:%s:line=%d) ", ANSI_COLOR_RED, __FILE__, __func__, __LINE__ );\
       do \
       {  \
           snprintf(&errMsg[strlen(errMsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       sprintf(errMsg, "%s%s\n", errMsg, ANSI_COLOR_RESET); \
       fprintf(stderr, "%s", errMsg); \
   } \
};
#endif

#ifndef RTSEIS_THROW_IA
#define RTSEIS_THROW_IA(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printError()) {\
       char errMsg[RTSEIS_MAXMSG_LEN]; \
       memset(errMsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char)); \
       sprintf(errMsg, "%s[ERROR]: (%s:%s:line=%d) ", ANSI_COLOR_RED, __FILE__, __func__, __LINE__ );\
       do \
       {  \
           snprintf(&errMsg[strlen(errMsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       sprintf(errMsg, "%s%s\n", errMsg, ANSI_COLOR_RESET); \
       throw std::invalid_argument(errMsg); \
   } \
};
#endif

#ifndef RTSEIS_WARNMSG 
#define RTSEIS_WARNMSG(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printErrorAndWarning()){ \
       char warnMsg[RTSEIS_MAXMSG_LEN]; \
       memset(warnMsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char));                           \
       sprintf(warnMsg, "%s[WARNING]: (%s:%s:line=%d) ", ANSI_COLOR_YELLOW, __FILE__, __func__, __LINE__ ); \
       do \
       {  \
          snprintf(&warnMsg[strlen(warnMsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       sprintf(warnMsg, "%s%s\n", warnMsg, ANSI_COLOR_RESET); \
       fprintf(stdout, "%s", warnMsg); \
   } \
};
#endif

#ifndef RTSEIS_INFOMSG 
#define RTSEIS_INFOMSG(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printErrorAndWarningAndInfo()) { \
       char infoMsg[RTSEIS_MAXMSG_LEN]; \
       memset(infoMsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char));                           \
       sprintf(infoMsg, "%s[INFO] (%s:line=%d) ", ANSI_COLOR_GREEN, __func__, __LINE__); \
       do \
       {  \
          snprintf(&infoMsg[strlen(infoMsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       sprintf(infoMsg, "%s%s\n", infoMsg, ANSI_COLOR_RESET); \
       fprintf(stdout ,"%s", infoMsg); \
    } \
};
#endif

#ifndef RTSEIS_DEBUGMSG 
#define RTSEIS_DEBUGMSG(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printAll()){ \
       char debugMsg[RTSEIS_MAXMSG_LEN]; \
       memset(debugMsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char));                           \
       sprintf(debugMsg, "%s[DEBUG] (%s:line=%d) ", ANSI_COLOR_BLUE, __func__, __LINE__); \
       do \
       {  \
          snprintf(&debugMsg[strlen(debugMsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       sprintf(debugMsg, "%s%s\n", debugMsg, ANSI_COLOR_RESET); \
       fprintf(stdout, "%s", debugMsg); \
    } \
};
#endif

#endif
#endif
