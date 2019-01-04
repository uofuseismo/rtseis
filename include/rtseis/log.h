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

#ifndef RTSEIS_ERRMSG 
#define RTSEIS_ERRMSG(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printError()) {\
       char errmsg[RTSEIS_MAXMSG_LEN]; \
       memset(errmsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char)); \
       sprintf(errmsg, "[ERROR]: (%s:%s:line=%d) ", __FILE__, __func__, __LINE__ );\
       do \
       {  \
           snprintf(&errmsg[strlen(errmsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       fprintf(stderr, "%s%s%s\n", ANSI_COLOR_RED, errmsg, ANSI_COLOR_RESET); \
   } \
};
#endif

#ifndef RTSEIS_WARNMSG 
#define RTSEIS_WARNMSG(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printErrorAndWarning()){ \
       char warnmsg[RTSEIS_MAXMSG_LEN]; \
       memset(warnmsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char));                           \
       sprintf(warnmsg, "[WARNING]: (%s:%s:line=%d) ", __FILE__, __func__, __LINE__ ); \
       do \
       {  \
          snprintf(&warnmsg[strlen(warnmsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       fprintf(stdout, "%s%s%s\n", ANSI_COLOR_YELLOW, warnmsg, ANSI_COLOR_RESET); \
   } \
};
#endif

#ifndef RTSEIS_INFOMSG 
#define RTSEIS_INFOMSG(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printErrorAndWarningAndInfo()) { \
       char infoMsg[RTSEIS_MAXMSG_LEN]; \
       memset(infoMsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char));                           \
       sprintf(infoMsg, "[INFO] (%s:line=%d) ", __func__, __LINE__); \
       do \
       {  \
          snprintf(&infoMsg[strlen(infoMsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       fprintf(stdout ,"%s%s%s\n", ANSI_COLOR_GREEN, infoMsg, ANSI_COLOR_RESET); \
    } \
};
#endif

#ifndef RTSEIS_DEBUGMSG 
#define RTSEIS_DEBUGMSG(fmt, ...) \
{ \
   if (rtseis_utils_verbosity_printAll()){ \
       char debugMsg[RTSEIS_MAXMSG_LEN]; \
       memset(debugMsg, 0, RTSEIS_MAXMSG_LEN*sizeof(char));                           \
       sprintf(debugMsg, "[DEBUG] (%s:line=%d) ", __func__, __LINE__); \
       do \
       {  \
          snprintf(&debugMsg[strlen(debugMsg)], RTSEIS_MAXMSG_LEN, fmt, __VA_ARGS__); \
       } while(0); \
       fprintf(stdout, "%s%s%s\n", ANSI_COLOR_BLUE, debugMsg, ANSI_COLOR_RESET); \
    } \
};
#endif

#endif
#endif
