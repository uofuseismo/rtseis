CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)
SET(headers config.h enums.h ippsHelper.h log.h modules.h verbosity.h
    modules/demean.h modules/detrend.h)
FILE(RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
INSTALL(FILES ${headers} DESTINATION ${_relPath} COMPONENT headers)
