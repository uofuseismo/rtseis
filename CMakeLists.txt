CMAKE_MINIMUM_REQUIRED(VERSION 3.2 FATAL_ERROR)

PROJECT(RTSeis)
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
ENABLE_TESTING()

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

SET(RTSEIS_VERSION_MAJOR 1)
SET(RTSEIS_VERSION_MINOR 1)
SET(RTSEIS_VERSION_PATCH 1)
SET(RTSEIS_VERSION ${RTSEIS_VERSION_MAJOR}.${RTSEIS_VERSION_MINOR}.${RTSEIS_VERSION_PATCH})
MESSAGE("Configuring RTSeis version: " ${RTSEIS_VERSION})
#set(CMAKE_BUILD_TYPE Release)
# Add modules for finding requisites
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

INCLUDE(FindOpenMP QUIET)
IF (OPENMP_FOUND)
   MESSAGE("OpenMP Found")
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
IF (COMPILER_SUPPORTS_CXX17)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
ELSE()
   CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
   IF (COMPILER_SUPPORTS_CXX14)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
   ELSE()
      CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
      IF (COMPILER_SUPPORTS_CXX11)
         SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      ELSE()
         MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
      ENDIF()
   ENDIF()
ENDIF()

IF (RTSEIS_WRAP_PYTHON)
   FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
   FIND_PACKAGE(PythonInterp 3)
   FIND_PACKAGE(PythonLibs 3 REQUIRED)
   MESSAGE(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
   MESSAGE(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
   MESSAGE(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
   MESSAGE(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
   MESSAGE(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")

   FIND_PACKAGE(pybind11 REQUIRED)
ENDIF()

INCLUDE_DIRECTORIES(
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   ${IPP_INCLUDE_DIR}
   ${MKL_INCLUDE_DIR}
   ${PYTHON_INCLUDE_DIRS}
   ${Boost_INCLUDE_DIRS} 
   ${CJSON_INCLUDE_DIR})

ADD_SUBDIRECTORY(include)

SET(LIBALL ${MKL_LIBRARY} ${IPP_LIBRARY} ${CJSON_LIBRARY})
SET(LIBALL_PYTHON rtseis_shared ${LIBALL} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

SET(PYTHON_SRC src/modules/wrap.cpp)
SET(UTILS_SRCS
    src/utils/verbosity.c
    src/utils/ipps.cpp
    src/utils/design/response.cpp
    src/utils/design/iir.cpp
    src/utils/design/fir.cpp
    src/utils/design/ba.cpp
    src/utils/design/sos.cpp
    src/utils/design/zpk.cpp
    src/utils/design/analogProtype.cpp
    src/utils/filters/downsample.cpp
    src/utils/filters/firfilter.cpp
    src/utils/filters/medianFilter.cpp
    src/utils/filters/sos.cpp
    src/utils/math/convolve.cpp
    src/utils/math/polynomial.cpp)
SET(IPPS_SRCS
    src/ipps/dft.c
    src/ipps/downsample.c 
    src/ipps/firfilter.c
    src/ipps/iirfilter.c
    src/ipps/medianFilter.c)
SET(MODULES_SRCS
    src/modules/detrend.cpp
    src/modules/demean.cpp
    src/modules/fir.cpp
    src/modules/iir.cpp)
SET(DATA_SRCS src/data/waveform.cpp)
SET(SRCS ${DATA_SRCS} ${IPPS_SRCS} ${UTILS_SRCS} ${MODULES_SRCS})

ADD_LIBRARY(rtseis_shared SHARED ${SRCS})
TARGET_LINK_LIBRARIES(rtseis_shared ${LIBALL})
IF (RTSEIS_WRAP_PYTHON)
   ADD_LIBRARY(pyrtseis MODULE python/pyrtseis.cpp)
   TARGET_LINK_LIBRARIES(pyrtseis PRIVATE pybind11::module)
   #PYTHON_ADD_MODULE(rtseis_python ${PYTHON_SRC})
   #TARGET_LINK_LIBRARIES(rtseis_python
   #                      rtseis_shared ${LIBALL})# ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
ENDIF()

#########################################################################################
#                                        Unit Tests                                     #
#########################################################################################

ADD_EXECUTABLE(testUtils
               testing/utils/utils.cpp
               testing/utils/polynomial.cpp
               testing/utils/convolve.cpp
               testing/utils/fir.cpp
               testing/utils/iir.cpp
               testing/utils/response.cpp
               testing/utils/filters.cpp)
ADD_EXECUTABLE(testModules
               testing/modules/modules.cpp
               testing/modules/detrend.cpp)

TARGET_LINK_LIBRARIES(testUtils   rtseis_shared ${LIBALL})
TARGET_LINK_LIBRARIES(testModules rtseis_shared ${LIBALL})
ADD_TEST(NAME utils
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND testUtils)
ADD_TEST(NAME modules
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND testModules)

#ADD_LIBRARY(pyrtseis SHARED src/modules/boost.cpp)
#TARGET_LINK_LIBRARIES(pyrtseis ${LIBALL} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
#SET_TARGET_PROPERTIES(pyrtseis PROPERTIES SUFFIX .so)
#SET_TARGET_PROPERTIES(pyrtseis PROPERTIES PREFIX "")

#------------------------------------------------------------------------------#
#                                   CPACK Packaging                            #
#------------------------------------------------------------------------------#
SET(CPACK_PACKAGE_NAME "rtseis")
SET(CPACK_PACKAGE_VENDOR "The Internet")
SET(CPACK_PACKAGE_CONTACT "bakerb845@gmail.com")
SET(CPACK_PACKAGE_LICENSE "MIT")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A real-time seismic signals processing library")
SET(CPACK_PACKAGE_VERSION_MAJOR ${RTSEIS_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${RTSEIS_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${RTSEIS_VERSION_PATCH})
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
#IF (WIN32 AND NOT UNIX)

#ELSE(WIN32 AND NOT UNIX)
  SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_SOURCE_GENERATOR TGZ)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "artspl-${RTSEIS_VERSION}" )
#ENDIF(WIN32 AND NOT UNIX)
INCLUDE(CPack)

#------------------------------------------------------------------------------#
#                              Configuration Details                           #
#------------------------------------------------------------------------------#
SET_PROPERTY(TARGET rtseis_shared PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(RTSeis_LIBRARY sspl)
SET(INCLUDE_INSTALL_DIR include/)
SET(LIB_INSTALL_DIR lib/)
SET(SYSCONFIG_INSTALL_DIR etc/rtseis/)
INCLUDE(CMakePackageConfigHelpers)
#CONFIGURE_PACKAGE_CONFIG_FILE(${CMAKE_MODULE_PATH}/FindRTSeis.cmake.in
#                              ${CMAKE_BINARY_DIR}/FindRTSeis.cmake
#                              INSTALL_DESTINATION ${LIB_INSTALL_DIR}/rtseis/cmake
#                              PATH_VARS INCLUDE_INSTALL_DIR SYSCONFIG_INSTALL_DIR)
#WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/RTSeisConfigVersion.cmake
#                                 VERSION ${RTSEIS_VERSION}
#                                 COMPATIBILITY SameMajorVersion)
#INSTALL(TARGETS ${TARGET_NAMES}
#        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_WRITE GROUP_READ WORLD_EXECUTE WORLD_WRITE WORLD_READ
#        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#INSTALL(FILES
#        ${CMAKE_BINARY_DIR}/FindRTSeis.cmake
#        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/RTSeisConfigVersion.cmake
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/rtseis/cmake)
