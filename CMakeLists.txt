cmake_minimum_required(VERSION 3.5)
project(RTSeis VERSION 0.1.1 LANGUAGES C CXX)
enable_testing()

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

#SET(RTSEIS_VERSION ${RTSEIS_VERSION_MAJOR}.${RTSEIS_VERSION_MINOR}.${RTSEIS_VERSION_PATCH})
MESSAGE("Configuring RTSeis version: " ${VERSION})
#set(CMAKE_BUILD_TYPE Release)
# Add modules for finding requisites
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

##########################################################################################
#                    Define the compiler standard, OpenMP, and Python                    #
##########################################################################################
include(FindOpenMP QUIET)
if (OPENMP_FOUND)
   message("OpenMP found")
   string(APPEND CMAKE_C_FLAGS   " ${OpenMP_C_FLAGS}")
   string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
endif()

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
if (COMPILER_SUPPORTS_CXX17)
   #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
   #SET(CMAKE_CXX_STANDARD 17)
   message("Using CXX 17")
else()
   CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
   if (COMPILER_SUPPORTS_CXX14)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
      #SET(CMAKE_CXX_STANDARD 14)
   else()
      CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
      if (COMPILER_SUPPORTS_CXX11)
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 
         #SET(CMAKE_CXX_STANDARD 11)
      else()
         message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
      endif()
   endif()
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message("Compling with CXX Flags: ${CMAKE_CXX_FLAGS}")

if (RTSEIS_WRAP_PYTHON)
   find_package(PythonInterp 3)
   find_package(PythonLibs 3 REQUIRED)
   message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
   message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
   message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
   find_package(pybind11 REQUIRED)
endif()

##########################################################################################
#                                       Include Directories                              #
##########################################################################################

include_directories(
   ${CMAKE_SOURCE_DIR}
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   ${IPP_INCLUDE_DIR}
   ${MKL_INCLUDE_DIR}
)

ADD_SUBDIRECTORY(include)

SET(LIBALL ${MKL_LIBRARY} ${IPP_LIBRARY})
SET(LIBALL_PYTHON rtseis ${LIBALL} ${PYTHON_LIBRARIES}) #${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

SET(PYTHON_SRC src/modules/wrap.cpp)
SET(UTILS_SRCS
    src/utils/verbosity.c
    src/utils/design/filterDesigner.cpp
    src/utils/design/response.cpp
    src/utils/design/iir.cpp
    src/utils/design/fir.cpp
    src/utils/design/analogProtype.cpp
    src/utils/design/windowFunctions.cpp
    src/utils/filterRepresentations/ba.cpp
    src/utils/filterRepresentations/fir.cpp
    src/utils/filterRepresentations/sos.cpp
    src/utils/filterRepresentations/zpk.cpp
    src/utils/filters/downsample.cpp
    src/utils/filters/firfilter.cpp
    src/utils/filters/firfilterMR.cpp
    src/utils/filters/iirFilter.cpp
    src/utils/filters/iiriirFilter.cpp
    src/utils/filters/medianFilter.cpp
    src/utils/filters/sos.cpp
    src/utils/math/convolve.cpp
    src/utils/math/interpolate.cpp
    src/utils/math/polynomial.cpp
    src/utils/math/vectorMath.cpp
    src/utils/normalization/signBit.cpp
    src/utils/transforms/dft.cpp
    src/utils/transforms/dftUtils.cpp)
#SET(IPPS_SRCS
#    src/ipps/dft.c
#    src/ipps/downsample.c 
#    src/ipps/firfilter.c
#    src/ipps/iirfilter.c
#    src/ipps/medianFilter.c)
SET(MODULES_SRCS
    src/modules/detrend.cpp
    src/modules/demean.cpp
    src/modules/classicSTALTA.cpp)
#SET(DATA_SRCS src/data/waveform.cpp)
SET(PROCESSING_SRCS 
    src/processing/singleChannel/postProcessing.cpp
    src/processing/singleChannel/taper.cpp
    src/processing/singleChannel/filter.cpp)
SET(SRCS ${DATA_SRCS} ${IPPS_SRCS} ${UTILS_SRCS} ${MODULES_SRCS} ${PROCESSING_SRCS})

# cmake -DBUILD_SHARED_LIBS=YES /path/to/source
set(BUILD_SHARED_LIBS YES)
add_library(rtseis SHARED ${SRCS})
target_link_libraries(rtseis
                      PRIVATE ${MKL_LIBRARY} ${IPP_LIBRARY})
#TARGET_COMPILE_FEATURES(rtseis PUBLIC CXX_STD_14)
#SET_TARGET_PROPERTIES(rtseis PROPERTIES
#                      CXX_STANDARD_REQUIRED YES
#                      CXX_EXTENSIONS NO)
set(RTSEIS_WRAP_PYTHON YES)
IF (RTSEIS_WRAP_PYTHON)
   add_library(pyrtseis MODULE
               python/pyrtseis.cpp
               python/filterRepresentations.cpp
               python/waveform.cpp)
   target_link_libraries(pyrtseis PRIVATE pybind11::module rtseis)
   ##PYTHON_ADD_MODULE(rtseis_python ${PYTHON_SRC})
   ##TARGET_LINK_LIBRARIES(rtseis_python
   ##                      rtseis ${LIBALL})# ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
ENDIF()

#########################################################################################
#                                        Unit Tests                                     #
#########################################################################################

ADD_EXECUTABLE(testUtils
               testing/utils/utils.cpp
               testing/utils/polynomial.cpp
               testing/utils/convolve.cpp
               testing/utils/interpolate.cpp
               testing/utils/fir.cpp
               testing/utils/iir.cpp
               testing/utils/response.cpp
               testing/utils/filters.cpp
               testing/utils/transforms.cpp
               testing/utils/normalization.cpp
               testing/utils/windowFunctions.cpp)
ADD_EXECUTABLE(testPPSC
               testing/postProcessing/singleChannel.cpp)
               #testing/modules/modules.cpp
               #testing/modules/detrend.cpp
               #testing/modules/classicSTALTA.cpp)
# The core library utilities - do these first
TARGET_LINK_LIBRARIES(testUtils PRIVATE rtseis ${MKL_LIBRARY} ${IPP_LIBRARY})
# Single-channel post-processing
TARGET_LINK_LIBRARIES(testPPSC  PRIVATE rtseis ${MKL_LIBRARY} ${IPP_LIBRARY})
#SET_PROPERTY(TARGET testUtils   PROPERTY CXX_STANDARD 17)
ADD_TEST(NAME utils
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND testUtils)
ADD_TEST(NAME singleChannlePostProcessing
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND testPPSC
         DEPENDS utils)

#ADD_LIBRARY(pyrtseis SHARED src/modules/boost.cpp)
#TARGET_LINK_LIBRARIES(pyrtseis ${LIBALL} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
#SET_TARGET_PROPERTIES(pyrtseis PROPERTIES SUFFIX .so)
#SET_TARGET_PROPERTIES(pyrtseis PROPERTIES PREFIX "")

#========================================================================================#
#                                      Installation                                      #
#========================================================================================#
include(GNUInstallDirs)
if (RTSEIS_WRAP_PYTHON)
   install(TARGETS rtseis pyrtseis
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
   install(TARGETS rtseis
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

#------------------------------------------------------------------------------#
#                                   CPACK Packaging                            #
#------------------------------------------------------------------------------#
SET(CPACK_PACKAGE_NAME "rtseis")
SET(CPACK_PACKAGE_VENDOR "The Internet")
SET(CPACK_PACKAGE_CONTACT "bakerb845@gmail.com")
SET(CPACK_PACKAGE_LICENSE "MIT")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A real-time seismic signals processing library")
SET(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
#IF (WIN32 AND NOT UNIX)

#ELSE(WIN32 AND NOT UNIX)
  SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_SOURCE_GENERATOR TGZ)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "artspl-${PROJECT_VERSION}" )
#ENDIF(WIN32 AND NOT UNIX)
include(CPack)

#------------------------------------------------------------------------------#
#                              Configuration Details                           #
#------------------------------------------------------------------------------#
#SET_PROPERTY(TARGET rtseis PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
#SET(RTSeis_LIBRARY rtseis)
#SET(INCLUDE_INSTALL_DIR include/)
#SET(LIB_INSTALL_DIR lib/)
#SET(SYSCONFIG_INSTALL_DIR etc/rtseis/)
#include(CMakePackageConfigHelpers)

include(GNUInstallDirs)

#CONFIGURE_PACKAGE_CONFIG_FILE(${CMAKE_MODULE_PATH}/FindRTSeis.cmake.in
#                              ${CMAKE_BINARY_DIR}/FindRTSeis.cmake
#                              INSTALL_DESTINATION ${LIB_INSTALL_DIR}/rtseis/cmake
#                              PATH_VARS INCLUDE_INSTALL_DIR SYSCONFIG_INSTALL_DIR)
#WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/RTSeisConfigVersion.cmake
#                                 VERSION ${RTSEIS_VERSION}
#                                 COMPATIBILITY SameMajorVersion)
#install(TARGETS rtseis
#        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
#           
#        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_WRITE GROUP_READ WORLD_EXECUTE WORLD_WRITE WORLD_READ
#        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#INSTALL(FILES
#        ${CMAKE_BINARY_DIR}/FindRTSeis.cmake
#        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/RTSeisConfigVersion.cmake
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/rtseis/cmake)
