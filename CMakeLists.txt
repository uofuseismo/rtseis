cmake_minimum_required(VERSION 3.5)
project(RTSeis VERSION 0.1.1 LANGUAGES C CXX)
enable_testing()

##########################################################################################
#              Define the compiler standard, OpenMP, and Python, and GTest               #
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
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
   #set(CMAKE_CXX_STANDARD 17)
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
message("Compling with CXX Flags: ${CMAKE_CXX_FLAGS}")

if (RTSEIS_WRAP_PYTHON)
   find_package(PythonInterp 3)
   find_package(PythonLibs 3 REQUIRED)
   message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
   message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
   message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
   find_package(pybind11 REQUIRED)
endif()

find_package(GTest REQUIRED)

##########################################################################################
#                                       Include Directories                              #
##########################################################################################

include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   ${IPP_INCLUDE_DIR}
   ${MKL_INCLUDE_DIR}
)
set(PRIVATE_INCLUDE_DEPENDS ${IPP_INCLUDE_DIR} ${MKL_INCLUDE_DIR})
message("Include directories:" ${PRIVATE_INCLUDE_DEPENDS})

#ADD_SUBDIRECTORY(include)

##########################################################################################
#                                        Sort Out Python Wrappers                        #
##########################################################################################
#set(RTSEIS_WRAP_PYTHON YES)
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

##########################################################################################
#                                      Define Library Dependencies                       #
##########################################################################################

configure_file(${CMAKE_SOURCE_DIR}/include/rtseis/config.h.in
               ${CMAKE_SOURCE_DIR}/include/rtseis/config.h)

set(LIBALL ${MKL_LIBRARY} ${IPP_LIBRARY})
set(LIBALL_PYTHON rtseis ${LIBALL} ${PYTHON_LIBRARIES})

#SET(PYTHON_SRC src/modules/wrap.cpp)
SET(UTILS_SRCS
    src/utilities/verbosity.c
    src/utilities/design/filterDesigner.cpp
    src/utilities/design/response.cpp
    src/utilities/design/iir.cpp
    src/utilities/design/fir.cpp
    src/utilities/design/analogProtype.cpp
    src/utilities/design/windowFunctions.cpp
    src/utilities/filterRepresentations/ba.cpp
    src/utilities/filterRepresentations/fir.cpp
    src/utilities/filterRepresentations/sos.cpp
    src/utilities/filterRepresentations/zpk.cpp
    src/utilities/filterImplementations/downsample.cpp
    src/utilities/filterImplementations/firfilter.cpp
    src/utilities/filterImplementations/multiRateFIRFilter.cpp
    src/utilities/filterImplementations/iirFilter.cpp
    src/utilities/filterImplementations/iiriirFilter.cpp
    src/utilities/filterImplementations/medianFilter.cpp
    src/utilities/filterImplementations/sos.cpp
    src/utilities/math/convolve.cpp
    src/utilities/math/interpolate.cpp
    src/utilities/math/polynomial.cpp
    src/utilities/math/vectorMath.cpp
    src/utilities/normalization/signBit.cpp
    src/utilities/transforms/dftRealToComplex.cpp
    src/utilities/transforms/dftUtils.cpp)
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
    src/postProcessing/singleChannel/waveform.cpp
    src/postProcessing/singleChannel/envelope.cpp
    src/postProcessing/singleChannel/taper.cpp
    )
SET(SRCS ${DATA_SRCS} ${IPPS_SRCS} ${UTILS_SRCS} ${MODULES_SRCS} ${PROCESSING_SRCS})

# cmake -DBUILD_SHARED_LIBS=YES /path/to/source
set(BUILD_SHARED_LIBS YES)
add_library(rtseis SHARED ${SRCS})
target_include_directories(rtseis
                           PRIVATE ${PRIVATE_INCLUDE_DEPENDS})
target_link_libraries(rtseis
                      PRIVATE ${LIBALL})
#TARGET_COMPILE_FEATURES(rtseis PUBLIC CXX_STD_14)
#SET_TARGET_PROPERTIES(rtseis PROPERTIES
#                      CXX_STANDARD_REQUIRED YES
#                      CXX_EXTENSIONS NO)
#set(RTSEIS_WRAP_PYTHON YES)

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
target_link_libraries(testUtils
                      PRIVATE rtseis ${MKL_LIBRARY} ${IPP_LIBRARY} ${GTEST_BOTH_LIBRARIES})
target_include_directories(testUtils
                           PRIVATE ${PRIVATE_INCLUDE_DEPENDS} ${GTEST_INCLUDE_DIRS})
# Single-channel post-processing
target_link_libraries(testPPSC
                      PRIVATE rtseis ${MKL_LIBRARY} ${IPP_LIBRARY} ${GTEST_BOTH_LIBRARIES})
#SET_PROPERTY(TARGET testUtils   PROPERTY CXX_STANDARD 17)
target_include_directories(testPPSC
                           PRIVATE ${PRIVATE_INCLUDE_DEPENDS} ${GTEST_INCLUDE_DIRS})

add_test(NAME utils
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND testUtils)
add_test(NAME singleChannlePostProcessing
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND testPPSC
         DEPENDS utils)

# Need to copy the answers to the test directory
file(COPY ${CMAKE_SOURCE_DIR}/testing/data DESTINATION .)
file(COPY ${CMAKE_SOURCE_DIR}/python/benchmarks.py DESTINATION .)

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

##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
#set(CPACK_PACKAGE_NAME "rtseis")
#set(CPACK_PACKAGE_VENDOR "The Internet")
#set(CPACK_PACKAGE_CONTACT "bakerb845@gmail.com")
#set(CPACK_PACKAGE_LICENSE "MIT")
#set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A real-time seismic signals processing library")
#set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
#set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
#set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
#set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
##IF (WIN32 AND NOT UNIX)
#
##ELSE(WIN32 AND NOT UNIX)
#  set(CPACK_GENERATOR "TGZ")
#  set(CPACK_SOURCE_GENERATOR TGZ)
#  set(CPACK_SOURCE_PACKAGE_FILE_NAME "rtseis-${PROJECT_VERSION}" )
##ENDIF(WIN32 AND NOT UNIX)
#include(CPack)

#------------------------------------------------------------------------------#
#                              Configuration Details                           #
#------------------------------------------------------------------------------#
#SET_PROPERTY(TARGET rtseis PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
#SET(RTSeis_LIBRARY rtseis)
#SET(INCLUDE_INSTALL_DIR include/)
#SET(LIB_INSTALL_DIR lib/)
#SET(SYSCONFIG_INSTALL_DIR etc/rtseis/)
#include(CMakePackageConfigHelpers)

#include(GNUInstallDirs)

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
