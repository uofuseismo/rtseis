cmake_minimum_required(VERSION 3.16)
project(RTSeis VERSION 2.0.0 LANGUAGES CXX)
enable_testing()

# Some options
option(WRAP_PYTHON "Compile the Python bindings" OFF)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
include(CheckCXXCompilerFlag)

if (${WRAP_PYTHON})
   message("Will build Python bindings")
   find_package(pybind11 REQUIRED)
endif()
#set(FindIPP_DIR ${CMAKE_SOURCE_DIR}/cmake)
#set(FindMKL_DIR ${CMAKE_SOURCE_DIR}/cmake)
find_package(GTest REQUIRED)
#find_package(FindIPP REQUIRED)
#find_package(FindMKL REQUIRED)

if (${FindIPP_FOUND})
   message("IPP found")
   add_compile_definitions(WITH_IPP)
else()
   message("Using native filter implementations")
endif()
if (${FindMKL_FOUND})
   message("MKL found")
   add_compile_definitions(WITH_MKL)
else()
   message("Using native FFT/interpolation implementations")
endif()


##########################################################################################
#                                       Include Directories                              #
##########################################################################################

#include_directories(
#   ${CMAKE_CURRENT_SOURCE_DIR}/include
#)
#set(PRIVATE_INCLUDE_DEPENDS ${IPP_INCLUDE_DIR} ${MKL_INCLUDE_DIR})
#message("Include directories:" ${PRIVATE_INCLUDE_DEPENDS})

#ADD_SUBDIRECTORY(include)

##########################################################################################
#                                        Sort Out Python Wrappers                        #
##########################################################################################
if (WRAP_PYTHON)
   add_library(pyrtseis MODULE
               python/pyrtseis.cpp
#               python/filterRepresentations.cpp
#               python/transforms.cpp
#               python/waveform.cpp
   )
   target_link_libraries(pyrtseis PRIVATE pybind11::module rtseis)
   set_target_properties(pyrtseis PROPERTIES
                         PREFIX ""
                         CXX_STANDARD 20
                         CXX_STANDARD_REQUIRED YES
                         CXX_EXTENSIONS NO)
endif()

##########################################################################################
#                                      Define Library Dependencies                       #
##########################################################################################

configure_file(${CMAKE_SOURCE_DIR}/include/private/version.hpp.in
               ${CMAKE_SOURCE_DIR}/include/rtseis/version.hpp)
#set(PUBLIC_HEADER_DIRECTORIES
#    ${CMAKE_SOURCE_DIR}/include/rtseis)
#set(PRIVATE_HEADER_DIRECTORIES
#    ${CMAKE_SOURCE_DIR}/include/private)
set(LIBALL ${MKL_LIBRARY} ${IPP_LIBRARY})
set(LIBALL_PYTHON rtseis ${LIBALL} ${PYTHON_LIBRARIES})

#SET(PYTHON_SRC src/modules/wrap.cpp)
SET(UTILS_SRCS
    src/system.cpp
    src/vector.cpp
    src/demean.cpp
    src/window.cpp
    src/utilities/version.cpp)
SET(SRCS ${DATA_SRCS} ${IPPS_SRCS} ${UTILS_SRCS} ${MODULES_SRCS} ${PROCESSING_SRCS})

set(BUILD_SHARED_LIBS YES)
add_library(rtseis SHARED ${SRCS})
target_include_directories(rtseis
                           PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                           #PUBLIC $<BUILD_INTERFACE:${PUBLIC_HEADER_DIRECTORIES}>
                           #PUBLIC $<INSTALL_INTERFACE:${PUBLIC_HEADER_DIRECTORIES}> 
                           #PRIVATE ${PRIVATE_INCLUDE_DEPENDS})
)
target_link_libraries(rtseis
                      PRIVATE ${MKL_LIBRARY} ${IPP_LIBRARY})
set_target_properties(rtseis PROPERTIES
                      CXX_STANDARD 20
                      CXX_STANDARD_REQUIRED YES 
                      CXX_EXTENSIONS NO) 

#########################################################################################
#                                        Unit Tests                                     #
#########################################################################################

ADD_EXECUTABLE(coreTests
               testing/utils/main.cpp
               testing/utils/demean.cpp
               testing/utils/window.cpp
)
target_link_libraries(coreTests PRIVATE rtseis ${GTEST_BOTH_LIBRARIES})
target_include_directories(coreTests
                           PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           PRIVATE ${GTEST_INCLUDE_DIRS})
add_test(NAME coreTests
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing
         COMMAND coreTests)

# Need to copy the answers to the test directory
file(COPY ${CMAKE_SOURCE_DIR}/testing/data DESTINATION .)
file(COPY ${CMAKE_SOURCE_DIR}/python/benchmarks.py DESTINATION .)

#========================================================================================#
#                                      Installation                                      #
#========================================================================================#
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/${PROJECT_NAME}Config.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)
if (${WRAP_PYTHON})
   install(TARGETS rtseis pyrtseis
           EXPORT ${PROJECT_NAME}-targets
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
           INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
           COMPONENT Runtime)
else()
   install(TARGETS rtseis
           EXPORT ${PROJECT_NAME}-targets
           RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
           LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
           ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
           INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
           PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
           COMPONENT Runtime)
endif()
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/rtseis
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
export(EXPORT ${PROJECT_NAME}-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
##########################################################################################
#                                     CPACK Packaging                                    #
##########################################################################################
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VENDOR "UUSS")
set(CPACK_PACKAGE_CONTACT "ben.baker@utah.edu")
set(CPACK_PACKAGE_LICENSE "MIT")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A real-time signals processing library geared towards seismic applications.")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set(CPACK_VERBATIM_VARIABLES YES)
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_COMPONENTS_ALL libraries headers)
if (WIN32)
    set(CPACK_GENERATOR ZIP WIX)
elseif (APPLE)
    set(CPACK_GENERATOR TGZ productbuild)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CPACK_GENERATOR TGZ RPM)
else()
    set(CPACK_GENERATOR TGZ)
endif()
set(CPACK_SOURCE_IGNORE_FILES
  /\\.git/
  \\.swp
  \\.orig
  /CMakeLists\\.txt\\.user
  /private/
)
include(CPack) # Put this last!
