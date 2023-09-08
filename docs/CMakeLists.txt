cmake_minimum_required(VERSION 3.16)
project(RTSeisDocs VERSION 0.1.0 LANGUAGES CXX)

find_package(Doxygen REQUIRED)

if (${DOXYGEN_FOUND})
   # set input and output files
   set(DOXYGEN_IN  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
   set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
   # request to configure the file
   configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
   message("Doxygen build started")
   file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/rtseis_layout.xml
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
   # Note: do not put "ALL" - this builds docs together with application EVERY TIME!
   add_custom_target(docs
                     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     COMMENT "Generating API documentation with Doxygen"
                     VERBATIM)
else (${DOXYGEN_FOUND})
   message("Doxygen need to be installed to generate the doxygen documentation")
endif()
