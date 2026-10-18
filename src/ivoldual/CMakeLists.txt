CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(IVOLDUAL
        VERSION 0.3.1
        LANGUAGES CXX)

#---------------------------------------------------------

IF (NOT DEFINED IVOL_DIR)
  GET_FILENAME_COMPONENT(IVOL_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IVOL_DIR ${IVOL_ABSOLUTE_PATH} CACHE PATH "IVOL directory")
ENDIF (NOT DEFINED IVOL_DIR)

SET(SRC_NAME "ivoldual")
SET(SRC_DIR "src/${SRC_NAME}")
SET(DOXYGEN_CONFIG_FILE "${SRC_NAME}_doxygen.config")
SET(DOXYGEN_CONFIG_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${DOXYGEN_CONFIG_FILE}")
SET(TAR_FILE "${SRC_NAME}.tar")
SET(ZIP_FILE "${SRC_NAME}.zip")
SET(LIB_DIR_LIST "${CMAKE_CURRENT_BINARY_DIR}/../lib"
                 "${CMAKE_CURRENT_BINARY_DIR}/../../lib"
                 "${IJK_DIR}/lib"
                 "$ENV{HOME}/lib"
                 "/opt/lib")

#---------------------------------------------------------


IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

# Possibly turn on usage of zlib compression (requires linking with libz)
# (i.e., programs compiled with ITKNrrdIO must also be compiled with zlib)
ADD_DEFINITIONS(-DTEEM_ZLIB=1)

# Turn on TEEM_BUILD so that the proper dll export def's are 
# used on windows builds.
ADD_DEFINITIONS(-DTEEM_BUILD=1)

# *** DISABLED FOR NOW ***
#ADD_COMPILE_OPTIONS(-Wall -Wno-unused-function)

LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(NrrdIO z)

SET(CMAKE_CXX_FLAGS "-std=c++20")

# *** DISABLED FOR NOW ***
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

ADD_EXECUTABLE(ivoldual ivoldual_main.cpp ivoldual.cpp ivoldualIO.cpp
                        ivoldual_datastruct.cpp ivoldual_triangulate.cpp
                        ivoldualtable.cpp
                        ivoldual_compute.cpp ivoldual_query.cpp
                        ivoldual_move.cpp ivoldual_reposition.cpp
			ivoldual_divide_hex.cpp
                        ../ijkdual/isodual.cxx
                        ../ijkdual/ijkdual_datastruct.cxx)

TARGET_INCLUDE_DIRECTORIES(ivoldual PRIVATE
                                    "${IVOL_DIR}/include"
                                    "../ijkdual"
                                    "../ijkdualtable")

#---------------------------------------------------------

FIND_LIBRARY(NRRDIO_LIB NAMES NrrdIO HINTS ${LIB_DIR_LIST})

IF (NOT NRRDIO_LIB)
   MESSAGE(FATAL_ERROR "Unable to find library NrrdIO.\n"
                       "Source code and CMakeLists.txt to create library NrrdIO is in ${IJK_DIR}/src/NrrdIO.\n"
                       "Library NrrdIO should be in one of the following directories:\n"
                       "${LIB_DIR_LIST}")
ELSE()
   MESSAGE(STATUS "Found NrrdIO library: ${NRRDIO_LIB}.")
   TARGET_LINK_DIRECTORIES(ivoldual PRIVATE ${LIB_DIR_LIST})
   TARGET_LINK_LIBRARIES(ivoldual NrrdIO z)
ENDIF()

#---------------------------------------------------------


ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY . COMMAND tar cvfh ivoldual.tar *.cxx *.h *.txx CMakeLists.txt ivoldual_doxygen.config)

ADD_CUSTOM_TARGET(doc COMMAND doxygen ivoldual_doxygen.config)
