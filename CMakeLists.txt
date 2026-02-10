
cmake_minimum_required(VERSION 3.16)
project(OpenDRTFilmPipelineOFX LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenFX headers assumed in external/openfx
set(OFX_ROOT ${CMAKE_SOURCE_DIR}/external/openfx)

include_directories(
  ${OFX_ROOT}/include
  ${OFX_ROOT}/Support/include
  ${OFX_ROOT}/Support/Library
  ${CMAKE_SOURCE_DIR}
)

file(GLOB OFX_SUPPORT_SOURCES
  ${OFX_ROOT}/Support/Library/*.cpp
)

add_library(OpenDRTFilmPipeline MODULE
  OpenDRTFilmPipeline.cpp
  luts_embedded.cpp
  ${OFX_SUPPORT_SOURCES}
)

# OFX modules should not be prefixed with "lib"
set_target_properties(OpenDRTFilmPipeline PROPERTIES PREFIX "")

if(APPLE)
  set_target_properties(OpenDRTFilmPipeline PROPERTIES SUFFIX ".ofx")
elseif(WIN32)
  set_target_properties(OpenDRTFilmPipeline PROPERTIES SUFFIX ".ofx")
else()
  set_target_properties(OpenDRTFilmPipeline PROPERTIES SUFFIX ".ofx")
endif()

