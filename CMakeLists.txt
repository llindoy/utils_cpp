cmake_minimum_required(VERSION 3.6) #we require FindBLAS for this to work
project(linalg LANGUAGES CXX)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)

add_library(utils INTERFACE)
target_include_directories(utils INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

set(UTILS_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "")



