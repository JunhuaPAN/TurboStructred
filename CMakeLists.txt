cmake_minimum_required (VERSION 2.6)

project(TurboStructured CXX)

set(CMAKE_CXX_COMPILER "mpicxx")
# set(CMAKE_CXX_COMPILER "g++")
SET(CMAKE_CXX_FLAGS                "-std=c++11 -w")
SET(CMAKE_CXX_FLAGS_DEBUG          "-g")
SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/libs")

#parse src/ directory for headers and source files 
file(GLOB_RECURSE src_SOURCES "src/*.cpp")
file(GLOB_RECURSE src_HEADERS "src/*.h")

set (src_INCLUDE_DIRS "")
foreach (_headerFile ${src_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND src_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES src_INCLUDE_DIRS)
include_directories(${src_INCLUDE_DIRS})

# external dependencies
include_directories ("${PROJECT_SOURCE_DIR}/libs/eigen")

# target_link_libraries(TurboStructured)

add_executable(TurboStructured "src/main.cpp")


