cmake_minimum_required(VERSION 3.15)

get_filename_component(PRJ_NAME ${CMAKE_SOURCE_DIR} NAME)

project(${PRJ_NAME}
    VERSION 0.1.0
    DESCRIPTION "A lightweight finite volume method CFD solver"
    LANGUAGES CXX
    HOMEPAGE_URL "https://github.com/yxcGi/CFD_FVM_Solver"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -O0 -Wall -Wextra -pedantic)
else()
    add_compile_options(-O3 -DNDEBUG)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

add_subdirectory(src)
add_subdirectory(example)
