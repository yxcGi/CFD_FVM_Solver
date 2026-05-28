# 所有源文件
file(GLOB_RECURSE CORE_SOURCES
    Geometry/*.cpp
    Math/*.cpp
)

# 所有头文件
file(GLOB_RECURSE CORE_HEADERS
    Geometry/*.h
    Geometry/*.hpp
    Math/*.h
    Math/*.hpp
    Fields/*.h
    Fields/*.hpp
    FVM/*.hpp
    thread_pool/*.hpp
    algorithm/*.hpp
)

# 添加动态库
add_library(core SHARED ${CORE_SOURCES})
# 添加头文件
target_include_directories(
    core
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/Geometry
    ${CMAKE_CURRENT_SOURCE_DIR}/Math
    ${CMAKE_CURRENT_SOURCE_DIR}/Fields
    ${CMAKE_CURRENT_SOURCE_DIR}/Fields/FieldOperators
    ${CMAKE_CURRENT_SOURCE_DIR}/FVM
    ${CMAKE_CURRENT_SOURCE_DIR}/thread_pool
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithm
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithm/SIMPLE
)