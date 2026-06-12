# 自动收集 example 目录及其子目录中的所有 .cpp 算例文件
file(GLOB_RECURSE EXAMPLE_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
    # 获取相对于 example/ 的路径
    file(RELATIVE_PATH EXAMPLE_REL_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${EXAMPLE_SOURCE}
    )

    # 去掉 .cpp 后缀
    string(REGEX REPLACE "\\.cpp$" "" EXAMPLE_TARGET_NAME ${EXAMPLE_REL_PATH})

    # 将目录分隔符替换成下划线
    # simple/cavity3d_unstructured.cpp -> simple_cavity3d_unstructured
    string(REPLACE "/" "_" EXAMPLE_TARGET_NAME ${EXAMPLE_TARGET_NAME})

    # 添加 example_ 前缀，避免和主程序或库重名
    set(EXAMPLE_TARGET_NAME "${EXAMPLE_TARGET_NAME}")

    add_executable(${EXAMPLE_TARGET_NAME}
        ${EXAMPLE_SOURCE}
    )

    target_link_libraries(${EXAMPLE_TARGET_NAME}
        PRIVATE
        core
    )

    target_compile_features(${EXAMPLE_TARGET_NAME}
        PRIVATE
        cxx_std_20
    )

    set_target_properties(${EXAMPLE_TARGET_NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
    )
endforeach()
