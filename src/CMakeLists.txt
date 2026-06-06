add_subdirectory(core)

add_executable(${PRJ_NAME}
    main.cpp
)

target_link_libraries(${PRJ_NAME}
    PRIVATE
    core
)