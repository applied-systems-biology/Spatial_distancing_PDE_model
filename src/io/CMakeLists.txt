add_library(io SHARED
        output_handler.cpp
)

add_library(abm::io ALIAS io)
target_include_directories(io PRIVATE ${PROJECT_SOURCE_DIR}/src)