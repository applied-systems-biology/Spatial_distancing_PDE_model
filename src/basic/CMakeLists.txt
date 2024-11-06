add_library(basic SHARED
        Coordinate3D.cpp
        Randomizer.cpp
        )
add_library(abm::basic ALIAS basic)
target_include_directories(basic PRIVATE ${PROJECT_SOURCE_DIR}/src)