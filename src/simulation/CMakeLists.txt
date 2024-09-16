add_library(simulation SHARED
        Agent.cpp
        AgentManager.cpp
        AgentProperties.cpp
        Algorithms.cpp
        Cell.cpp
        CellFactory.cpp
        Morphology.cpp
        Simulator.cpp
        Site.cpp
        morphology/MorphologyElement.cpp
        morphology/SphereRepresentation.cpp
        morphology/SphericalMorphology.cpp
        site/CuboidSite.cpp
        MoleculeManager.cpp)
add_library(abm::simulation ALIAS simulation)
target_include_directories(simulation PRIVATE ${PROJECT_SOURCE_DIR}/src)
