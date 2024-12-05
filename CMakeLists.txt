# Minimum CMake version required to build this project
cmake_minimum_required(VERSION 3.10)

#Project name
project(FakeSFG)

# Add the ROOT components you need
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Ensure ROOT uses the correct flags
add_definitions(${ROOT_DEFINITIONS})

# Set the C++ standard to use
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add library target
add_library(FakeSFG SHARED src/FakeSFG.cxx)
target_link_libraries(FakeSFG PUBLIC ${ROOT_LIBRARIES})

# Add executable target
add_executable(MakeMatchingHits main.cpp src/FakeSFG.cxx)
target_link_libraries(MakeMatchingHits PUBLIC FakeSFG ${ROOT_LIBRARIES})
