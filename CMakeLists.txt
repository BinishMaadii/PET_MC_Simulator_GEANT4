cmake_minimum_required(VERSION 3.16)
project(PET_MC_Simulator CXX)

# ============================================================================
# C++ Standard
# ============================================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()

# ============================================================================
# Find GEANT4 (Required for realistic physics)
# ============================================================================

find_package(Geant4 REQUIRED)

if(Geant4_FOUND)
    message(STATUS "✓ GEANT4 Found: ${Geant4_VERSION}")
    message(STATUS "  Location: ${Geant4_DIR}")
    include(${Geant4_USE_FILE})
else()
    message(FATAL_ERROR "GEANT4 not found! Install with:")
    message(FATAL_ERROR "  cd /Users/binishbatool/projects/GEANT4/")
    message(FATAL_ERROR "  git clone https://github.com/Geant4/geant4.git")
    message(FATAL_ERROR "  cd geant4 && git checkout geant4-11.2")
    message(FATAL_ERROR "  mkdir build && cd build")
    message(FATAL_ERROR "  cmake -DCMAKE_BUILD_TYPE=Release \\")
    message(FATAL_ERROR "        -DGEANT4_INSTALL_DATA=ON \\")
    message(FATAL_ERROR "        -DCMAKE_INSTALL_PREFIX=/Users/binishbatool/projects/GEANT4/geant4-install ..")
    message(FATAL_ERROR "  make -j4 && make install")
endif()

# ============================================================================
# Find ROOT (Required for visualization)
# ============================================================================

find_package(ROOT REQUIRED COMPONENTS Hist Graf Core)

if(ROOT_FOUND)
    message(STATUS "✓ ROOT Found: ${ROOT_VERSION}")
    message(STATUS "  Location: ${ROOT_DIR}")
    include(${ROOT_USE_FILE})
else()
    message(STATUS "ROOT not found. Install with: brew install root")
    message(STATUS "Or build from source:")
    message(STATUS "  cd ~/Developer/physics-libraries")
    message(STATUS "  git clone https://github.com/root-project/root.git")
    message(STATUS "  cd root && mkdir build && cd build")
    message(STATUS "  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=~/Developer/physics-libraries/root-install ..")
    message(STATUS "  make -j4 && make install")
endif()

# ============================================================================
# Include Directories
# ============================================================================

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${Geant4_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS})

# ============================================================================
# Source Files
# ============================================================================

set(SOURCES
    src/main.cpp
    src/stubs.cpp
    src/Physics/StoppingPower.cpp
    src/Transport/Particle.cpp
)

# ============================================================================
# Create Executable
# ============================================================================

add_executable(PET_Simulator ${SOURCES})

# ============================================================================
# Link Libraries
# ============================================================================

target_link_libraries(PET_Simulator PRIVATE
    ${Geant4_LIBRARIES}
    ${ROOT_LIBRARIES}
    Threads::Threads
    m
)

# ============================================================================
# Set Output Directory
# ============================================================================

set_target_properties(PET_Simulator PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# ============================================================================
# Status Messages
# ============================================================================

message(STATUS "")
message(STATUS "════════════════════════════════════════════════════")
message(STATUS "  PET MC SIMULATOR - BUILD CONFIGURATION")
message(STATUS "════════════════════════════════════════════════════")
message(STATUS "System:           ${CMAKE_SYSTEM_NAME}")
message(STATUS "C++ Compiler:     ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard:     ${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type:       ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "Physics Engine:   GEANT4 ${Geant4_VERSION}")
message(STATUS "Visualization:    ROOT ${ROOT_VERSION}")
message(STATUS "")
message(STATUS "Output Binary:    ${CMAKE_BINARY_DIR}/bin/PET_Simulator")
message(STATUS "ROOT Output:      PET_Simulation_Results.root")
message(STATUS "")
message(STATUS "Build Command:    cmake --build ${CMAKE_BINARY_DIR} -j4")
message(STATUS "Run Command:      ${CMAKE_BINARY_DIR}/bin/PET_Simulator")
message(STATUS "Visualize:        root PET_Simulation_Results.root")
message(STATUS "════════════════════════════════════════════════════")
message(STATUS "")
