cmake_minimum_required(VERSION 3.24)
include(CMakePrintHelpers)

# Project name
project(YQ_InterpolatorAndExtrapolator VERSION 1.0 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

# Define the source directory and source file
set(SRC_DIR src)
set(SRC ${SRC_DIR}/main.cpp)

# Find YAML-CPP
find_package(yaml-cpp REQUIRED)
if (yaml-cpp_FOUND)
    message(STATUS "Found yaml-cpp.")
    include_directories(${YAML_CPP_INCLUDE_DIR})
endif()


# Detect the operating system and set compiler flags accordingly
if (APPLE)
    set(CMAKE_CXX_COMPILER "clang++")
    set(YAML_CPP_LIBRARIES yaml-cpp::yaml-cpp)
    set(OpenMP_CXX_FLAGS "-Xclang -fopenmp")
    # Set OpenMP flags for Apple Clang
    if (NOT DEFINED ${OpenMP_CXX_INCLUDE_DIR})
        set(OpenMP_CXX_INCLUDE_DIR "/opt/homebrew/opt/libomp/include")
    endif()
    if (NOT DEFINED ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_CXX_LIB_NAMES "libomp")
    endif()
    if (NOT DEFINED ${OpenMP_libomp_LIBRARY})
        set(OpenMP_libomp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    endif()
    # Try to find OpenMP
    if (USE_OPENMP)
        find_package(OpenMP)
        if (OpenMP_CXX_FOUND)
            set(OPENMP_LINK OpenMP::OpenMP_CXX)
        else()
            message(STATUS "OpenMP not found. Compiling without OpenMP support.")
            set(OPENMP_LINK "")
        endif()
    else()
        message(STATUS "OpenMP disabled. Compiling without OpenMP support.")
    endif()
elseif (UNIX AND NOT APPLE)
    set(CMAKE_CXX_COMPILER "g++")
    if (USE_OPENMP)
        find_package(OpenMP)
        if (OpenMP_CXX_FOUND)
            set(OpenMP_CXX_FLAGS "-fopenmp")
        endif()
    else()
        message(STATUS "OpenMP disabled. Compiling without OpenMP support.")
    endif()
else()
    message(FATAL_ERROR "Unsupported operating system.")
endif()

# Add the executable
add_executable(${PROJECT_NAME} ${SRC})

# Link yaml-cpp
target_link_libraries(${PROJECT_NAME} ${YAML_CPP_LIBRARIES})

# Link OpenMP
target_link_libraries(${PROJECT_NAME} ${OPENMP_LINK})
if (USE_OPENMP AND OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP enabled. Compiling with OpenMP support.")
    add_definitions(-DOPENMP_AVAILABLE)
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR})