cmake_minimum_required(VERSION 3.5)
project(SAM)

# Set C++ standard explicitly
set(CMAKE_CXX_STANDARD 11)

# Set C++ compiler flags
if (MSVC)
    # Use MSVC-specific flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /O2 /MDd /Z7 /EHsc /openmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /DEBUG /INCREMENTAL:NO")
else()
    # Use general flags for other compilers
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g -O3 -fopenacc -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# Include directories for Eigen
include_directories(${CMAKE_SOURCE_DIR}/packages/Eigen.3.3.3/build/native/include)

# Add all .cpp files in the source directory
file(GLOB SOURCES "*.cpp")

# Specify the target executable and its dependencies
add_executable(${PROJECT_NAME} ${SOURCES})

# Manually set OpenMP flags for MSVC if not found automatically
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
else()
    if (MSVC)
        set(OpenMP_CXX_FLAGS "/openmp")
        set(OpenMP_CXX_LIB_NAMES "vcomp")
        set(OpenMP_CXX_LIBRARIES "vcomp.lib")
        set(OpenMP_EXE_LINKER_FLAGS "")
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    target_link_libraries(${PROJECT_NAME} ${OpenMP_CXX_LIBRARIES})
endif()

# Setting up OpenACC link flags if necessary
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenacc")
endif()
