cmake_minimum_required(VERSION 4.0)

project(VascularCircuit0D VERSION 0.1.0 LANGUAGES CXX)

#====================================================================================
# CMake configuration... Credit to Dr. Brendan Kochunas' NERS 570 SpMV project!
#====================================================================================

# Disable in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR
    "In-source builds are not allowed."
    " Create a separate directory for build files and delete CMakeCache.txt.")
endif()

# If no build type is specified, default to Debug 
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, defaulting to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# Set the C++ standard
set(SPMV_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD ${SPMV_CXX_STANDARD} CACHE STRING "Default C++ standard")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require C++ standard")

# Set module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Check that the compiler supports the C++ standard
set(MIN_GCC_VERSION 11.0) # NOTE: changed from the SpMV project... I think 11 has better C++20 support?
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${MIN_GCC_VERSION})
    message(FATAL_ERROR "GCC version must be at least ${MIN_GCC_VERSION}")
  endif()
else()
  message(WARNING 
          "You are using an unsupported compiler! "
          "Please use GCC >= ${MIN_GCC_VERSION}")
endif()
# ALSO note that I'm using Apple Clang on my Mac (I think the entire team use Mac), so this check is not comprehensive.

#===============================================================================
# Compiler flags
#===============================================================================

include(CheckCXXCompilerFlag)

# A macro to check if a compiler flag is supported, and if so, add it
# to the C++ compiler flags
macro(add_cxx_compiler_flag FLAG)
  # SFLAG is the flag without the hyphen, e.g., "-Wall" -> "Wall"
  string(REGEX REPLACE "-" "" SFLAG ${FLAG})
  # Check if the compiler supports the flag
  check_cxx_compiler_flag(${FLAG} COMPILER_SUPPORT_${SFLAG})
  # If the compiler supports the flag, add it to the C++ compiler flags
  if(COMPILER_SUPPORT_${SFLAG})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
  endif()
endmacro()

# Set warning flags
add_cxx_compiler_flag("-Wall")
add_cxx_compiler_flag("-Wextra")
add_cxx_compiler_flag("-Wcast-align")
add_cxx_compiler_flag("-Wcast-qual")
add_cxx_compiler_flag("-Wconversion")
add_cxx_compiler_flag("-Wdouble-promotion")
add_cxx_compiler_flag("-Wfloat-equal")
add_cxx_compiler_flag("-Wimplicit-fallthrough")
add_cxx_compiler_flag("-Wmissing-noreturn")
add_cxx_compiler_flag("-Wnon-virtual-dtor")
add_cxx_compiler_flag("-Wshadow")
add_cxx_compiler_flag("-Wunused")
add_cxx_compiler_flag("-Wvla")
add_cxx_compiler_flag("-Wzero-as-null-pointer-constant")
add_cxx_compiler_flag("-Wunreachable-code")
add_cxx_compiler_flag("-Wmissing-noreturn")
add_cxx_compiler_flag("-Wtautological-compare")
add_cxx_compiler_flag("-Wundef")
#add_cxx_compiler_flag("-Werror")
add_cxx_compiler_flag("-Wold-style-cast")
add_cxx_compiler_flag("-pedantic")

#===============================================================================
# Create the library
#===============================================================================

# Source files
#set(VASCULAR_CIRCUIT_SOURCES
#    "src/solver_serial.cpp"
#    )


add_library(vascular_circuit0d_lib
    src/solver_core.cpp
)
target_include_directories(vascular_circuit0d_lib 
    PUBLIC 
        include
)

add_executable(solver_serial.exe src/solver_main.cpp)

target_link_libraries(solver_serial.exe
    PRIVATE
        vascular_circuit0d_lib
)
#===============================================================================
# Testing configuration
#===============================================================================

include(CTest)

if (BUILD_TESTING)
    include(FetchContent)

    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.11.0
    )

    FetchContent_MakeAvailable(Catch2)

    add_executable(vascular_circuit_tests
        tests/test_serial_solver.cpp
    )

    target_link_libraries(vascular_circuit_tests
        PRIVATE
            vascular_circuit0d_lib
            Catch2::Catch2WithMain
    )
    include(Catch)
    catch_discover_tests(vascular_circuit_tests)
endif()

# Copy input files to build directory for easier execution
configure_file(${CMAKE_SOURCE_DIR}/input/input.inp ${CMAKE_BINARY_DIR}/input/input.inp COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/input/AAo.csv ${CMAKE_BINARY_DIR}/input/AAo.csv COPYONLY)