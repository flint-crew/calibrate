cmake_minimum_required(VERSION 3.10)
project(calibrate VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -march=native")

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)
set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore 3.6 REQUIRED COMPONENTS casa ms tables measures)

# Check if casacore was found
if(NOT CASACORE_FOUND)
    message(FATAL_ERROR "Casacore not found. Please install casacore or specify its location with -DCMAKE_PREFIX_PATH=/path/to/casacore")
endif()

# Find other required libraries
find_library(GSL_LIBRARY gsl)
find_library(GSL_CBLAS_LIBRARY gslcblas)
find_library(BOOST_THREAD_LIBRARY boost_thread)
find_library(BOOST_SYSTEM_LIBRARY boost_system)

if(NOT GSL_LIBRARY OR NOT GSL_CBLAS_LIBRARY)
    message(FATAL_ERROR "GSL libraries not found. Please install gsl")
endif()

if(NOT BOOST_THREAD_LIBRARY OR NOT BOOST_SYSTEM_LIBRARY)
    message(FATAL_ERROR "Boost libraries not found. Please install boost")
endif()

# Include directories - use the ones from FindCasacore
include_directories(${CASACORE_INCLUDE_DIRS})

# Common source files for calibrate and addmodel
set(COMMON_SOURCES
    predicter.cpp
    mspredicter.cpp
    calibrator.cpp
    calibrationmethod.cpp
    polynomialfitter.cpp
    model/model.cpp
    nlplfitter.cpp
    progressbar.cpp
)

# Create object library for common sources
add_library(calibrate_common OBJECT ${COMMON_SOURCES})

# Set target properties for the object library
set_target_properties(calibrate_common PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

# Main calibrate executable
add_executable(calibrate calibrate.cpp)
target_link_libraries(calibrate calibrate_common)
target_link_libraries(calibrate
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${CASACORE_LIBRARIES}
    ${BOOST_THREAD_LIBRARY}
    ${BOOST_SYSTEM_LIBRARY}
)

# Addmodel executable
add_executable(addmodel addmodel.cpp)
target_link_libraries(addmodel calibrate_common)
target_link_libraries(addmodel
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${CASACORE_LIBRARIES}
    ${BOOST_THREAD_LIBRARY}
    ${BOOST_SYSTEM_LIBRARY}
)

# Applysolutions executable (only needs CASA libraries)
add_executable(applysolutions applysolutions.cpp)
target_link_libraries(applysolutions
    ${CASACORE_LIBRARIES}
)

# Installation
install(TARGETS calibrate addmodel applysolutions
    RUNTIME DESTINATION bin
)

# Print configuration info
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Casacore include dirs: ${CASACORE_INCLUDE_DIRS}")
message(STATUS "Casacore libraries: ${CASACORE_LIBRARIES}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# Custom target for cleaning object files
add_custom_target(clean-objs
    COMMAND ${CMAKE_COMMAND} -E remove $<TARGET_OBJECTS:calibrate_common>
    COMMENT "Cleaning object files"
)
