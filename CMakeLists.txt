enable_testing()

cmake_minimum_required(VERSION 3.22.0)

project(STG VERSION 0.1.0)
set(STG_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add conan packages
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
include(${CMAKE_BINARY_DIR}/conanbuildinfo_multi.cmake)
conan_basic_setup(TARGETS)

# Enable testing with catch2
include(CTest)

# include(${CONAN_BUILD_DIRS_CATCH2_DEBUG}/Catch.cmake)
include(${CONAN_BUILD_DIRS_CATCH2}/Catch.cmake)

# Find threads library
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_compile_options(-Wall)

# -fsanitize=thread)
# add_link_options(-fsanitize=thread)
include_directories(CONAN_PKG::boost
    CONAN_PKG::fmt
    CONAN_PKG::range-v3)

# Cross includes
set(STG_TENSOR_MATRIX_INCLUDE_DIR ${STG_ROOT_DIR}/tensor/include)
set(STG_STATISTICS_INCLUDE_DIR ${STG_ROOT_DIR}/statistics/include)
set(STG_RNGENERATOR_INCLUDE_DIR ${STG_ROOT_DIR}/random_generator/include)
set(STG_MESH_INCLUDE_DIR ${STG_ROOT_DIR}/mesh/include)
set(STG_VELOCITY_FIELD_INCLUDE_DIR ${STG_ROOT_DIR}/velocity_field/include)
set(STG_GENERATORS_INCLUDE_DIR ${STG_ROOT_DIR}/stg_generators/include)
set(STG_GAUSSIAN_METHOD_INCLUDE_DIR ${STG_ROOT_DIR}/gaussian_method/include)
set(STG_APPLICATION_INCLUDE_DIR ${STG_ROOT_DIR}/app/include)
set(STG_UTILITY_INCLUDE_DIR ${STG_ROOT_DIR}/utility/include)
set(STG_SANDBOX_INCLUDE_DIR ${STG_ROOT_DIR}/sandbox/include)

add_compile_definitions(D_GLIBCXX_USE_CXX11_ABI=0)

# Subdirectories with code
add_subdirectory(random_generator)
add_subdirectory(tensor)
add_subdirectory(statistics)
add_subdirectory(mesh)
add_subdirectory(velocity_field)
add_subdirectory(stg_generators)
add_subdirectory(gaussian_method)
add_subdirectory(app)
add_subdirectory(sandbox)
add_subdirectory(utility)

# Просмотреть все пакеты конан
get_property(importTargets DIRECTORY "${CMAKE_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
message(STATUS "${importTargets}")
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests_mesh)
file(INSTALL ./test_resources DESTINATION ${CMAKE_BINARY_DIR}/bin)
