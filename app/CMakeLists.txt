set(STG_APPLICATION_LIB stg-application)


set(STG_APPLICATION_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_APPLICATION_ROOT_DIR}/src/*.cpp
  ${STG_APPLICATION_ROOT_DIR}/include/stg/*)

add_library(${STG_APPLICATION_LIB} STATIC ${SOURCES})

target_link_libraries(${STG_APPLICATION_LIB} PUBLIC
  ${STG_MESH_LIB}
  ${STG_STATISTICS_LIB}
  ${STG_TENSOR_MATRIX_LIB}
  ${STG_GENERATORS_LIB}
  ${STG_VELOCITY_FIELD_LIB}
  ${STG_GENERATORS_LIB}
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)
target_link_directories(${STG_APPLICATION_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::range-v3)
target_include_directories(${STG_APPLICATION_LIB} PUBLIC
  ${STG_MESH_INCLUDE_DIR}
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_VELOCITY_FIELD_INCLUDE_DIR}
  ${STG_GENERATORS_INCLUDE_DIR}
  ${STG_APPLICATION_INCLUDE_DIR}
  ${STG_GENERATORS_INCLUDE_DIR}
  ${STG_RNGENERATOR_INCLUDE_DIR}
  ${STG_STATISTICS_INCLUDE_DIR}
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)

###########################
###### SPECTRAL MAIN ######
###########################

set(SPECTRAL_METHOD spectral_method)
add_executable(${SPECTRAL_METHOD} spectral_main.cpp kriging_main.cpp)
target_link_libraries(${SPECTRAL_METHOD} PUBLIC
  ${STG_APPLICATION_LIB})
target_link_directories(${SPECTRAL_METHOD} PUBLIC
  ${STG_MESH_INCLUDE_DIR}
  CONAN_PKG::boost
  CONAN_PKG::range-v3)
target_include_directories(${SPECTRAL_METHOD} PUBLIC
  ${STG_MESH_INCLUDE_DIR}
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_VELOCITY_FIELD_INCLUDE_DIR}
  ${STG_GENERATORS_INCLUDE_DIR}
  ${STG_APPLICATION_INCLUDE_DIR})

#####################
###### TESTING ######
#####################

set(APPLICATION_TESTS stg_application_tests)

file(GLOB_RECURSE TESTSOURCES
  ${STG_APPLICATION_ROOT_DIR}/tests/*.cpp
  ${STG_APPLICATION_ROOT_DIR}/tests/*.hpp)

add_executable(${APPLICATION_TESTS} ${TESTSOURCES})

target_link_libraries(${APPLICATION_TESTS} PUBLIC
  ${STG_APPLICATION_LIB}
  ${STG_MESH_LIB}
  CONAN_PKG::boost
  CONAN_PKG::range-v3
  CONAN_PKG::catch2)

target_link_directories(${APPLICATION_TESTS} PUBLIC
  ${STG_VELOCITY_FIELD_INCLUDE_DIR}
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

target_include_directories(${APPLICATION_TESTS} PUBLIC
  ${STG_VELOCITY_FIELD_INCLUDE_DIR}
  ${STG_MESH_INCLUDE_DIR}
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

catch_discover_tests(${APPLICATION_TESTS})

add_executable(kriging kriging_main.cpp include/stg/kriging/kriging_1d_analysis.hpp)

target_link_libraries(kriging PUBLIC
        ${STG_APPLICATION_LIB}
        ${STG_MESH_LIB}
        CONAN_PKG::boost
        CONAN_PKG::range-v3
        CONAN_PKG::catch2)

target_link_directories(kriging PUBLIC
        ${STG_VELOCITY_FIELD_INCLUDE_DIR}
        CONAN_PKG::catch2
        CONAN_PKG::range-v3
        CONAN_PKG::boost)

target_include_directories(kriging PUBLIC
        ${STG_VELOCITY_FIELD_INCLUDE_DIR}
        ${STG_MESH_INCLUDE_DIR}
        CONAN_PKG::catch2
        CONAN_PKG::range-v3
        CONAN_PKG::boost)