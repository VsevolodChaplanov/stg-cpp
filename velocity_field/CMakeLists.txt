set(STG_VELOCITY_FIELD_LIB stg-velocity-field)

set(STG_VELOCITY_FIELD_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_VELOCITY_FIELD_ROOT_DIR}/src/*.cpp
  ${STG_VELOCITY_FIELD_INCLUDE_DIR})

add_library(${STG_VELOCITY_FIELD_LIB} STATIC ${SOURCES} include/velocity_field/velocity_field_1d.hpp include/velocity_field/velocity_field_1d.hpp include/velocity_field/velocity_samples_1d.hpp)

target_link_libraries(${STG_VELOCITY_FIELD_LIB} PUBLIC
  ${STG_MESH_LIB}
  ${STG_TENSOR_MATRIX_LIB}
  ${STG_GENERATORS_LIB}
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)
target_link_directories(${STG_VELOCITY_FIELD_LIB} PUBLIC
  ${STG_MESH_INCLUDE_DIR}
  CONAN_PKG::boost
  CONAN_PKG::range-v3)
target_include_directories(${STG_VELOCITY_FIELD_LIB} PUBLIC
  ${STG_MESH_INCLUDE_DIR}
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_VELOCITY_FIELD_INCLUDE_DIR}
  ${STG_GENERATORS_INCLUDE_DIR}
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)

#####################
###### TESTING ######
#####################

set(VELOCITY_FIELD_TESTS stg_velocity_field_tests)

file(GLOB_RECURSE TESTSOURCES
  ${STG_VELOCITY_FIELD_ROOT_DIR}/tests/*.cpp
  ${STG_VELOCITY_FIELD_ROOT_DIR}/tests/*.hpp)

add_executable(${VELOCITY_FIELD_TESTS} ${TESTSOURCES})

target_link_libraries(${VELOCITY_FIELD_TESTS} PUBLIC
  ${STG_VELOCITY_FIELD_LIB}
  ${STG_MESH_LIB}
  CONAN_PKG::boost
  CONAN_PKG::range-v3
  CONAN_PKG::catch2)

target_link_directories(${VELOCITY_FIELD_TESTS} PUBLIC
  ${STG_VELOCITY_FIELD_INCLUDE_DIR}
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

target_include_directories(${VELOCITY_FIELD_TESTS} PUBLIC
  ${STG_VELOCITY_FIELD_INCLUDE_DIR}
  ${STG_MESH_INCLUDE_DIR}
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

catch_discover_tests(${VELOCITY_FIELD_TESTS})