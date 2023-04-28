set(STG_MESH_LIB stg-mesh)

set(STG_MESH_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_MESH_ROOT_DIR}/src/*.cpp
  ${STG_MESH_INCLUDE_DIR}/statistics/*.hpp)

add_library(${STG_MESH_LIB} STATIC ${SOURCES})

target_link_libraries(${STG_MESH_LIB} PUBLIC
  ${STG_TENSOR_MATRIX_LIB}
  CONAN_PKG::fmt
  CONAN_PKG::boost
  CONAN_PKG::range-v3)
target_link_directories(${STG_MESH_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)
target_include_directories(${STG_MESH_LIB} PUBLIC
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_MESH_INCLUDE_DIR}
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)

#####################
###### TESTING ######
#####################

set(MESH_TESTS stg_mesh_tests)

file(GLOB_RECURSE TESTSOURCES
  ${STG_MESH_ROOT_DIR}/tests/*.cpp
  ${STG_MESH_ROOT_DIR}/tests/*.hpp)

add_executable(${MESH_TESTS} ${TESTSOURCES})

target_link_libraries(${MESH_TESTS} PUBLIC
  ${STG_TENSOR_MATRIX_LIB}
  ${STG_MESH_LIB}
  CONAN_PKG::boost
  CONAN_PKG::range-v3
  CONAN_PKG::catch2)

target_link_directories(${MESH_TESTS} PUBLIC
  ${STG_MESH_INCLUDE_DIR}/
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

target_include_directories(${MESH_TESTS} PUBLIC
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_MESH_INCLUDE_DIR}
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

catch_discover_tests(${MESH_TESTS})