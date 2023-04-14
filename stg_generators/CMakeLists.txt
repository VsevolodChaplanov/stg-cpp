set(STG_GENERATORS_LIB stg-generators)

set(STG_GENERATORS_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_GENERATORS_ROOT_DIR}/src/*.cpp
  ${STG_GENERATORS_INCLUDE_DIR})

add_library(${STG_GENERATORS_LIB} STATIC ${SOURCES})

target_link_libraries(${STG_GENERATORS_LIB} PUBLIC
  ${STG_TENSOR_MATRIX_LIB}
  ${STG_MESH_LIB}
  ${STG_RN_GENERATOR_LIB}
  CONAN_PKG::boost
  CONAN_PKG::range-v3)
target_link_directories(${STG_GENERATORS_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::range-v3)
target_include_directories(${STG_GENERATORS_LIB} PUBLIC
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_MESH_INCLUDE_DIR}
  ${STG_RNGENERATOR_INCLUDE_DIR}
  ${STG_GENERATORS_INCLUDE_DIR}
  CONAN_PKG::boost
  CONAN_PKG::range-v3)

#####################
###### TESTING ######
#####################

set(GENERATORS_TESTS stg_generators_tests)

file(GLOB_RECURSE TESTSOURCES
  ${STG_GENERATORS_ROOT_DIR}/tests/*.cpp
  ${STG_GENERATORS_ROOT_DIR}/tests/*.hpp)

add_executable(${GENERATORS_TESTS} ${TESTSOURCES})

target_link_libraries(${GENERATORS_TESTS} PUBLIC
  ${STG_GENERATORS_LIB}
  CONAN_PKG::boost
  CONAN_PKG::range-v3
  CONAN_PKG::catch2)

target_link_directories(${GENERATORS_TESTS} PUBLIC
  ${STG_GENERATORS_INCLUDE_DIR}
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

target_include_directories(${GENERATORS_TESTS} PUBLIC
  ${STG_GENERATORS_INCLUDE_DIR}
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

catch_discover_tests(${GENERATORS_TESTS})