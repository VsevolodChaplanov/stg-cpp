set(STG_RN_GENERATOR_LIB stg-rn-generator)

set(STG_RNGENERATOR_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_RNGENERATOR_ROOT_DIR}/src/*.cpp
  ${STG_RNGENERATOR_ROOT_DIR}/include/stg_random/*.hpp)

add_library(${STG_RN_GENERATOR_LIB} STATIC ${SOURCES})

target_link_libraries(${STG_RN_GENERATOR_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::armadillo
  CONAN_PKG::range-v3)
target_link_directories(${STG_RN_GENERATOR_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::armadillo
  CONAN_PKG::range-v3)

target_include_directories(${STG_RN_GENERATOR_LIB} PUBLIC
  ${STG_MESH_INCLUDE_DIR}
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_RNGENERATOR_ROOT_DIR}/include/
  ${STG_RNGENERATOR_ROOT_DIR}/include/stg_random
  CONAN_PKG::boost
  CONAN_PKG::armadillo
  CONAN_PKG::range-v3)

# ####################
# ##### TESTING ######
# ####################
set(RNGENERATOR_TEST stg_rn_generator_tests)

file(GLOB TESTSOURCES
  ${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp
  ${CMAKE_CURRENT_LIST_DIR}/tests/*.hpp)

add_executable(${RNGENERATOR_TEST} ${TESTSOURCES})

target_link_libraries(${RNGENERATOR_TEST} PUBLIC
  ${STG_RN_GENERATOR_LIB}
  CONAN_PKG::boost
  CONAN_PKG::armadillo
  CONAN_PKG::catch2)

target_link_directories(${RNGENERATOR_TEST} PUBLIC
  ${STG_RNGENERATOR_ROOT_DIR}/include/
  ${STG_RNGENERATOR_ROOT_DIR}/include/stg_random/
  CONAN_PKG::catch2
  CONAN_PKG::armadillo
  CONAN_PKG::boost)

target_include_directories(${RNGENERATOR_TEST} PUBLIC
  ${STG_MESH_INCLUDE_DIR}
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_RNGENERATOR_ROOT_DIR}/include/
  ${STG_RNGENERATOR_ROOT_DIR}/include/stg_random
  CONAN_PKG::boost
  CONAN_PKG::armadillo)

catch_discover_tests(${RNGENERATOR_TEST})