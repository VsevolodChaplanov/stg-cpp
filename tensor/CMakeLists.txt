set(STG_TENSOR_MATRIX_LIB stg-tensor)

set(STG_TENSOR_MATRIX_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_TENSOR_MATRIX_ROOT_DIR}/src/*.cpp
  ${STG_TENSOR_MATRIX_INCLUDE_DUR}/stg_tensor/*.hpp)

add_library(${STG_TENSOR_MATRIX_LIB} STATIC ${SOURCES})

target_link_libraries(${STG_TENSOR_MATRIX_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::range-v3)
target_link_directories(${STG_TENSOR_MATRIX_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::range-v3)
target_include_directories(${STG_TENSOR_MATRIX_LIB} PUBLIC
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}/stg_tensor/
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}/
  CONAN_PKG::boost
  CONAN_PKG::range-v3)

# ####################
# ##### TESTING ######
# ####################
set(TENSOR_TESTS stg_tensor_tests)

file(GLOB TESTSOURCES
  ${STG_TENSOR_MATRIX_ROOT_DIR}/tests/*.cpp
  ${STG_TENSOR_MATRIX_ROOT_DIR}/tests/*.hpp)

add_executable(${TENSOR_TESTS} ${TESTSOURCES})

target_link_libraries(${TENSOR_TESTS} PUBLIC
  ${STG_TENSOR_MATRIX_LIB}
  CONAN_PKG::boost
  CONAN_PKG::catch2)

target_link_directories(${TENSOR_TESTS} PUBLIC
  ${STG_STATISTICS_INCLUDE_DIR}/stg_tensor/
  CONAN_PKG::catch2
  CONAN_PKG::boost)

target_include_directories(${TENSOR_TESTS} PUBLIC
  ${STG_STATISTICS_INCLUDE_DIR}/stg_tensor/
  CONAN_PKG::catch2
  CONAN_PKG::boost)

catch_discover_tests(${TENSOR_TESTS})