set(STG_GAUSSIAN_LIB stg-gaussian)

set(STG_GAUSSIAN_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_GAUSSIAN_ROOT_DIR}/src/*.cpp
  ${STG_GAUSSIAN_METHOD_INCLUDE_DIR})

add_library(${STG_GAUSSIAN_LIB} STATIC ${SOURCES})

target_link_libraries(${STG_GAUSSIAN_LIB} PUBLIC
  ${STG_MESH_LIB}
  ${STG_TENSOR_MATRIX_LIB}
  CONAN_PKG::fmt
  CONAN_PKG::fftw
  CONAN_PKG::boost
  CONAN_PKG::range-v3
  CONAN_PKG::armadillo)
target_link_directories(${STG_GAUSSIAN_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::fftw
  CONAN_PKG::range-v3
  CONAN_PKG::armadillo)
target_include_directories(${STG_GAUSSIAN_LIB} PUBLIC
  ${STG_MESH_INCLUDE_DIR}
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_GAUSSIAN_METHOD_INCLUDE_DIR}
  ${STG_GAUSSIAN_METHOD_INCLUDE_DIR}/kriging/
  ${STG_GAUSSIAN_METHOD_INCLUDE_DIR}/fourier/
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::fftw
  CONAN_PKG::range-v3
  CONAN_PKG::armadillo)

# ####################
# ##### TESTING ######
# ####################
set(GAUSSIAN_TESTS stg_gaussian_tests)

file(GLOB_RECURSE TESTSOURCES
  ${STG_GAUSSIAN_ROOT_DIR}/tests/*.cpp
  ${STG_GAUSSIAN_ROOT_DIR}/tests/*.hpp)

add_executable(${GAUSSIAN_TESTS} ${TESTSOURCES})

target_link_libraries(${GAUSSIAN_TESTS} PUBLIC
  ${STG_TENSOR_MATRIX_LIB}
  ${STG_GAUSSIAN_LIB}
  CONAN_PKG::boost
  CONAN_PKG::range-v3
  CONAN_PKG::fftw
  CONAN_PKG::catch2
  CONAN_PKG::armadillo)

target_link_directories(${GAUSSIAN_TESTS} PUBLIC
  ${STG_GAUSSIAN_METHOD_INCLUDE_DIR}/
  CONAN_PKG::catch2
  CONAN_PKG::fftw
  CONAN_PKG::range-v3
  CONAN_PKG::boost
  CONAN_PKG::armadillo)

target_include_directories(${GAUSSIAN_TESTS} PUBLIC
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_GAUSSIAN_METHOD_INCLUDE_DIR}
  ${STG_GAUSSIAN_METHOD_INCLUDE_DIR}/kriging/
  ${STG_GAUSSIAN_METHOD_INCLUDE_DIR}/fourier/
  CONAN_PKG::catch2
  CONAN_PKG::fftw
  CONAN_PKG::range-v3
  CONAN_PKG::boost
  CONAN_PKG::armadillo)

catch_discover_tests(${GAUSSIAN_TESTS})