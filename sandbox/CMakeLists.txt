set(STG_SANDBOX stg-sandbox)

set(STG_SANDBOX_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_SANDBOX_ROOT_DIR}/*
  ${STG_SANDBOX_INCLUDE_DIR}/*)

add_executable(${STG_SANDBOX} ${SOURCES})

target_link_libraries(${STG_SANDBOX} PUBLIC
  ${STG_TENSOR_MATRIX_LIB}
  ${STG_STATISTICS_LIB}
  CONAN_PKG::fmt

  # CONAN_LIB::fftw
  CONAN_PKG::boost
  CONAN_PKG::range-v3
  CONAN_PKG::armadillo)
target_link_directories(${STG_SANDBOX} PUBLIC
  CONAN_PKG::boost

  # CONAN_LIB::fftw
  CONAN_PKG::range-v3
  CONAN_PKG::armadillo)
target_include_directories(${STG_SANDBOX} PUBLIC
  ${STG_TENSOR_MATRIX_INCLUDE_DIR}
  ${STG_STATISTICS_INCLUDE_DIR}
  CONAN_PKG::armadillo
  CONAN_PKG::boost

  # CONAN_PKG::fftw
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)
