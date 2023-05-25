set(STG_UTILITY_LIB stg-utility)

set(STG_UTILITY_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES
  ${STG_UTILITY_ROOT_DIR}/src/*.cpp
  ${STG_UTILITY_ROOT_DIR}/include/*)

add_library(${STG_UTILITY_LIB} STATIC ${SOURCES})

target_link_libraries(${STG_UTILITY_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)
target_link_directories(${STG_UTILITY_LIB} PUBLIC
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)
target_include_directories(${STG_UTILITY_LIB} PUBLIC
  ${STG_UTILITY_INCLUDE_DIR}
  ${STG_UTILITY_INCLUDE_DIR}/stg_coro_future/
  ${STG_UTILITY_INCLUDE_DIR}/stg_thread_pool/
  CONAN_PKG::boost
  CONAN_PKG::fmt
  CONAN_PKG::range-v3)

# ####################
# ##### TESTING ######
# ####################
set(UTILITY_TESTS stg-utility-tests)

file(GLOB_RECURSE TESTSOURCES
  ${STG_UTILITY_ROOT_DIR}/tests/*.cpp
  ${STG_UTILITY_ROOT_DIR}/tests/*.hpp)

add_executable(${UTILITY_TESTS} ${TESTSOURCES})

target_link_libraries(${UTILITY_TESTS} PUBLIC
  ${STG_UTILITY_LIB}
  CONAN_PKG::boost
  CONAN_PKG::range-v3
  CONAN_PKG::catch2)

target_link_directories(${UTILITY_TESTS} PUBLIC
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

target_include_directories(${UTILITY_TESTS} PUBLIC
  ${STG_UTILITY_INCLUDE_DIR}
  CONAN_PKG::catch2
  CONAN_PKG::range-v3
  CONAN_PKG::boost)

catch_discover_tests(${UTILITY_TESTS})