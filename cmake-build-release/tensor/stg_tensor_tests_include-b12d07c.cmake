if(EXISTS "/home/vsevolod/coding/stg-cpp/cmake-build-release/tensor/stg_tensor_tests_tests-b12d07c.cmake")
  include("/home/vsevolod/coding/stg-cpp/cmake-build-release/tensor/stg_tensor_tests_tests-b12d07c.cmake")
else()
  add_test(stg_tensor_tests_NOT_BUILT-b12d07c stg_tensor_tests_NOT_BUILT-b12d07c)
endif()
