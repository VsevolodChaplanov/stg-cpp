add_test( [==[Scenario: Creating fourier space mesh tests]==] [==[/home/vsevolod/coding/C++/diplom/std_v3-cpp/cmake-build-release/bin/stg_gaussian_tests]==] [==[Scenario: Creating fourier space mesh tests]==]  )
set_tests_properties( [==[Scenario: Creating fourier space mesh tests]==] PROPERTIES WORKING_DIRECTORY [==[/home/vsevolod/coding/C++/diplom/std_v3-cpp/cmake-build-release/gaussian_method]==])
add_test( [==[Scenario: Direct fourier transform tests]==] [==[/home/vsevolod/coding/C++/diplom/std_v3-cpp/cmake-build-release/bin/stg_gaussian_tests]==] [==[Scenario: Direct fourier transform tests]==]  )
set_tests_properties( [==[Scenario: Direct fourier transform tests]==] PROPERTIES WORKING_DIRECTORY [==[/home/vsevolod/coding/C++/diplom/std_v3-cpp/cmake-build-release/gaussian_method]==])
set( stg_gaussian_tests_TESTS [==[Scenario: Creating fourier space mesh tests]==] [==[Scenario: Direct fourier transform tests]==])