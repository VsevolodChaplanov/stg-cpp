
#################
###  BOOST
#################
set(CONAN_BOOST_ROOT_RELEASE "/home/vsevolod/.conan/data/boost/1.80.0/_/_/package/83fcf0af66d2d8b7b7d60ef6c249b03c5b0dc744")
set(CONAN_INCLUDE_DIRS_BOOST_RELEASE "/home/vsevolod/.conan/data/boost/1.80.0/_/_/package/83fcf0af66d2d8b7b7d60ef6c249b03c5b0dc744/include")
set(CONAN_LIB_DIRS_BOOST_RELEASE "/home/vsevolod/.conan/data/boost/1.80.0/_/_/package/83fcf0af66d2d8b7b7d60ef6c249b03c5b0dc744/lib")
set(CONAN_BIN_DIRS_BOOST_RELEASE )
set(CONAN_RES_DIRS_BOOST_RELEASE )
set(CONAN_SRC_DIRS_BOOST_RELEASE )
set(CONAN_BUILD_DIRS_BOOST_RELEASE )
set(CONAN_FRAMEWORK_DIRS_BOOST_RELEASE )
set(CONAN_LIBS_BOOST_RELEASE boost_contract boost_coroutine boost_fiber_numa boost_fiber boost_context boost_graph boost_iostreams boost_json boost_locale boost_log_setup boost_log boost_math_c99 boost_math_c99f boost_math_c99l boost_math_tr1 boost_math_tr1f boost_math_tr1l boost_nowide boost_program_options boost_random boost_regex boost_stacktrace_addr2line boost_stacktrace_backtrace boost_stacktrace_basic boost_stacktrace_noop boost_timer boost_type_erasure boost_thread boost_chrono boost_container boost_date_time boost_unit_test_framework boost_prg_exec_monitor boost_test_exec_monitor boost_exception boost_wave boost_filesystem boost_atomic boost_wserialization boost_serialization)
set(CONAN_PKG_LIBS_BOOST_RELEASE boost_contract boost_coroutine boost_fiber_numa boost_fiber boost_context boost_graph boost_iostreams boost_json boost_locale boost_log_setup boost_log boost_math_c99 boost_math_c99f boost_math_c99l boost_math_tr1 boost_math_tr1f boost_math_tr1l boost_nowide boost_program_options boost_random boost_regex boost_stacktrace_addr2line boost_stacktrace_backtrace boost_stacktrace_basic boost_stacktrace_noop boost_timer boost_type_erasure boost_thread boost_chrono boost_container boost_date_time boost_unit_test_framework boost_prg_exec_monitor boost_test_exec_monitor boost_exception boost_wave boost_filesystem boost_atomic boost_wserialization boost_serialization)
set(CONAN_SYSTEM_LIBS_BOOST_RELEASE dl rt pthread)
set(CONAN_FRAMEWORKS_BOOST_RELEASE )
set(CONAN_FRAMEWORKS_FOUND_BOOST_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_BOOST_RELEASE "-DBOOST_STACKTRACE_ADDR2LINE_LOCATION=\"/usr/bin/addr2line\""
			"-DBOOST_STACKTRACE_USE_ADDR2LINE"
			"-DBOOST_STACKTRACE_USE_BACKTRACE"
			"-DBOOST_STACKTRACE_USE_NOOP")
set(CONAN_BUILD_MODULES_PATHS_BOOST_RELEASE )
# COMPILE_DEFINITIONS are equal to CONAN_DEFINES without -D, for targets
set(CONAN_COMPILE_DEFINITIONS_BOOST_RELEASE "BOOST_STACKTRACE_ADDR2LINE_LOCATION=\"/usr/bin/addr2line\""
			"BOOST_STACKTRACE_USE_ADDR2LINE"
			"BOOST_STACKTRACE_USE_BACKTRACE"
			"BOOST_STACKTRACE_USE_NOOP")

set(CONAN_C_FLAGS_BOOST_RELEASE "")
set(CONAN_CXX_FLAGS_BOOST_RELEASE "")
set(CONAN_SHARED_LINKER_FLAGS_BOOST_RELEASE "")
set(CONAN_EXE_LINKER_FLAGS_BOOST_RELEASE "")

# For modern cmake targets we use the list variables (separated with ;)
set(CONAN_C_FLAGS_BOOST_RELEASE_LIST "")
set(CONAN_CXX_FLAGS_BOOST_RELEASE_LIST "")
set(CONAN_SHARED_LINKER_FLAGS_BOOST_RELEASE_LIST "")
set(CONAN_EXE_LINKER_FLAGS_BOOST_RELEASE_LIST "")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_BOOST_RELEASE "${CONAN_FRAMEWORKS_BOOST_RELEASE}" "_BOOST" "_RELEASE")
# Append to aggregated values variable
set(CONAN_LIBS_BOOST_RELEASE ${CONAN_PKG_LIBS_BOOST_RELEASE} ${CONAN_SYSTEM_LIBS_BOOST_RELEASE} ${CONAN_FRAMEWORKS_FOUND_BOOST_RELEASE})


#################
###  CATCH2
#################
set(CONAN_CATCH2_ROOT_RELEASE "/home/vsevolod/.conan/data/catch2/3.1.0/_/_/package/9950e1afe939202602a2de8a1b8a986781821358")
set(CONAN_INCLUDE_DIRS_CATCH2_RELEASE "/home/vsevolod/.conan/data/catch2/3.1.0/_/_/package/9950e1afe939202602a2de8a1b8a986781821358/include")
set(CONAN_LIB_DIRS_CATCH2_RELEASE "/home/vsevolod/.conan/data/catch2/3.1.0/_/_/package/9950e1afe939202602a2de8a1b8a986781821358/lib")
set(CONAN_BIN_DIRS_CATCH2_RELEASE )
set(CONAN_RES_DIRS_CATCH2_RELEASE )
set(CONAN_SRC_DIRS_CATCH2_RELEASE )
set(CONAN_BUILD_DIRS_CATCH2_RELEASE "/home/vsevolod/.conan/data/catch2/3.1.0/_/_/package/9950e1afe939202602a2de8a1b8a986781821358/lib/cmake/Catch2")
set(CONAN_FRAMEWORK_DIRS_CATCH2_RELEASE )
set(CONAN_LIBS_CATCH2_RELEASE Catch2Main Catch2)
set(CONAN_PKG_LIBS_CATCH2_RELEASE Catch2Main Catch2)
set(CONAN_SYSTEM_LIBS_CATCH2_RELEASE m)
set(CONAN_FRAMEWORKS_CATCH2_RELEASE )
set(CONAN_FRAMEWORKS_FOUND_CATCH2_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_CATCH2_RELEASE )
set(CONAN_BUILD_MODULES_PATHS_CATCH2_RELEASE )
# COMPILE_DEFINITIONS are equal to CONAN_DEFINES without -D, for targets
set(CONAN_COMPILE_DEFINITIONS_CATCH2_RELEASE )

set(CONAN_C_FLAGS_CATCH2_RELEASE "")
set(CONAN_CXX_FLAGS_CATCH2_RELEASE "")
set(CONAN_SHARED_LINKER_FLAGS_CATCH2_RELEASE "")
set(CONAN_EXE_LINKER_FLAGS_CATCH2_RELEASE "")

# For modern cmake targets we use the list variables (separated with ;)
set(CONAN_C_FLAGS_CATCH2_RELEASE_LIST "")
set(CONAN_CXX_FLAGS_CATCH2_RELEASE_LIST "")
set(CONAN_SHARED_LINKER_FLAGS_CATCH2_RELEASE_LIST "")
set(CONAN_EXE_LINKER_FLAGS_CATCH2_RELEASE_LIST "")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_CATCH2_RELEASE "${CONAN_FRAMEWORKS_CATCH2_RELEASE}" "_CATCH2" "_RELEASE")
# Append to aggregated values variable
set(CONAN_LIBS_CATCH2_RELEASE ${CONAN_PKG_LIBS_CATCH2_RELEASE} ${CONAN_SYSTEM_LIBS_CATCH2_RELEASE} ${CONAN_FRAMEWORKS_FOUND_CATCH2_RELEASE})


#################
###  FMT
#################
set(CONAN_FMT_ROOT_RELEASE "/home/vsevolod/.conan/data/fmt/9.1.0/_/_/package/300cb578a6fb15627a390ac70b832e82fd040283")
set(CONAN_INCLUDE_DIRS_FMT_RELEASE "/home/vsevolod/.conan/data/fmt/9.1.0/_/_/package/300cb578a6fb15627a390ac70b832e82fd040283/include")
set(CONAN_LIB_DIRS_FMT_RELEASE "/home/vsevolod/.conan/data/fmt/9.1.0/_/_/package/300cb578a6fb15627a390ac70b832e82fd040283/lib")
set(CONAN_BIN_DIRS_FMT_RELEASE )
set(CONAN_RES_DIRS_FMT_RELEASE )
set(CONAN_SRC_DIRS_FMT_RELEASE )
set(CONAN_BUILD_DIRS_FMT_RELEASE )
set(CONAN_FRAMEWORK_DIRS_FMT_RELEASE )
set(CONAN_LIBS_FMT_RELEASE fmt)
set(CONAN_PKG_LIBS_FMT_RELEASE fmt)
set(CONAN_SYSTEM_LIBS_FMT_RELEASE m)
set(CONAN_FRAMEWORKS_FMT_RELEASE )
set(CONAN_FRAMEWORKS_FOUND_FMT_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_FMT_RELEASE )
set(CONAN_BUILD_MODULES_PATHS_FMT_RELEASE )
# COMPILE_DEFINITIONS are equal to CONAN_DEFINES without -D, for targets
set(CONAN_COMPILE_DEFINITIONS_FMT_RELEASE )

set(CONAN_C_FLAGS_FMT_RELEASE "")
set(CONAN_CXX_FLAGS_FMT_RELEASE "")
set(CONAN_SHARED_LINKER_FLAGS_FMT_RELEASE "")
set(CONAN_EXE_LINKER_FLAGS_FMT_RELEASE "")

# For modern cmake targets we use the list variables (separated with ;)
set(CONAN_C_FLAGS_FMT_RELEASE_LIST "")
set(CONAN_CXX_FLAGS_FMT_RELEASE_LIST "")
set(CONAN_SHARED_LINKER_FLAGS_FMT_RELEASE_LIST "")
set(CONAN_EXE_LINKER_FLAGS_FMT_RELEASE_LIST "")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_FMT_RELEASE "${CONAN_FRAMEWORKS_FMT_RELEASE}" "_FMT" "_RELEASE")
# Append to aggregated values variable
set(CONAN_LIBS_FMT_RELEASE ${CONAN_PKG_LIBS_FMT_RELEASE} ${CONAN_SYSTEM_LIBS_FMT_RELEASE} ${CONAN_FRAMEWORKS_FOUND_FMT_RELEASE})


#################
###  RANGE-V3
#################
set(CONAN_RANGE-V3_ROOT_RELEASE "/home/vsevolod/.conan/data/range-v3/0.12.0/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9")
set(CONAN_INCLUDE_DIRS_RANGE-V3_RELEASE "/home/vsevolod/.conan/data/range-v3/0.12.0/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/include")
set(CONAN_LIB_DIRS_RANGE-V3_RELEASE )
set(CONAN_BIN_DIRS_RANGE-V3_RELEASE )
set(CONAN_RES_DIRS_RANGE-V3_RELEASE )
set(CONAN_SRC_DIRS_RANGE-V3_RELEASE )
set(CONAN_BUILD_DIRS_RANGE-V3_RELEASE )
set(CONAN_FRAMEWORK_DIRS_RANGE-V3_RELEASE )
set(CONAN_LIBS_RANGE-V3_RELEASE )
set(CONAN_PKG_LIBS_RANGE-V3_RELEASE )
set(CONAN_SYSTEM_LIBS_RANGE-V3_RELEASE )
set(CONAN_FRAMEWORKS_RANGE-V3_RELEASE )
set(CONAN_FRAMEWORKS_FOUND_RANGE-V3_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_RANGE-V3_RELEASE )
set(CONAN_BUILD_MODULES_PATHS_RANGE-V3_RELEASE )
# COMPILE_DEFINITIONS are equal to CONAN_DEFINES without -D, for targets
set(CONAN_COMPILE_DEFINITIONS_RANGE-V3_RELEASE )

set(CONAN_C_FLAGS_RANGE-V3_RELEASE "")
set(CONAN_CXX_FLAGS_RANGE-V3_RELEASE "")
set(CONAN_SHARED_LINKER_FLAGS_RANGE-V3_RELEASE "")
set(CONAN_EXE_LINKER_FLAGS_RANGE-V3_RELEASE "")

# For modern cmake targets we use the list variables (separated with ;)
set(CONAN_C_FLAGS_RANGE-V3_RELEASE_LIST "")
set(CONAN_CXX_FLAGS_RANGE-V3_RELEASE_LIST "")
set(CONAN_SHARED_LINKER_FLAGS_RANGE-V3_RELEASE_LIST "")
set(CONAN_EXE_LINKER_FLAGS_RANGE-V3_RELEASE_LIST "")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_RANGE-V3_RELEASE "${CONAN_FRAMEWORKS_RANGE-V3_RELEASE}" "_RANGE-V3" "_RELEASE")
# Append to aggregated values variable
set(CONAN_LIBS_RANGE-V3_RELEASE ${CONAN_PKG_LIBS_RANGE-V3_RELEASE} ${CONAN_SYSTEM_LIBS_RANGE-V3_RELEASE} ${CONAN_FRAMEWORKS_FOUND_RANGE-V3_RELEASE})


#################
###  FFTW
#################
set(CONAN_FFTW_ROOT_RELEASE "/home/vsevolod/.conan/data/fftw/3.3.10/_/_/package/113f6963bbd263d1f9a23fe9431f2476563a0626")
set(CONAN_INCLUDE_DIRS_FFTW_RELEASE "/home/vsevolod/.conan/data/fftw/3.3.10/_/_/package/113f6963bbd263d1f9a23fe9431f2476563a0626/include")
set(CONAN_LIB_DIRS_FFTW_RELEASE "/home/vsevolod/.conan/data/fftw/3.3.10/_/_/package/113f6963bbd263d1f9a23fe9431f2476563a0626/lib")
set(CONAN_BIN_DIRS_FFTW_RELEASE )
set(CONAN_RES_DIRS_FFTW_RELEASE )
set(CONAN_SRC_DIRS_FFTW_RELEASE )
set(CONAN_BUILD_DIRS_FFTW_RELEASE )
set(CONAN_FRAMEWORK_DIRS_FFTW_RELEASE )
set(CONAN_LIBS_FFTW_RELEASE fftw3_omp fftw3)
set(CONAN_PKG_LIBS_FFTW_RELEASE fftw3_omp fftw3)
set(CONAN_SYSTEM_LIBS_FFTW_RELEASE m pthread)
set(CONAN_FRAMEWORKS_FFTW_RELEASE )
set(CONAN_FRAMEWORKS_FOUND_FFTW_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_FFTW_RELEASE )
set(CONAN_BUILD_MODULES_PATHS_FFTW_RELEASE )
# COMPILE_DEFINITIONS are equal to CONAN_DEFINES without -D, for targets
set(CONAN_COMPILE_DEFINITIONS_FFTW_RELEASE )

set(CONAN_C_FLAGS_FFTW_RELEASE "")
set(CONAN_CXX_FLAGS_FFTW_RELEASE "")
set(CONAN_SHARED_LINKER_FLAGS_FFTW_RELEASE "")
set(CONAN_EXE_LINKER_FLAGS_FFTW_RELEASE "")

# For modern cmake targets we use the list variables (separated with ;)
set(CONAN_C_FLAGS_FFTW_RELEASE_LIST "")
set(CONAN_CXX_FLAGS_FFTW_RELEASE_LIST "")
set(CONAN_SHARED_LINKER_FLAGS_FFTW_RELEASE_LIST "")
set(CONAN_EXE_LINKER_FLAGS_FFTW_RELEASE_LIST "")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_FFTW_RELEASE "${CONAN_FRAMEWORKS_FFTW_RELEASE}" "_FFTW" "_RELEASE")
# Append to aggregated values variable
set(CONAN_LIBS_FFTW_RELEASE ${CONAN_PKG_LIBS_FFTW_RELEASE} ${CONAN_SYSTEM_LIBS_FFTW_RELEASE} ${CONAN_FRAMEWORKS_FOUND_FFTW_RELEASE})


#################
###  ZLIB
#################
set(CONAN_ZLIB_ROOT_RELEASE "/home/vsevolod/.conan/data/zlib/1.2.13/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6")
set(CONAN_INCLUDE_DIRS_ZLIB_RELEASE "/home/vsevolod/.conan/data/zlib/1.2.13/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/include")
set(CONAN_LIB_DIRS_ZLIB_RELEASE "/home/vsevolod/.conan/data/zlib/1.2.13/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/lib")
set(CONAN_BIN_DIRS_ZLIB_RELEASE )
set(CONAN_RES_DIRS_ZLIB_RELEASE )
set(CONAN_SRC_DIRS_ZLIB_RELEASE )
set(CONAN_BUILD_DIRS_ZLIB_RELEASE "/home/vsevolod/.conan/data/zlib/1.2.13/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/")
set(CONAN_FRAMEWORK_DIRS_ZLIB_RELEASE )
set(CONAN_LIBS_ZLIB_RELEASE z)
set(CONAN_PKG_LIBS_ZLIB_RELEASE z)
set(CONAN_SYSTEM_LIBS_ZLIB_RELEASE )
set(CONAN_FRAMEWORKS_ZLIB_RELEASE )
set(CONAN_FRAMEWORKS_FOUND_ZLIB_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_ZLIB_RELEASE )
set(CONAN_BUILD_MODULES_PATHS_ZLIB_RELEASE )
# COMPILE_DEFINITIONS are equal to CONAN_DEFINES without -D, for targets
set(CONAN_COMPILE_DEFINITIONS_ZLIB_RELEASE )

set(CONAN_C_FLAGS_ZLIB_RELEASE "")
set(CONAN_CXX_FLAGS_ZLIB_RELEASE "")
set(CONAN_SHARED_LINKER_FLAGS_ZLIB_RELEASE "")
set(CONAN_EXE_LINKER_FLAGS_ZLIB_RELEASE "")

# For modern cmake targets we use the list variables (separated with ;)
set(CONAN_C_FLAGS_ZLIB_RELEASE_LIST "")
set(CONAN_CXX_FLAGS_ZLIB_RELEASE_LIST "")
set(CONAN_SHARED_LINKER_FLAGS_ZLIB_RELEASE_LIST "")
set(CONAN_EXE_LINKER_FLAGS_ZLIB_RELEASE_LIST "")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_ZLIB_RELEASE "${CONAN_FRAMEWORKS_ZLIB_RELEASE}" "_ZLIB" "_RELEASE")
# Append to aggregated values variable
set(CONAN_LIBS_ZLIB_RELEASE ${CONAN_PKG_LIBS_ZLIB_RELEASE} ${CONAN_SYSTEM_LIBS_ZLIB_RELEASE} ${CONAN_FRAMEWORKS_FOUND_ZLIB_RELEASE})


#################
###  BZIP2
#################
set(CONAN_BZIP2_ROOT_RELEASE "/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2")
set(CONAN_INCLUDE_DIRS_BZIP2_RELEASE "/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2/include")
set(CONAN_LIB_DIRS_BZIP2_RELEASE "/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2/lib")
set(CONAN_BIN_DIRS_BZIP2_RELEASE "/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2/bin")
set(CONAN_RES_DIRS_BZIP2_RELEASE )
set(CONAN_SRC_DIRS_BZIP2_RELEASE )
set(CONAN_BUILD_DIRS_BZIP2_RELEASE "/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2/")
set(CONAN_FRAMEWORK_DIRS_BZIP2_RELEASE )
set(CONAN_LIBS_BZIP2_RELEASE bz2)
set(CONAN_PKG_LIBS_BZIP2_RELEASE bz2)
set(CONAN_SYSTEM_LIBS_BZIP2_RELEASE )
set(CONAN_FRAMEWORKS_BZIP2_RELEASE )
set(CONAN_FRAMEWORKS_FOUND_BZIP2_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_BZIP2_RELEASE )
set(CONAN_BUILD_MODULES_PATHS_BZIP2_RELEASE )
# COMPILE_DEFINITIONS are equal to CONAN_DEFINES without -D, for targets
set(CONAN_COMPILE_DEFINITIONS_BZIP2_RELEASE )

set(CONAN_C_FLAGS_BZIP2_RELEASE "")
set(CONAN_CXX_FLAGS_BZIP2_RELEASE "")
set(CONAN_SHARED_LINKER_FLAGS_BZIP2_RELEASE "")
set(CONAN_EXE_LINKER_FLAGS_BZIP2_RELEASE "")

# For modern cmake targets we use the list variables (separated with ;)
set(CONAN_C_FLAGS_BZIP2_RELEASE_LIST "")
set(CONAN_CXX_FLAGS_BZIP2_RELEASE_LIST "")
set(CONAN_SHARED_LINKER_FLAGS_BZIP2_RELEASE_LIST "")
set(CONAN_EXE_LINKER_FLAGS_BZIP2_RELEASE_LIST "")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_BZIP2_RELEASE "${CONAN_FRAMEWORKS_BZIP2_RELEASE}" "_BZIP2" "_RELEASE")
# Append to aggregated values variable
set(CONAN_LIBS_BZIP2_RELEASE ${CONAN_PKG_LIBS_BZIP2_RELEASE} ${CONAN_SYSTEM_LIBS_BZIP2_RELEASE} ${CONAN_FRAMEWORKS_FOUND_BZIP2_RELEASE})


#################
###  LIBBACKTRACE
#################
set(CONAN_LIBBACKTRACE_ROOT_RELEASE "/home/vsevolod/.conan/data/libbacktrace/cci.20210118/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6")
set(CONAN_INCLUDE_DIRS_LIBBACKTRACE_RELEASE "/home/vsevolod/.conan/data/libbacktrace/cci.20210118/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/include")
set(CONAN_LIB_DIRS_LIBBACKTRACE_RELEASE "/home/vsevolod/.conan/data/libbacktrace/cci.20210118/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/lib")
set(CONAN_BIN_DIRS_LIBBACKTRACE_RELEASE )
set(CONAN_RES_DIRS_LIBBACKTRACE_RELEASE )
set(CONAN_SRC_DIRS_LIBBACKTRACE_RELEASE )
set(CONAN_BUILD_DIRS_LIBBACKTRACE_RELEASE "/home/vsevolod/.conan/data/libbacktrace/cci.20210118/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/")
set(CONAN_FRAMEWORK_DIRS_LIBBACKTRACE_RELEASE )
set(CONAN_LIBS_LIBBACKTRACE_RELEASE backtrace)
set(CONAN_PKG_LIBS_LIBBACKTRACE_RELEASE backtrace)
set(CONAN_SYSTEM_LIBS_LIBBACKTRACE_RELEASE )
set(CONAN_FRAMEWORKS_LIBBACKTRACE_RELEASE )
set(CONAN_FRAMEWORKS_FOUND_LIBBACKTRACE_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_LIBBACKTRACE_RELEASE )
set(CONAN_BUILD_MODULES_PATHS_LIBBACKTRACE_RELEASE )
# COMPILE_DEFINITIONS are equal to CONAN_DEFINES without -D, for targets
set(CONAN_COMPILE_DEFINITIONS_LIBBACKTRACE_RELEASE )

set(CONAN_C_FLAGS_LIBBACKTRACE_RELEASE "")
set(CONAN_CXX_FLAGS_LIBBACKTRACE_RELEASE "")
set(CONAN_SHARED_LINKER_FLAGS_LIBBACKTRACE_RELEASE "")
set(CONAN_EXE_LINKER_FLAGS_LIBBACKTRACE_RELEASE "")

# For modern cmake targets we use the list variables (separated with ;)
set(CONAN_C_FLAGS_LIBBACKTRACE_RELEASE_LIST "")
set(CONAN_CXX_FLAGS_LIBBACKTRACE_RELEASE_LIST "")
set(CONAN_SHARED_LINKER_FLAGS_LIBBACKTRACE_RELEASE_LIST "")
set(CONAN_EXE_LINKER_FLAGS_LIBBACKTRACE_RELEASE_LIST "")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_LIBBACKTRACE_RELEASE "${CONAN_FRAMEWORKS_LIBBACKTRACE_RELEASE}" "_LIBBACKTRACE" "_RELEASE")
# Append to aggregated values variable
set(CONAN_LIBS_LIBBACKTRACE_RELEASE ${CONAN_PKG_LIBS_LIBBACKTRACE_RELEASE} ${CONAN_SYSTEM_LIBS_LIBBACKTRACE_RELEASE} ${CONAN_FRAMEWORKS_FOUND_LIBBACKTRACE_RELEASE})


### Definition of global aggregated variables ###

set(CONAN_DEPENDENCIES_RELEASE boost catch2 fmt range-v3 fftw zlib bzip2 libbacktrace)

set(CONAN_INCLUDE_DIRS_RELEASE "/home/vsevolod/.conan/data/boost/1.80.0/_/_/package/83fcf0af66d2d8b7b7d60ef6c249b03c5b0dc744/include"
			"/home/vsevolod/.conan/data/catch2/3.1.0/_/_/package/9950e1afe939202602a2de8a1b8a986781821358/include"
			"/home/vsevolod/.conan/data/fmt/9.1.0/_/_/package/300cb578a6fb15627a390ac70b832e82fd040283/include"
			"/home/vsevolod/.conan/data/range-v3/0.12.0/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/include"
			"/home/vsevolod/.conan/data/fftw/3.3.10/_/_/package/113f6963bbd263d1f9a23fe9431f2476563a0626/include"
			"/home/vsevolod/.conan/data/zlib/1.2.13/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/include"
			"/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2/include"
			"/home/vsevolod/.conan/data/libbacktrace/cci.20210118/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/include" ${CONAN_INCLUDE_DIRS_RELEASE})
set(CONAN_LIB_DIRS_RELEASE "/home/vsevolod/.conan/data/boost/1.80.0/_/_/package/83fcf0af66d2d8b7b7d60ef6c249b03c5b0dc744/lib"
			"/home/vsevolod/.conan/data/catch2/3.1.0/_/_/package/9950e1afe939202602a2de8a1b8a986781821358/lib"
			"/home/vsevolod/.conan/data/fmt/9.1.0/_/_/package/300cb578a6fb15627a390ac70b832e82fd040283/lib"
			"/home/vsevolod/.conan/data/fftw/3.3.10/_/_/package/113f6963bbd263d1f9a23fe9431f2476563a0626/lib"
			"/home/vsevolod/.conan/data/zlib/1.2.13/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/lib"
			"/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2/lib"
			"/home/vsevolod/.conan/data/libbacktrace/cci.20210118/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/lib" ${CONAN_LIB_DIRS_RELEASE})
set(CONAN_BIN_DIRS_RELEASE "/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2/bin" ${CONAN_BIN_DIRS_RELEASE})
set(CONAN_RES_DIRS_RELEASE  ${CONAN_RES_DIRS_RELEASE})
set(CONAN_FRAMEWORK_DIRS_RELEASE  ${CONAN_FRAMEWORK_DIRS_RELEASE})
set(CONAN_LIBS_RELEASE boost_contract boost_coroutine boost_fiber_numa boost_fiber boost_context boost_graph boost_iostreams boost_json boost_locale boost_log_setup boost_log boost_math_c99 boost_math_c99f boost_math_c99l boost_math_tr1 boost_math_tr1f boost_math_tr1l boost_nowide boost_program_options boost_random boost_regex boost_stacktrace_addr2line boost_stacktrace_backtrace boost_stacktrace_basic boost_stacktrace_noop boost_timer boost_type_erasure boost_thread boost_chrono boost_container boost_date_time boost_unit_test_framework boost_prg_exec_monitor boost_test_exec_monitor boost_exception boost_wave boost_filesystem boost_atomic boost_wserialization boost_serialization Catch2Main Catch2 fmt fftw3_omp fftw3 z bz2 backtrace ${CONAN_LIBS_RELEASE})
set(CONAN_PKG_LIBS_RELEASE boost_contract boost_coroutine boost_fiber_numa boost_fiber boost_context boost_graph boost_iostreams boost_json boost_locale boost_log_setup boost_log boost_math_c99 boost_math_c99f boost_math_c99l boost_math_tr1 boost_math_tr1f boost_math_tr1l boost_nowide boost_program_options boost_random boost_regex boost_stacktrace_addr2line boost_stacktrace_backtrace boost_stacktrace_basic boost_stacktrace_noop boost_timer boost_type_erasure boost_thread boost_chrono boost_container boost_date_time boost_unit_test_framework boost_prg_exec_monitor boost_test_exec_monitor boost_exception boost_wave boost_filesystem boost_atomic boost_wserialization boost_serialization Catch2Main Catch2 fmt fftw3_omp fftw3 z bz2 backtrace ${CONAN_PKG_LIBS_RELEASE})
set(CONAN_SYSTEM_LIBS_RELEASE dl rt m pthread ${CONAN_SYSTEM_LIBS_RELEASE})
set(CONAN_FRAMEWORKS_RELEASE  ${CONAN_FRAMEWORKS_RELEASE})
set(CONAN_FRAMEWORKS_FOUND_RELEASE "")  # Will be filled later
set(CONAN_DEFINES_RELEASE "-DBOOST_STACKTRACE_ADDR2LINE_LOCATION=\"/usr/bin/addr2line\""
			"-DBOOST_STACKTRACE_USE_ADDR2LINE"
			"-DBOOST_STACKTRACE_USE_BACKTRACE"
			"-DBOOST_STACKTRACE_USE_NOOP" ${CONAN_DEFINES_RELEASE})
set(CONAN_BUILD_MODULES_PATHS_RELEASE  ${CONAN_BUILD_MODULES_PATHS_RELEASE})
set(CONAN_CMAKE_MODULE_PATH_RELEASE "/home/vsevolod/.conan/data/catch2/3.1.0/_/_/package/9950e1afe939202602a2de8a1b8a986781821358/lib/cmake/Catch2"
			"/home/vsevolod/.conan/data/zlib/1.2.13/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/"
			"/home/vsevolod/.conan/data/bzip2/1.0.8/_/_/package/3cfc45772763dad1237052f26c1fe8b2bae3f7d2/"
			"/home/vsevolod/.conan/data/libbacktrace/cci.20210118/_/_/package/2a19826344ff00be1c04403f2f8e7008ed3a7cc6/" ${CONAN_CMAKE_MODULE_PATH_RELEASE})

set(CONAN_CXX_FLAGS_RELEASE " ${CONAN_CXX_FLAGS_RELEASE}")
set(CONAN_SHARED_LINKER_FLAGS_RELEASE " ${CONAN_SHARED_LINKER_FLAGS_RELEASE}")
set(CONAN_EXE_LINKER_FLAGS_RELEASE " ${CONAN_EXE_LINKER_FLAGS_RELEASE}")
set(CONAN_C_FLAGS_RELEASE " ${CONAN_C_FLAGS_RELEASE}")

# Apple Frameworks
conan_find_apple_frameworks(CONAN_FRAMEWORKS_FOUND_RELEASE "${CONAN_FRAMEWORKS_RELEASE}" "" "_RELEASE")
# Append to aggregated values variable: Use CONAN_LIBS instead of CONAN_PKG_LIBS to include user appended vars
set(CONAN_LIBS_RELEASE ${CONAN_LIBS_RELEASE} ${CONAN_SYSTEM_LIBS_RELEASE} ${CONAN_FRAMEWORKS_FOUND_RELEASE})
