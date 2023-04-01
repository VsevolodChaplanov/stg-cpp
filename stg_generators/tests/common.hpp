#ifndef STG_STG_RANDOM_COMMON_HPP
#define STG_STG_RANDOM_COMMON_HPP

#include <iostream>
#include <stg_generators.hpp>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// random number generator seed
constexpr std::size_t seed = 42;
using namespace stg;
using namespace stg::generators;
using namespace stg::generator_engines;
#endif