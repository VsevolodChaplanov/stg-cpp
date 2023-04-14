#ifndef STG_STATISTICS_COMMON_HPP
#define STG_STATISTICS_COMMON_HPP

#include <random>
#include <array>
#include <stg_tensor/tensor.hpp>
#include <stg_rangom.hpp>
#include <statistics.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/catch_approx.hpp>
#include <range/v3/all.hpp>

using namespace stg::tensor;
using namespace stg::statistics;
using namespace Catch::Matchers;

// random number generator seed
constexpr std::size_t seed = 42;

#endif