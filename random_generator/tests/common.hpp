#ifndef STG_STG_RANDOM_COMMON_HPP
#define STG_STG_RANDOM_COMMON_HPP

#include <iostream>
#include <memory>
#include <stg_random/generator_engines.hpp>
#include <stg_random/irn_generator.hpp>
#include <stg_random/generator_concept.hpp>
#include <stg_random/rn_generator_impl.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// random number generator seed
constexpr std::size_t seed = 42;

#endif