#ifndef STG_MESH_COMMON_HPP
#define STG_MESH_COMMON_HPP

#include <random>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <fem.hpp>
#include <geometry/geometry.hpp>
#include <velocity_field/velocity_field.hpp>

using namespace stg;
using namespace stg::mesh;
using namespace stg::field;
using namespace Catch::Matchers;

static const double eps = 1.e-6;

#endif //STG_COMMON_HPP
