#ifndef STG_MESH_COMMON_HPP
#define STG_MESH_COMMON_HPP

#include <random>
#include <array>
#include <range/v3/experimental/utility/generator.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <spherical_mesh/sphere_mesh.hpp>
#include <fem.hpp>
#include <mesh_builders.hpp>
#include <vtk_parser/rectilinear_grid_parser.hpp>
#include <geometry/geometry.hpp>

using namespace stg;
using namespace stg::mesh;
using namespace stg::tensor;
using namespace Catch::Matchers;

#endif //STG_COMMON_HPP
