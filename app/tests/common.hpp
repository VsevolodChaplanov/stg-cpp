#ifndef STG_MESH_COMMON_HPP
#define STG_MESH_COMMON_HPP

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_container_properties.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <fem.hpp>
#include <geometry/geometry.hpp>
#include <mesh_builders.hpp>
#include <mocks/mock_spectral.hpp>
#include <random>
#include <spherical_mesh/sphere_mesh.hpp>
#include <stg/kriging.hpp>
#include <stg/spectral_method.hpp>

namespace fs = std::filesystem;
using namespace std::literals;
using namespace Catch::Matchers;
using namespace stg;
using namespace stg::mesh;
using namespace stg::tensor;
using namespace stg::spectral;
using namespace stg::kriging;

#endif//STG_COMMON_HPP
