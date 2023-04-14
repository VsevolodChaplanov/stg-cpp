#include "common.hpp"

struct VoxelTestsFixture {
  constexpr static inline double eps = 1.e-6;
  const std::vector<double> vertices {
    0,   0,   0,
    0.5, 0,   0,
    0.5, 0.5, 0,
    0,   0.5, 0,
    0,   0,   0.5,
    0.5, 0,   0.5,
    0.5, 0.5, 0.5,
    0,   0.5, 0.5,
  };

  const std::vector<double> vertices_2 {
    0.5, 0.5, 0.5,
    1.0, 0.5, 0.5,
    1.0, 1.0, 0.5,
    0.5, 1.0, 0.5,
    0.5, 0.5, 1.0,
    1.0, 0.5, 1.0,
    1.0, 1.0, 1.0,
    0.5, 1.0, 1.0,
  };

  const std::vector<size_t> indices { 0, 1, 2, 3, 4, 5, 6, 7 };

  const VoxelFiniteElement<double> voxel_element{ vertices, indices };
  const VoxelFiniteElement<double> voxel_element_2{ vertices_2, indices };

  const std::vector<stg::Vector<double>> gradients {
    {0,    0,    0},
    {0,    0,    0},
    {0,    0,    0.25},
    {0,    0,    0},
    {0,    0,    0},
    {0,    0.25, 0},
    {0.25, 0.25, 0.25},
    {0.25, 0,    0},
  };

  const std::vector<stg::Vector<double>> gradients_2 {
    {0.25, 0.25, 0.25},
    {0.25, 0.5, 0.5},
    {0.5, 0.5, 1},
    {0.5, 0.25,0.5},
    {0.5, 0.5, 0.25},
    {0.5, 1, 0.5},
    {1, 1, 1},
    {1, 0.5, 0.5}};
};

SCENARIO_METHOD(VoxelTestsFixture, "Check logic of voxel finite element") {
  GIVEN("Assembled voxel finite element at base of CS") {
    THEN("Check for constant function") {
      std::vector<double> values(8, 1.);

      THEN("Check integrating") {
        const double volume = voxel_element.integrate(values.cbegin(), values.cend());
        CHECK_THAT(volume, WithinRel(0.125, eps));
      }

      THEN("Check gradients") {
        const auto test_gradients = voxel_element.gradients(values.cbegin(), values.cend());
        for (const auto& vec_grad: test_gradients) {
          CHECK(std::fabs(vec_grad.get<0>()) < eps);
          CHECK(std::fabs(vec_grad.get<1>()) < eps);
          CHECK(std::fabs(vec_grad.get<2>()) < eps);
        }
      }

      THEN("Check interpolation value") {
        auto interp_value = voxel_element.interpolate({0.5, 0.5, 0.5}, values);
        CHECK_THAT(interp_value, WithinRel(1., eps));
        interp_value = voxel_element.interpolate({0.3, 0.6, 0.9}, values);
        CHECK_THAT(interp_value, WithinRel(1., eps));
      }
    }

    THEN("Check for function set values") {
      std::vector<double> values{0, 0, 0, 0, 0, 0, 0.125, 0.};

      THEN("Check integrating values") {
        const double volume = voxel_element.integrate(values.cbegin(), values.cend());
        CHECK_THAT(volume, WithinRel(0.00195313, 1.e-4));
      }

      THEN("Check gradients") {
        const auto test_gradients = voxel_element.gradients(values.cbegin(), values.cend());
        REQUIRE(test_gradients.size() == gradients.size());

        for (size_t i = 0; i < 8; ++i) {
          const auto test_grad = test_gradients[i];
          const auto exact_grad = gradients[i];

          CHECK_THAT(test_grad.get<0>(), WithinRel(exact_grad.get<0>(), eps));
          CHECK_THAT(test_grad.get<1>(), WithinRel(exact_grad.get<1>(), eps));
          CHECK_THAT(test_grad.get<2>(), WithinRel(exact_grad.get<2>(), eps));
        }
      }
    }
  }

  GIVEN("Assembled voxel finite element not at base of CS") {
    THEN("Check for constant function") {
      std::vector<double> values(8, 1.);

      THEN("Check integrating") {
        const double volume = voxel_element_2.integrate(values.cbegin(), values.cend());
        CHECK_THAT(volume, WithinRel(0.125, eps));
      }

      THEN("Check gradients") {
        const auto test_gradients = voxel_element_2.gradients(values.cbegin(), values.cend());
        for (const auto& vec_grad: test_gradients) {
          CHECK(std::fabs(vec_grad.get<0>()) < eps);
          CHECK(std::fabs(vec_grad.get<1>()) < eps);
          CHECK(std::fabs(vec_grad.get<2>()) < eps);
        }
      }
    }

    THEN("Check for function set values") {
      std::vector<double> values{0.125, 0.25, 0.5, 0.25, 0.25, 0.5, 1., 0.5};

      THEN("Check integrating values") {
        const double volume = voxel_element_2.integrate(values.cbegin(), values.cend());
        CHECK_THAT(volume, WithinRel(0.0527344, 1.e-4));
      }

      THEN("Check gradients") {
        const auto test_gradients = voxel_element_2.gradients(values.cbegin(), values.cend());
        REQUIRE(test_gradients.size() == gradients.size());

        for (size_t i = 0; i < 8; ++i) {
          const auto test_grad = test_gradients[i];
          const auto exact_grad = gradients_2[i];

          CHECK_THAT(test_grad.get<0>(), WithinRel(exact_grad.get<0>(), eps));
          CHECK_THAT(test_grad.get<1>(), WithinRel(exact_grad.get<1>(), eps));
          CHECK_THAT(test_grad.get<2>(), WithinRel(exact_grad.get<2>(), eps));
        }
      }
    }
  }
}