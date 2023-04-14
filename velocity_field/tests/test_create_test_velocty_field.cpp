#include "common.hpp"

TEST_CASE("Create velocity field and check public methods", "[VelocityField]") {
  std::vector<double> vx {0, 1, 2};
  std::vector<double> vy {3, 4, 5};
  std::vector<double> vz {6, 7, 8};

  VelocityField<double> test_field{vx, vy, vz};

  REQUIRE(test_field.size() == 3);

  const auto field_view = test_field.values_view();

  const auto [vx_test, vy_test, vz_test] = field_view[1];

  CHECK_THAT(vx_test, WithinRel(vx[1], eps));
  CHECK_THAT(vy_test, WithinRel(vy[1], eps));
  CHECK_THAT(vz_test, WithinRel(vz[1], eps));

  const auto test_vx = test_field.vx_view();

  CHECK_THAT(test_vx[0], WithinRel(vx[0], eps));
  CHECK_THAT(test_vx[1], WithinRel(vx[1], eps));
  CHECK_THAT(test_vx[2], WithinRel(vx[2], eps));

  const auto first_value = test_field.value(0);

  CHECK_THAT(first_value.get<0>(), WithinRel(vx[0], eps));
  CHECK_THAT(first_value.get<1>(), WithinRel(vy[0], eps));
  CHECK_THAT(first_value.get<2>(), WithinRel(vz[0], eps));

  test_field.set_value({3, 3, 3}, 2);
  const auto third_value = test_field.value(2);

  CHECK_THAT(third_value.get<0>(), WithinRel(3, eps));
  CHECK_THAT(third_value.get<1>(), WithinRel(3, eps));
  CHECK_THAT(third_value.get<2>(), WithinRel(3, eps));
}