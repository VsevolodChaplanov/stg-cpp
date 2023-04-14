#include "common.hpp"

const double eps = 1.e-8;

SCENARIO("Geometry tests") {
  stg::Point<double> first_point{0.5, 0.5, 0.5};
  stg::Point<double> second_point{0.7, 0.7, 0.7};
  first_point += second_point;

  CHECK(first_point.get<0>() == 1.2);
  CHECK(first_point.get<1>() == 1.2);
  CHECK(first_point.get<2>() == 1.2);
}

SCENARIO("Cross product tests") {
  const stg::Point<double> first_point{0.5, -0.5, 0.5};
  const stg::Point<double> second_point{0.7, 0.7, 0.3};

  auto result = cross_product(first_point, second_point);
  CHECK_THAT(result.get<0>(), Catch::Matchers::WithinRel(-0.5, eps));
  CHECK_THAT(result.get<1>(), Catch::Matchers::WithinRel(0.2, eps));
  CHECK_THAT(result.get<2>(), Catch::Matchers::WithinRel(0.7, eps));

  const stg::Point<double> third_point{0, -0.5, 3.};
  const stg::Point<double> fourth_point{-0.7, 1.4, -2.};

  result = cross_product(third_point, fourth_point);
  CHECK_THAT(result.get<0>(), Catch::Matchers::WithinRel(-3.2, eps));
  CHECK_THAT(result.get<1>(), Catch::Matchers::WithinRel(-2.1, eps));
  CHECK_THAT(result.get<2>(), Catch::Matchers::WithinRel(-0.35, eps));
}

SCENARIO("Dot product tests") {
  const stg::Point<double> first_point{0.5, -0.5, 0.5};
  const stg::Point<double> second_point{0.7, 0.7, 0.3};

  const auto result = dot_product(first_point, second_point);
  CHECK_THAT(result, Catch::Matchers::WithinRel(0.15, eps));
}

