#include "common.hpp"

struct SphericalIntegrationFixture {
  const double eps = 1.e-4;

  const double r = 1.;
  const std::size_t n_theta = 1000;
  const std::size_t n_phi = 1000;
  SphericalMesh<double> spherical_mesh{ r, n_theta, n_phi };
  const std::vector<Point<double>> central_points = spherical_mesh.elements_centers();

  static const auto inline sym_func = [] (Point<double> point) {
    const double x = point.get<0>();
    const double y = point.get<1>();
    const double z = point.get<2>();

    return x * y * z;
  };

  static const auto inline sym_func_2 = [] (Point<double> point) {
    const double x = point.get<0>();
    const double y = point.get<1>();
    const double z = point.get<2>();

    return x * x + y * y + z * z;
  };

  static const auto inline non_sym_func = [] (Point<double> point) {
    const double x = point.get<0>();
    const double y = point.get<1>();
    const double z = point.get<2>();

    return x * x * x + y * y + z;
  };
};

SCENARIO_METHOD(SphericalIntegrationFixture, "Spherical integration routine tests") {
  GIVEN("Spherical mesh") {

    const size_t size = central_points.size();

    REQUIRE(size == (n_phi - 1) *  (n_theta - 1));

    THEN("Integrate 1 and get sphere square") {
      std::vector<double> values_to_integrate(size, 1);

      const auto integral = spherical_mesh.integrate(values_to_integrate.cbegin(), values_to_integrate.cend());
      CHECK_THAT(integral, Catch::Matchers::WithinRel(4 * std::numbers::pi * r * r, eps));
    }

    THEN("Integrate symmetric function to get 0") {
      std::vector<double> values_to_integrate;
      values_to_integrate.reserve(size);
      for (const size_t index : ranges::views::iota(0ul, size)) {
        values_to_integrate.push_back(sym_func(central_points[index]));
      }

      const auto integral = spherical_mesh.integrate(values_to_integrate.cbegin(), values_to_integrate.cend());
      CHECK_THAT(integral, Catch::Matchers::WithinAbs(0, eps));
    }

    THEN("Integrate sphere function and get sphere square") {
      std::vector<double> values_to_integrate;
      values_to_integrate.reserve(size);
      for (const size_t index : ranges::views::iota(0ul, size)) {
        values_to_integrate.push_back(sym_func_2(central_points[index]));
      }

      const auto integral = spherical_mesh.integrate(values_to_integrate.cbegin(), values_to_integrate.cend());
      CHECK_THAT(integral, Catch::Matchers::WithinRel(4 * std::numbers::pi * r * r, eps));
    }

    THEN("Integrate non symmetric function") {
      std::vector<double> values_to_integrate;
      values_to_integrate.reserve(size);
      for (const size_t index : ranges::views::iota(0ul, size)) {
        values_to_integrate.push_back(non_sym_func(central_points[index]));
      }

      const auto integral = spherical_mesh.integrate(values_to_integrate.cbegin(), values_to_integrate.cend());
      CHECK_THAT(integral, Catch::Matchers::WithinRel(4.18879, eps));
    }
  }
}