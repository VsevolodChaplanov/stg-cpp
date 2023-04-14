#include "common.hpp"

struct DirectFFTTestFixture {
  const double eps = 1.e-6;
  const std::size_t n = 11;

  const double l = 2;

  CubeMeshBuilder<double> builder{l, n};
  std::shared_ptr<CubeRelationTable<double>> space = builder.build_relation_table();

  const VtkSaver saver;
  const std::string filename
    = fmt::format("tests_mesh/vtk_fourier_space_{}_{}_{}.vtk", n, n, n);
};

SCENARIO_METHOD(DirectFFTTestFixture, "Creating fourier space mesh tests") {

}

SCENARIO_METHOD(DirectFFTTestFixture, "Direct fourier transform tests") {
  const auto func = [](double x, double y, double z) {
    double v1 = 1. / (std::abs(x) + 1);
    double v2 = 1. / (std::abs(y) + 1);
    double v3 = 1. / (std::abs(z) + 1);

    return std::pow(v1, 6) * std::pow(v2, 4) * std::pow(v3, 5);
  };

  GIVEN("Values to get fourier image") {
    std::vector<double> values(n * n * n);
    for (size_t k = 0; k < n; ++k) {
      for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i) {
          const auto vertex = space->vertex(i, j, k);
          values[space->lin_index(i, j, k)] = func(vertex.get<0>(), vertex.get<1>(), vertex.get<2>());
        }
      }
    }

    std::vector<double> fourier_image = fourier3(space, values);

    CHECK_THAT(fourier_image[5 + 5 * n + 5 * n * n], WithinRel(0.0005079819819736192, eps));
    CHECK_THAT(fourier_image[6 + 5 * n + 5 * n * n], WithinRel(0.00036019289899672674, eps));
    CHECK_THAT(fourier_image[7 + 5 * n + 5 * n * n], WithinRel(0.000216620583015753, eps));
    CHECK_THAT(fourier_image[4 + 4 * n + 4 * n * n], WithinRel(0.00011965752189560347, eps));
    CHECK_THAT(fourier_image[4 + 5 * n + 7 * n * n], WithinRel(0.0001205415104388087, eps));
  }
}
