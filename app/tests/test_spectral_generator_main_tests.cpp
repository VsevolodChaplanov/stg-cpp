#include "common.hpp"

struct SpectralAppFixture {
  const double edge_len = 2.;
  const size_t edge_points_n = 11;
  const std::string filename = fmt::format("tests_mesh/spectral_method_{:.1}_{}.vtk", edge_len, edge_points_n);

  SpectralMethodApplicationImpl<double> spectral_generation_method{
    edge_len, edge_points_n, mock::MakeDefaultSpectralConfig(), 100, 4};
};

SCENARIO_METHOD(SpectralAppFixture, "Simple generation procedure to check algos") {
  spectral_generation_method.generate_on_mesh(0);
  spectral_generation_method.generate_samples_on_mesh(0);
  spectral_generation_method.generate_correlations_relative_to_space_center();

  spectral_generation_method.save_to_vtk(filename);
}
