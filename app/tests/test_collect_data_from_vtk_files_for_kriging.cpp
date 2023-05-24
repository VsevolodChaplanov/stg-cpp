#include "common.hpp"

struct DataLoadingFixture {
  static inline const double eps = 1.e-6;
  const Tensor<double> exact_first_cov_tensor = {
    -9.19963e-05,
    -8.35613e-05,
    -8.35613e-05,
    -8.35613e-05,
    -9.19963e-05,
    -8.35613e-05,
    -8.35613e-05,
    -8.35613e-05,
    -9.19963e-05
  };

  const Tensor<double> exact_first_phi_tensor = {
    3.69811e-07,
   -1.84906e-07,
   -1.84906e-07,
   -1.84906e-07,
    3.69811e-07,
   -1.84906e-07,
   -1.84906e-07,
   -1.84906e-07,
    3.69811e-07,
  };

  static inline const auto floating_point_eq = [] (auto it1, auto it2, auto it3, auto it4) {
    REQUIRE(std::distance(it1, it2) == std::distance(it3, it4));
    for (; it1 != it2; ++it1, ++it3) {
      CHECK_THAT(*it1, WithinRel(*it3, 1.e-6));
    }
  };
};

SCENARIO_METHOD(DataLoadingFixture, "Load generated values from vtk files") {
  stg::kriging::DataLoader data_loader{"test_resources/"};
  const auto real_space_mesh = data_loader.load_mesh<double>("r_11.vtk");
  const auto covariation_data = data_loader.load_covariation_data<double>();
  const auto velocity_samples = data_loader.load_velocity_samples<double>();
  const auto phi_data = data_loader.load_fert<double>();

  CHECK(real_space_mesh->n_vertices() == 51 * 51 * 51);
  CHECK(real_space_mesh->n_elements() == 50 * 50 * 50);

  const auto& first_cov_tensor = covariation_data[0];
  floating_point_eq(first_cov_tensor.cbegin(), first_cov_tensor.cend(),
                    exact_first_cov_tensor.cbegin(), exact_first_cov_tensor.cend());
  const auto& first_phi_data = phi_data[0];
  floating_point_eq(first_phi_data.cbegin(), first_phi_data.cend(),
                    exact_first_phi_tensor.cbegin(), exact_first_phi_tensor.cend());
}