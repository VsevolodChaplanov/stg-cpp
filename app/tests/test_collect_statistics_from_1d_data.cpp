#include "common.hpp"

struct StochasticGaussianMethod1D {
  const std::string path_str{"test_resources/1d_kriging/"};
  const std::array<std::size_t, 6> eigen_cuts = {100, 120, 140, 160, 180, 200};
};

SCENARIO_METHOD(StochasticGaussianMethod1D, "Collect data from stochastic gaussian method") {
  for (const auto eig_cut : eigen_cuts) {
    const std::string eig_cut_str = std::to_string(eig_cut);
    const std::string current_dir = "eigcut"s + eig_cut_str + "/"s;
    const std::string work_dir = path_str + current_dir;

    DataLoader data_loader{work_dir};
    KrigingAnalysis1D<double> kriging1d{data_loader};

    kriging1d.calculate_covariations();

    std::cout << kriging1d.std_btw_covariations() << std::endl;
  }
}