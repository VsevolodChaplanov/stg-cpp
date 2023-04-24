#include "common.hpp"
#include <fmt/format.h>
#include <fmt/printf.h>
#include <iomanip>

struct StochasticGaussianMethod1D {
  const std::string path_str{"/home/vsevolod/coding/C++/diplom/data/kriging/covcut0/"};
  const std::string save_path_str{"/home/vsevolod/coding/C++/diplom/data/kriging/calc_covcut0/"};
  const std::array<std::size_t, 12> eigen_cuts = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600};
};

SCENARIO_METHOD(StochasticGaussianMethod1D, "Collect data from stochastic gaussian method") {
  fmt::print("|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|\n",
             "cut", "mean sqrt dev", "mean abs dev", "mean dev",
             "sum sqr dev", "sum abs dev", "sum dev");
  for (const auto eig_cut : eigen_cuts) {
    const std::string eig_cut_str = std::to_string(eig_cut);
    const std::string current_dir = "eigcut"s + eig_cut_str + "/"s;
    const std::string work_dir = path_str + current_dir;

    DataLoader data_loader{work_dir};
    KrigingAnalysis1D<double> kriging1d{data_loader};

    kriging1d.calculate_covariations();

/*    fmt::print("|{:^16.4e}|{:^16.4e}|{:^16.4e}|{:^16.4e}|{:^16.4e}|{:^16.4e}|{:^16.4e}|\n",
      eig_cut_str,
      kriging1d.mean_sqrt_deviation(),
      kriging1d.mean_abs_deviation(),
      kriging1d.mean_deviation(),
      kriging1d.sum_sqrt_deviation(),
      kriging1d.sum_abs_deviation(),
      kriging1d.sum_deviation());*/
    kriging1d.save_calculated_covariations(save_path_str + "eigcut"s + eig_cut_str + ".vtk"s);
  }
}