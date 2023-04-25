#include "common.hpp"
#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/ostream.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/std.h>
#include <iomanip>
#include <fstream>

struct StochasticGaussianMethod1D {
  const std::size_t cov_cut_percent = 5;
  const std::string path_str{fmt::format("/home/vsevolod/coding/kriging/covcut{}/", cov_cut_percent)};
  const std::string save_path_str{fmt::format("/home/vsevolod/coding/kriging/calc_covcut{}/", cov_cut_percent)};
  const std::array<std::size_t, 11> eigen_cuts = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600};
};

SCENARIO_METHOD(StochasticGaussianMethod1D, "Collect data from stochastic gaussian method") {
  auto file_with_integrals = fmt::output_file("/home/vsevolod/coding/kriging/integrals.txt",
                                              fmt::file::CREATE | fmt::file::WRONLY | fmt::file::APPEND);

  file_with_integrals.print("covariation cut {}\n", cov_cut_percent);
  file_with_integrals.print(
    "|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}\n",
    "cut", "mean sqrt dev", "mean abs dev", "mean dev",
    "sum sqr dev", "sum abs dev", "sum dev", "peak diff", "int sqr diff");
  for (const auto eig_cut : eigen_cuts) {
    const std::string eig_cut_str = std::to_string(eig_cut);
    const std::string current_dir = "eigcut"s + eig_cut_str + "/"s;
    const std::string work_dir = path_str + current_dir;

    DataLoader data_loader{work_dir};
    KrigingAnalysis1D<double> kriging1d{data_loader};

    kriging1d.calculate_covariations();

    file_with_integrals.print(
      "|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}\n",
      eig_cut_str,
      kriging1d.mean_sqrt_deviation(),
      kriging1d.mean_abs_deviation(),
      kriging1d.mean_deviation(),
      kriging1d.sum_sqrt_deviation(),
      kriging1d.sum_abs_deviation(),
      kriging1d.sum_deviation(),
      kriging1d.peaks_deviation(),
      kriging1d.integrate_square_diff());

    const std::string save_filename = save_path_str + "eigcut"s + eig_cut_str + ".vtk"s;
    kriging1d.save_calculated_covariations(save_filename);
  }
}