#include "common.hpp"
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/os.h>
#include <fmt/ostream.h>
#include <fmt/printf.h>
#include <fmt/std.h>
#include <fstream>
#include <iomanip>

struct StochasticGaussianMethod1D {
    const std::string path_str{fmt::format("/home/vsevolod/coding/kriging/new_params/")};
    const std::string save_path_str{fmt::format("/home/vsevolod/coding/kriging/calc_new_params/")};
};

SCENARIO_METHOD(StochasticGaussianMethod1D, "Collect data from stochastic gaussian method in single folder") {
    auto file_with_integrals = fmt::output_file("/home/vsevolod/coding/kriging/integrals.txt",
                                                fmt::file::CREATE | fmt::file::WRONLY | fmt::file::APPEND);

    file_with_integrals.print("covariation cut {}\n and max_iter = 10000, method = la", 0.01);
    file_with_integrals.print(
            "|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}|{:^16}\n",
            "cut", "mean sqrt dev", "mean abs dev", "mean dev",
            "sum sqr dev", "sum abs dev", "sum dev", "peak diff", "int sqr diff");

    stg::spectral::DataLoader data_loader{std::filesystem::path(path_str)};
    KrigingAnalysis1D<double> kriging1d{data_loader};

    kriging1d.calculate_covariations();

    file_with_integrals.print(
            "|{:^16}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}|{:^16.4}\n",
            600,
            kriging1d.mean_sqrt_deviation(),
            kriging1d.mean_abs_deviation(),
            kriging1d.mean_deviation(),
            kriging1d.sum_sqrt_deviation(),
            kriging1d.sum_abs_deviation(),
            kriging1d.sum_deviation(),
            kriging1d.peaks_deviation(),
            kriging1d.integrate_square_diff());

    const std::string save_filename = save_path_str + "eigcut600.vtk"s;
    kriging1d.save_calculated_covariations(save_filename);
}