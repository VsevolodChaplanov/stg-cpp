#include "common.hpp"
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <fmt/os.h>
#include <ranges>

double energy_function(double kappa) {

    double k_0 = 1;
    double sigma_k = 0.2;
    return std::pow(std::numbers::e, -1. / 2. * (kappa - k_0) * (kappa - k_0) / (sigma_k * sigma_k));
    // double logkappa = log10(kappa);
    // double logE;
    // if (logkappa < 0.0) {
    //     logE = 2 * logkappa - 1;
    // } else if (logkappa < 3.0) {
    //     logE = -5.0 / 3.0 * logkappa - 1;
    // } else {
    //     logE = -3 * logkappa + 3;
    // }
    // return std::pow(10, logE);
};

SCENARIO("Save spectra function") {
    auto file = fmt::output_file("energy.csv", fmt::file::CREATE | fmt::file::WRONLY);

    const double k_start = 0;
    const double k_end = 10;
    const std::size_t n = 100;
    const double dk = (k_end - k_start) / (n - 1);

    for (const auto value: std::views::iota(0ul, n) | std::views::transform([&](auto index) mutable {
                               return k_start + dk * index;
                           })) {
        file.print("{}, ", energy_function(value));
    }
}