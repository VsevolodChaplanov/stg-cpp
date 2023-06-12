#include "common.hpp"
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/os.h>
#include <fmt/ostream.h>
#include <fmt/printf.h>
#include <fmt/std.h>
#include <fstream>
#include <iomanip>

struct SequantialAnalysisStochasticGaussianMethod1D {
    Params params{8000, 21, 10, 10000};
    const std::string path_str{fmt::format("/home/vsevolod/coding/kriging/kriging_1d/eigcut0_amount100000_alleigvals")};
    const std::string save_path_str{fmt::format("/home/vsevolod/coding/kriging/kriging_1d/eigcut0_amount100000_alleigvals_calculated/")};
};

SCENARIO_METHOD(SequantialAnalysisStochasticGaussianMethod1D, "Do analysis in sequence") {
    auto file_with_integrals = fmt::output_file("/home/vsevolod/coding/kriging/kriging_1d/eigcut0_amount100000_alleigvals_calculated/integrals_seq.txt",
                                                fmt::file::CREATE | fmt::file::WRONLY | fmt::file::APPEND);
    file_with_integrals.print("covariation cut {}\n and max_iter = 10000, method = la", 0.00);
    file_with_integrals.print(
            "{:^16}\t&\t{:^16}\t&\t{:^16}\n",
            "samples", "peak diff", "int sqr diff");

    for (const std::size_t amount: {100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000}) {
        stg::spectral::DataLoader data_loader{path_str};
        SequentialAnalysis<double> kriging1d{std::move(data_loader), params};
        kriging1d.calc_covariations_for_amount(amount);
        file_with_integrals.print(
                "{:^16}\t&\t{:^16}\t&\t{:^16}\n",
                amount,
                kriging1d.peaks_deviation(),
                kriging1d.integrate_square_diff());

        const std::string save_filename = save_path_str + fmt::format("cov_samples_amount{}.vtk", amount);
        kriging1d.save_calculated_covariations(save_filename);
    }
}