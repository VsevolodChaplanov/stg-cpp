#include "common.hpp"
#include "stg/kriging/data_loader.hpp"
#include "stg/kriging/sequential_kriging_1d_analysis.hpp"
#include "stg/spectral_method/sequential_analysis.hpp"
#include <array>
#include <fmt/core.h>
#include <string>

struct SequentialCorrelationsTestsFixture {
    const std::string dir = "./spectral_result/";
    const std::string mesh_file = "velocity_field_0.vtk";
    const std::string save = "./spectral_result/covariations.vtk";
};

SCENARIO_METHOD(SequentialCorrelationsTestsFixture, "Calculate squentially covariation function") {

    stg::spectral::DataLoader loader{dir};
    SequentialCorrelations<double> analyser{std::move(loader), mesh_file};

    analyser.calc_covariations_for_amount(10000);
    analyser.save_calculated_covariations(save);
}

SCENARIO_METHOD(SequentialCorrelationsTestsFixture, "Calculate correlations and spectra for Kraichnan method") {
    stg::spectral::DataLoader loader{dir};
    SequentialCorrelations<double> analyser{std::move(loader), "covariations.vtk", SequentialCorrelations<double>::WithoutGeneration{}};
    const auto mesh_size = analyser.size();
    CHECK(mesh_size > 1);


    analyser.convert_covariations_to_correlations();
    analyser.save_calculated_correlations("./spectral_result/correlations.vtk");


    auto file_with_integrals = fmt::output_file("./spectral_result/energies.csv",
                                                fmt::file::CREATE | fmt::file::WRONLY | fmt::file::APPEND);

    std::array<double, 50> k_moduluses{};
    double k = 0;
    double dk = 0.5;
    for (auto& k_mod: k_moduluses) {
        k_mod = k;
        k += dk;
    }
    file_with_integrals.print("[");
    for (const auto k_mod: k_moduluses) {
        const auto [Ex, Ey, Ez] = analyser.energy_for_wavenumber(k_mod, 30, 30);
        file_with_integrals.print("[{},{},{}]\n,", Ex, Ey, Ez);
    }
    file_with_integrals.print("]");
}