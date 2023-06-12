#include "common.hpp"
#include "stg/kriging/sequential_kriging_1d_analysis.hpp"
#include "stg/spectral_method/data_loader.hpp"
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

    analyser.convert_covariations_to_correlations();
    analyser.save_calculated_correlations(dir + "correlations.vtk");

    analyser.make_correlations_symmetric();
    analyser.save_calculated_symm_correlations(dir + "correlations_symm.vtk");

    analyser.calculate_energies_by_symm_correlations<0, 0>(dir + "fert_11_s.vtk", dir + "energy_11_s.vtk");
    analyser.calculate_energies_by_symm_correlations<0, 0>(dir + "fert_22_s.vtk", dir + "energy_22_s.vtk");
    analyser.calculate_energies_by_symm_correlations<0, 0>(dir + "fert_33_s.vtk", dir + "energy_33_s.vtk");


    analyser.make_covariance_symmetric();
    analyser.save_calculated_symm_covariance(dir + "covariance_symm.vtk");
    analyser.calculate_energies_by_symm_covariance<0, 0>(dir + "fert_11_s_cov.vtk", dir + "energy_11_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<0, 1>(dir + "fert_12_s_cov.vtk", dir + "energy_12_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<0, 2>(dir + "fert_13_s_cov.vtk", dir + "energy_13_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<1, 0>(dir + "fert_21_s_cov.vtk", dir + "energy_21_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<1, 1>(dir + "fert_22_s_cov.vtk", dir + "energy_22_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<1, 2>(dir + "fert_23_s_cov.vtk", dir + "energy_23_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<2, 0>(dir + "fert_31_s_cov.vtk", dir + "energy_31_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<2, 1>(dir + "fert_32_s_cov.vtk", dir + "energy_32_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<2, 2>(dir + "fert_33_s_cov.vtk", dir + "energy_33_s_cov.vtk");
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

struct SequentialCorrelationsTestsSuperpositionMethodFixture {
    const std::string dir = "./spectral_superpos_result/";
    const std::string mesh_file = "velocity_field_0.vtk";
    const std::string save = "./spectral_superpos_result/covariations.vtk";
};

SCENARIO_METHOD(SequentialCorrelationsTestsSuperpositionMethodFixture, "Calculate squentially covariation function") {

    stg::spectral::DataLoader loader{dir};
    SequentialCorrelations<double> analyser{std::move(loader), mesh_file};

    analyser.calc_covariations_for_amount(10000);
    analyser.save_calculated_covariations(save);

    analyser.convert_covariations_to_correlations();
    analyser.save_calculated_correlations(dir + "correlations.vtk");

    analyser.make_correlations_symmetric();
    analyser.save_calculated_symm_correlations(dir + "correlations_symm.vtk");

    analyser.calculate_energies_by_symm_correlations<0, 0>(dir + "fert_11_s.vtk", dir + "energy_11_s.vtk");
    analyser.calculate_energies_by_symm_correlations<0, 0>(dir + "fert_22_s.vtk", dir + "energy_22_s.vtk");
    analyser.calculate_energies_by_symm_correlations<0, 0>(dir + "fert_33_s.vtk", dir + "energy_33_s.vtk");


    analyser.make_covariance_symmetric();
    analyser.save_calculated_symm_covariance(dir + "covariance_symm.vtk");
    analyser.calculate_energies_by_symm_covariance<0, 0>(dir + "fert_11_s_cov.vtk", dir + "energy_11_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<0, 1>(dir + "fert_12_s_cov.vtk", dir + "energy_12_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<0, 2>(dir + "fert_13_s_cov.vtk", dir + "energy_13_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<1, 0>(dir + "fert_21_s_cov.vtk", dir + "energy_21_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<1, 1>(dir + "fert_22_s_cov.vtk", dir + "energy_22_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<1, 2>(dir + "fert_23_s_cov.vtk", dir + "energy_23_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<2, 0>(dir + "fert_31_s_cov.vtk", dir + "energy_31_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<2, 1>(dir + "fert_32_s_cov.vtk", dir + "energy_32_s_cov.vtk");
    analyser.calculate_energies_by_symm_covariance<2, 2>(dir + "fert_33_s_cov.vtk", dir + "energy_33_s_cov.vtk");
}