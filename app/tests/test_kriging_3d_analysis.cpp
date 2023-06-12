#include "common.hpp"
#include "stg/kriging/kriging_3d_analysis.hpp"
#include "stg/spectral_method/data_loader.hpp"
#include "stg/spectral_method/sequential_analysis.hpp"
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <fmt/core.h>
#include <stg/spectral_method/validation_data.hpp>
#include <string>

struct Kriging3DCorrelationsTestsFixture {
    const double l = 10;
    const std::size_t n = 17;
    const std::size_t eigcut = 1000;
    const std::string dir = "./kriging_result/";
    const std::string cov_mesh_file = "r_11.vtk";
    const std::string fert_mesh_file = "phi_11.vtk";
};

SCENARIO_METHOD(Kriging3DCorrelationsTestsFixture, "Calculate statistic data for kriging simulation") {
    stg::spectral::DataLoader loader{dir};
    Kriging3DAnalysis<double> analyser{std::move(loader), cov_mesh_file, fert_mesh_file, l, n, eigcut};

    analyser.calc_covariations_for_amount(10000);
    analyser.save_calculated_covariations(dir + "covariance.vtk");

    // analyser.convert_covariations_to_correlations();
    // analyser.save_calculated_correlations(dir + "correlations.vtk");

    // const auto E_xx = analyser.calculate_energies<0, 0>(dir + "fert_11.vtk", dir + "energy_11.vtk");
    // const auto E_yy = analyser.calculate_energies<1, 1>(dir + "fert_22.vtk", dir + "energy_22.vtk");
    // const auto E_zz = analyser.calculate_energies<2, 2>(dir + "fert_33.vtk", dir + "energy_33.vtk");

    analyser.make_correlations_symmetric();
    analyser.save_calculated_symm_correlations(dir + "symm_correlations.vtk");
    const auto E_xx_s = analyser.calculate_energies_by_symm_correlations<0, 0>(dir + "fert_11_s.vtk", dir + "energy_11_s.vtk");
    const auto E_yy_s = analyser.calculate_energies_by_symm_correlations<1, 1>(dir + "fert_22_s.vtk", dir + "energy_22_s.vtk");
    const auto E_zz_s = analyser.calculate_energies_by_symm_correlations<2, 2>(dir + "fert_33_s.vtk", dir + "energy_33_s.vtk");
}

SCENARIO_METHOD(Kriging3DCorrelationsTestsFixture, "Calculate energy tensor and enregy function data for kriging simulation") {
    stg::spectral::DataLoader loader{dir};
    Kriging3DAnalysis<double> analyser{std::move(loader), cov_mesh_file, fert_mesh_file, l, n, eigcut};

    analyser.calc_covariations_for_amount(10000);
    analyser.save_calculated_covariations(dir + "covariance.vtk");

    analyser.convert_covariations_to_correlations();
    analyser.save_calculated_correlations(dir + "correlations.vtk");

    analyser.make_correlations_symmetric();
    analyser.save_calculated_symm_correlations(dir + "symm_correlations.vtk");
    const auto E_xx_s = analyser.calculate_energies_by_symm_correlations<0, 0>(dir + "fert_11_s.vtk", dir + "energy_11_s.vtk");
    const auto E_yy_s = analyser.calculate_energies_by_symm_correlations<1, 1>(dir + "fert_22_s.vtk", dir + "energy_22_s.vtk");
    const auto E_zz_s = analyser.calculate_energies_by_symm_correlations<2, 2>(dir + "fert_33_s.vtk", dir + "energy_33_s.vtk");

    analyser.make_covariance_symmetric();
    analyser.save_calculated_symm_covariance(dir + "symm_covariance.vtk");
    const auto E_xx_s_cov = analyser.calculate_energies_by_symm_covariance<0, 0>(dir + "fert_11_s_cov.vtk", dir + "energy_11_s_cov.vtk");
    const auto E_xy_s_cov = analyser.calculate_energies_by_symm_covariance<0, 1>(dir + "fert_12_s_cov.vtk", dir + "energy_12_s_cov.vtk");
    const auto E_xz_s_cov = analyser.calculate_energies_by_symm_covariance<0, 2>(dir + "fert_13_s_cov.vtk", dir + "energy_13_s_cov.vtk");
    const auto E_yx_s_cov = analyser.calculate_energies_by_symm_covariance<1, 0>(dir + "fert_21_s_cov.vtk", dir + "energy_21_s_cov.vtk");
    const auto E_yy_s_cov = analyser.calculate_energies_by_symm_covariance<1, 1>(dir + "fert_22_s_cov.vtk", dir + "energy_22_s_cov.vtk");
    const auto E_yz_s_cov = analyser.calculate_energies_by_symm_covariance<1, 2>(dir + "fert_23_s_cov.vtk", dir + "energy_23_s_cov.vtk");
    const auto E_zx_s_cov = analyser.calculate_energies_by_symm_covariance<2, 0>(dir + "fert_31_s_cov.vtk", dir + "energy_31_s_cov.vtk");
    const auto E_zy_s_cov = analyser.calculate_energies_by_symm_covariance<2, 1>(dir + "fert_32_s_cov.vtk", dir + "energy_32_s_cov.vtk");
    const auto E_zz_s_cov = analyser.calculate_energies_by_symm_covariance<2, 2>(dir + "fert_33_s_cov.vtk", dir + "energy_33_s_cov.vtk");
}

double E(double kappa) {
    double logkappa = log10(kappa);
    double logE;
    if (logkappa < 0.0) {
        logE = 2 * logkappa - 1;
    } else if (logkappa < 3.0) {
        logE = -5.0 / 3.0 * logkappa - 1;
    } else {
        logE = -3 * logkappa + 3;
    }
    return std::pow(10, logE);
};

SCENARIO_METHOD(Kriging3DCorrelationsTestsFixture, "Create exact energy disctrete function") {
    DataLoader loader{dir};
    EnergyFunction<double> energy_function{loader, "phi_11.vtk", &E};
    energy_function.fill();
    energy_function.save_to_file(dir + "energy_function.vtk");
}