#include "common.hpp"
#include "stg/spectral_method/des_generator_impl.hpp"
#include "stg/spectral_method/kraichnan_method_impl.hpp"
#include "stg/spectral_method/sequential_analysis.hpp"
#include "stg_thread_pool/stg_thread_pool.hpp"
#include <boost/asio/co_spawn.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <fmt/core.h>
#include <memory>
#include <numbers>
#include <ranges>
#include <thread>

struct KraichnanSpectralMethodApplicationFixture {
    const std::string directory_with_files = "./des_spectral_result/";
    const double time = 0;
    const double cube_edge_length = 9 * pi;
    const std::size_t cube_edge_n = 64;
    const std::size_t samples_n = 4000;
    const std::size_t fourier_n = 5000;
    const double k_0 = std::numbers::pi / 3;
    const double v_0 = 1.;
    const double w_0 = 1.;
};

double EnergySpectra(double k) {
    double logkappa = log10(k);
    double logE;
    if (logkappa < 0.0) {
        logE = 2 * logkappa - 1;
    } else if (logkappa < 3.0) {
        logE = -5.0 / 3.0 * logkappa - 1;
    } else {
        logE = -3 * logkappa + 3;
    }
    return std::pow(10, logE);
    // constexpr auto alpha = 1.453;
    // constexpr auto nu = 1.;
    // constexpr auto usub = 0.25;
    // const auto k_e = 40 * std::sqrt(5. / 12.);
    // const auto L = 0.746834 / k_e;
    // const auto epsilon = usub * usub * usub / L;
    // const auto k_eta = std::pow(epsilon, 1. / 4.) * std::pow(nu, -3. / 4.);
    // const auto expvalue = std::pow(std::numbers::e, -2. * (k / k_e) * (k / k_eta));
    // const auto numerator = std::pow(k / k_e, 4);
    // const auto denumerator = std::pow(1 + k / k_e, 17. / 6.);

    // const auto E = alpha * usub * usub / k_e * numerator / denumerator * expvalue;
    // return E;
}

using namespace std::chrono_literals;

SCENARIO_METHOD(KraichnanSpectralMethodApplicationFixture, "Generate velocity samples by Direct energy spectral method") {
    utility::ThreadPool pool{};

    for (const std::size_t isample: std::views::iota(0ull, samples_n)) {
        pool.post([isample, this] {
            const std::size_t seed = isample * 42;
            thread_local const auto energy_lambda = [](double k) {
                double logkappa = log10(k);
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
            auto des_generator = std::make_shared<DesGeneratorImpl<double>>(cube_edge_length, cube_edge_n, samples_n, fourier_n, w_0, std::move(energy_lambda), seed);
            des_generator->generate_sample(time);
            des_generator->save_sample(directory_with_files,
                                       fmt::format("velocity_field_{}.vtk", isample));
        });
    }
    fmt::print("done");
    pool.get_context().join();
}

struct SequentialCorrelationsTestsDESMethodFixture {
    const std::string dir = "./des_spectral_result/";
    const std::string mesh_file = "velocity_field_0.vtk";
    const std::string save = "./des_spectral_result/covariations.vtk";
};

SCENARIO_METHOD(SequentialCorrelationsTestsDESMethodFixture, "Calculate squentially covariation function for direct energy spectra method") {

    stg::spectral::DataLoader loader{dir};
    SequentialCorrelations<double> analyser{std::move(loader), mesh_file};

    analyser.calc_covariations_for_amount(4000);
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