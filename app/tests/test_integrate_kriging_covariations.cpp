#include "common.hpp"
#include "stg/spectral_method/data_loader.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <fmt/core.h>
#include <fmt/os.h>
#include <ranges>


struct IntegratorFixture {
    const std::string work_dir = "./kriging_result/";
    const std::string cov_filename = "r_11.vtk";
    const std::string fert_filename = "phi_11.vtk";
};


SCENARIO_METHOD(IntegratorFixture, "Integrate covariations obtained form kriging method") {
    IntegrateCovariations1D<double> integrator{stg::spectral::DataLoader{work_dir}, cov_filename};

    const double k_start = 0;
    const double k_end = 5;
    const std::size_t n = 10;
    const double dk = (k_end - k_start) / (n - 1);

    auto out = fmt::output_file(fmt::format("{}{}", work_dir, "energies.csv"), fmt::file::CREATE | fmt::file::WRONLY | fmt::file::APPEND);

    for (const std::size_t i: std::views::iota(0ul, n)) {
        const double k_mod = k_start + dk * i;
        const double e_value = integrator.integrate_for(k_mod);

        out.print("{},{},\n", k_mod, e_value);
    }
}

SCENARIO_METHOD(IntegratorFixture, "Integrate covariation and obtain same fert values obtained form kriging method") {
    IntegrateFert1D<double> fert_integrator{stg::spectral::DataLoader{work_dir}, cov_filename, fert_filename};

    fert_integrator.integrate_covariations();
    fert_integrator.save_calculated_fert(work_dir + "fert_11.vtk"s);
}

SCENARIO_METHOD(IntegratorFixture, "Integrate covariation and obtain same fert values obtained form kriging method and check with zero r") {
    IntegrateFert1D<double> fert_integrator{stg::spectral::DataLoader{work_dir}, cov_filename, fert_filename};

    const auto [fert, r] = fert_integrator.integrate_pure();
    CHECK_THAT(r, WithinRel(fert, 1.e-4));
}

SCENARIO_METHOD(IntegratorFixture, "Integrate covariation and obtain same fert values obtained form kriging method using fourier tranform") {
    IntegrateFert1D<double> fert_integrator{stg::spectral::DataLoader{work_dir}, cov_filename, fert_filename};

    const auto [fert, r] = fert_integrator.apply_fourier_to_covariations();
    fert_integrator.save_calculated_fert(work_dir + "fert_11_fourier.vtk"s);
    CHECK_THAT(r, WithinRel(fert, 1.e-4));
}