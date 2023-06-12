#include "common.hpp"
#include "stg/spectral_method/kraichnan_method_impl.hpp"
#include "stg_thread_pool/stg_thread_pool.hpp"
#include <catch2/catch_test_macros.hpp>
#include <fmt/core.h>
#include <numbers>
#include <ranges>

struct KraichnanSpectralMethodApplicationFixture {
    const std::string directory_with_files = "./spectral_result/";
    const double time = 0;
    const double cube_edge_length = 5;
    const std::size_t cube_edge_n = 31;
    const std::size_t samples_n = 10000;
    const std::size_t fourier_n = 1000;
    const double v_0 = 1;
    const double k_0 = 2 * std::numbers::pi / cube_edge_length;
    const double w_0 = k_0 * v_0;
    const double T = 2 * std::numbers::pi / w_0;
    const double full_time = 100 * T;
    const double dt = full_time / samples_n;
};

SCENARIO_METHOD(KraichnanSpectralMethodApplicationFixture, "Generate velocity samples by Kraichnan method in time") {
    utility::ThreadPool pool{};
    for (const std::size_t isample: std::views::iota(0ull, samples_n)) {
        pool.post([isample, this] {
            const double time_moment = dt * isample;
            const std::size_t seed = isample * 42;
            KraichanMethodImpl<double> kraichnan_generator{cube_edge_length, cube_edge_n, samples_n, fourier_n, k_0, w_0, seed};
            kraichnan_generator.generate_sample(time_moment);
            kraichnan_generator.save_samples(directory_with_files,
                                             fmt::format("velocity_field_{}.vtk", isample));
        });
    }
}

SCENARIO_METHOD(KraichnanSpectralMethodApplicationFixture, "Generate velocity samples by Kraichnan method in time if generator class is same") {
    utility::ThreadPool pool{};
    // const std::size_t seed = 42;

    for (const std::size_t isample: std::views::iota(0ull, samples_n)) {
        pool.post([isample, this]() mutable {
            KraichanMethodImpl<double> kraichnan_generator{cube_edge_length, cube_edge_n, samples_n, fourier_n, k_0, w_0, 42};
            const double time_moment = dt * isample;
            kraichnan_generator.generate_sample(time_moment);
            kraichnan_generator.save_samples(directory_with_files,
                                             fmt::format("velocity_field_{}.vtk", isample));
        });
    }
}