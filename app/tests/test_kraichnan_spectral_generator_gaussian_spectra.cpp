#include "common.hpp"
#include "stg/spectral_method/kraichnan_method_impl.hpp"
#include "stg_thread_pool/stg_thread_pool.hpp"
#include <catch2/catch_test_macros.hpp>
#include <fmt/core.h>
#include <numbers>
#include <ranges>

struct KraichnanSpectralMethodApplicationGaussianSpectraFixture {
    const std::string directory_with_files = "./spectral_result/";
    const double time = 0;
    const double cube_edge_length = 10;
    const std::size_t cube_edge_n = 21;
    const std::size_t samples_n = 10000;
    const std::size_t fourier_n = 1000;
    const double k_0 = std::numbers::pi / 3;
    const double w_0 = 1;
};

SCENARIO_METHOD(KraichnanSpectralMethodApplicationGaussianSpectraFixture, "Generate velocity samples by Kraichnan method with gaussian spectra") {
    utility::ThreadPool pool{};
    for (const std::size_t isample: std::views::iota(0ull, samples_n)) {
        pool.post([isample, this] {
            const std::size_t seed = isample * 42;
            KraichanMethodGaussianSpectraImpl<double> kraichnan_generator{cube_edge_length, cube_edge_n, samples_n, fourier_n, k_0, w_0, seed};
            kraichnan_generator.generate_sample(time);
            kraichnan_generator.save_samples(directory_with_files,
                                             fmt::format("velocity_field_{}.vtk", isample));
        });
    }
    pool.get_context().join();
}