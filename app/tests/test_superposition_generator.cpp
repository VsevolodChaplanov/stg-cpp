#include "common.hpp"
#include "stg/spectral_method/kraichnan_method_impl.hpp"
#include "stg/spectral_method/superposition_generator.hpp"
#include "stg_generators/kraichnan_appr_generator.hpp"
#include "stg_thread_pool/stg_thread_pool.hpp"
#include <boost/asio/co_spawn.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <fmt/core.h>
#include <memory>
#include <numbers>
#include <ranges>
#include <thread>
#include <vector>

struct KraichnanSpectralMethodApplicationFixture {
    const std::string directory_with_files = "./spectral_superpos_result/";
    const double time = 0;
    const double cube_edge_length = 10;
    const std::size_t cube_edge_n = 31;
    const std::size_t samples_n = 10000;
    std::vector<NodeParameters<double>> nodes{
            NodeParameters<double>{.k_0 = std::numbers::pi / 3, .v_0 = 1., .nfourier = 1000, .seed = 42},
            NodeParameters<double>{.k_0 = std::numbers::pi * 2 / 3, .v_0 = 1., .nfourier = 1000, .seed = 42},
            NodeParameters<double>{.k_0 = std::numbers::pi, .v_0 = 1., .nfourier = 1000, .seed = 42},
            NodeParameters<double>{.k_0 = std::numbers::pi * 4 / 3, .v_0 = 1., .nfourier = 1000, .seed = 42},
            NodeParameters<double>{.k_0 = std::numbers::pi * 5 / 3, .v_0 = 1., .nfourier = 1000, .seed = 42},
    };
};

using namespace std::chrono_literals;

SCENARIO_METHOD(KraichnanSpectralMethodApplicationFixture, "Generate velocity samples by superposition method method") {
    utility::ThreadPool pool{};
    for (const std::size_t isample: std::views::iota(0ull, samples_n)) {
        pool.post([isample, this] {
            const std::size_t seed = isample * 42;
            SpectralApproximatorGenerator generator{nodes};
            auto kraichnan_generator = std::make_shared<SuperpositionGenerator<double>>(cube_edge_length, cube_edge_n, generator);
            kraichnan_generator->generate_sample(time);
            kraichnan_generator->save_sample(directory_with_files,
                                             fmt::format("velocity_field_{}.vtk", isample));
        });
    }
    fmt::print("done");
    pool.get_context().join();
}