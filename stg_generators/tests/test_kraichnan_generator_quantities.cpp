#include "common.hpp"
#include "geometry/geometry.hpp"
#include "stg_generators/kraichnan_spectral_generator.hpp"
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>
#include <functional>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/view/transform.hpp>

struct KraichnanGeneratorTestsFixture {
    static constexpr inline std::size_t nsamples = 10000;
    static constexpr inline std::size_t init_seed = 42;

    const std::size_t nfourier = 1000;

    const double k0 = 1.;
    const double w0 = 1.;
    const double time = 0;

    const Vector<double> zero_point{1, 1, 1};
};

SCENARIO_METHOD(KraichnanGeneratorTestsFixture, "Kraichnan generator statistical parameters") {
    std::array<Vector<double>, nsamples> generated_values{};
    for (std::size_t i = 0; auto& value: generated_values) {
        const auto seed = init_seed * i;
        KraichanGeneratorDeltaFunction<double> generator{nfourier, k0, w0, seed};
        auto&& generated_vector = generator(zero_point, time);
        value = std::move(generated_vector);
        i++;
    }

    Vector<double> mean{0, 0, 0};
    for (const auto value: generated_values) {
        mean += value;
    }
    mean = mean / static_cast<double>(nsamples);

    CHECK_THAT(mean.get<0>(), Catch::Matchers::WithinRel(-18.265085507537123, 1.e-4));
    CHECK_THAT(mean.get<1>(), Catch::Matchers::WithinRel(-4.6558832116794804, 1.e-4));
    CHECK_THAT(mean.get<2>(), Catch::Matchers::WithinRel(23.083140961377627, 1.e-4));

    Vector<double> variance_vector{0, 0, 0};
    for (const auto& value: generated_values) {
        const auto diff_vector = value - mean;
        const auto x = diff_vector.get<0>();
        const auto y = diff_vector.get<1>();
        const auto z = diff_vector.get<2>();
        const Vector<double> squared{x * x, y * y, z * z};
        variance_vector += squared;
    }
    variance_vector = variance_vector / static_cast<double>(nsamples);

    const auto std_x = std::sqrt(variance_vector.get<0>());
    const auto std_y = std::sqrt(variance_vector.get<1>());
    const auto std_z = std::sqrt(variance_vector.get<2>());

    CHECK_THAT(std_x, Catch::Matchers::WithinRel(24.216065077138651, 1.e-4));
    CHECK_THAT(std_y, Catch::Matchers::WithinRel(23.710892361597871, 1.e-4));
    CHECK_THAT(std_z, Catch::Matchers::WithinRel(23.85119477289448, 1.e-4));

    Vector<double> variance_phys{0, 0, 0};
    for (const auto& value: generated_values) {
        const auto x = value.get<0>();
        const auto y = value.get<1>();
        const auto z = value.get<2>();
        const Vector<double> squared{x * x, y * y, z * z};
        variance_phys += squared;
    }
    variance_phys = variance_phys / static_cast<double>(nsamples);

    const auto std_x_phys = std::sqrt(variance_phys.get<0>());
    const auto std_y_phys = std::sqrt(variance_phys.get<1>());
    const auto std_z_phys = std::sqrt(variance_phys.get<2>());

    CHECK_THAT(std_x_phys, Catch::Matchers::WithinRel(30.332015370196821, 1.e-4));
    CHECK_THAT(std_y_phys, Catch::Matchers::WithinRel(24.163684840356566, 1.e-4));
    CHECK_THAT(std_z_phys, Catch::Matchers::WithinRel(33.192030500368219, 1.e-4));
}