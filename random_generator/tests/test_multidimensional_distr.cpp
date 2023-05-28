#include "common.hpp"
#include "multidim_gaussian_genrator.hpp"
#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>
#include <ranges>

SCENARIO("Create random values using multidimensional distribution") {
    const double eps = 1.e-6;
    arma::vec3 mean{0, 0, 0};
    arma::mat33 covariance{1, 0, 0, 0, 1, 0, 0, 0, 1};

    stg::random::MultidimenasionalGaussianGenerator<double> dist{mean, covariance, 42};

    const auto generated_value = dist();

    CHECK_THAT(generated_value.get<0>(), WithinRel(-1.4436212261, eps));
    CHECK_THAT(generated_value.get<1>(), WithinRel(1.0039611555, eps));
    CHECK_THAT(generated_value.get<2>(), WithinRel(0.1236618641, eps));

    double sum = 0;
    for (const std::size_t index: std::views::iota(0ull, 1'000'000ull)) {
        const auto generated_vec = dist();
        sum += generated_vec.get<0>();
    }
    sum /= 1'000'000ull;

    CHECK_THAT(sum, WithinRel(-0.0004840747, eps));

    double std = 0;
    for (const std::size_t index: std::views::iota(0ull, 1'000'000ull)) {
        const auto generated_vec = dist();
        std += generated_vec.get<0>() * generated_vec.get<0>();
    }
    std /= 1'000'000ull;

    CHECK_THAT(std, WithinRel(1.0010585787, eps));

    std = 0;
    for (const std::size_t index: std::views::iota(0ull, 1'000'000ull)) {
        const auto generated_vec = dist();
        std += generated_vec.get<0>() * generated_vec.get<1>();
    }
    std /= 1'000'000ull;

    CHECK_THAT(std, WithinRel(-0.000316637, eps));
}