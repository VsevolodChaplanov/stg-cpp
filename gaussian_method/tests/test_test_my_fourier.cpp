#include "common.hpp"
#include "ft.hpp"
#include "rtable/cube_relation_table.hpp"
#include "space.hpp"
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <cstddef>
#include <geometry/geometry.hpp>
#include <iostream>
#include <numbers>
#include <ranges>
#include <vector>


TEST_CASE("Calculate own fourier transformations") {
    const std::size_t n = 21;
    const double edge_length = 3.;
    std::vector<double> values(n * n * n, 1.);
    CubeMeshBuilder<double> builder{edge_length, 21};
    const auto fe_mesh_ = builder.build();
    const double dh = fe_mesh_->relation_table()->h();
    //     double hk = 2 * M_PI / L;
    // double Lk = hk * N;
    // return {N, Lk};
    const double dk = 2 * std::numbers::pi / edge_length;
    CubeMeshBuilder<double> fourier_mesh_builder{dk * n, n};
    const auto fourier_fe_mesh_ = fourier_mesh_builder.build();

    CHECK(fe_mesh_->n_vertices() == fourier_fe_mesh_->n_vertices());

    const std::size_t nvert = fe_mesh_->n_vertices();
    std::vector<std::complex<double>> direct_fft_values(nvert, 0.);

    for (const std::size_t ivert: std::views::iota(0ul, nvert)) {
        const auto wave_vector = fourier_fe_mesh_->relation_table()->vertex(ivert);
        const auto value = fe_mesh_->integrate_fourier2(values, wave_vector);
        direct_fft_values[ivert] = value;
    }


    std::vector<std::complex<double>> values_to_check(nvert, 0.);
    for (const std::size_t ivert: std::views::iota(0ul, nvert)) {
        const auto point = fe_mesh_->relation_table()->vertex(ivert);
        const auto value = fourier_fe_mesh_->integrate_inverse_fourier2(direct_fft_values, point);
        values_to_check[ivert] = value;
    }

    std::cout << std::endl;
}

TEST_CASE("Calculate own fourier transformations based on own functions") {
    const std::size_t n = 21;
    const double edge_length = 3.;
    CubeMeshBuilder<double> builder{edge_length, 21};
    const auto fe_mesh_ = builder.build();
    const double dh = fe_mesh_->relation_table()->h();
    const double dk = 2 * std::numbers::pi / edge_length;
    const std::size_t nvert = fe_mesh_->n_vertices();
    std::vector<double> values(nvert, 1.);
    CubeMeshBuilder<double> fourier_mesh_builder{dk * n, n};
    const auto fourier_fe_mesh_ = fourier_mesh_builder.build();

    CHECK(fe_mesh_->n_vertices() == fourier_fe_mesh_->n_vertices());

    for (const std::size_t index: std::views::iota(0ul, nvert)) {
        auto func = [](double x, double y, double z) {
            double v1 = 1.0 / (std::abs(x) + 1);
            double v2 = 1.0 / (std::abs(y) + 1);
            double v3 = 1.0 / (std::abs(z) + 1);

            return std::pow(v1, 6) * std::pow(v2, 4) * std::pow(v3, 5);
        };
        const auto vertex = fe_mesh_->relation_table()->vertex(index);
        values[index] = func(vertex.get<0>(), vertex.get<1>(), vertex.get<2>());
    }

    const auto result_direct = fourier3(*fe_mesh_->relation_table(), values);
    const auto result_inverse = inverse_fourier3(*fourier_fe_mesh_->relation_table(), result_direct);

    std::cout << std::endl;
}