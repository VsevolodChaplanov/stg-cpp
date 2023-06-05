#include "common.hpp"
#include "ft.hpp"
#include "mesh_builders/mesh_builders.hpp"
#include "rtable/cube_relation_table.hpp"
#include "space.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
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
    // const double dh = fe_mesh_->relation_table()->h();
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
    // const double dh = fe_mesh_->relation_table()->h();
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

    for (const std::size_t index: std::views::iota(0ul, nvert)) {
        CHECK_THAT(result_inverse[index], WithinRel(values[index], 1.e-4));
    }

    std::cout << std::endl;
}

TEST_CASE("Test base fourier routine") {
    auto func = [](double x, double y, double z) {
        double v1 = 1.0 / (std::abs(x) + 1);
        double v2 = 1.0 / (std::abs(y) + 1);
        double v3 = 1.0 / (std::abs(z) + 1);

        return std::pow(v1, 6) * std::pow(v2, 4) * std::pow(v3, 5);
    };
    double eps = 1.e-4;

    // 1d partition
    size_t N = 11;

    // physical space boundaries
    double Lx = 2;

    PhysicalSpace ps(N, Lx);
    FourierSpace fs = ps.fourier_space();

    // fill 3d function array
    std::vector<double> fx(N * N * N);

    for (size_t k = 0; k < N; ++k)
        for (size_t j = 0; j < N; ++j)
            for (size_t i = 0; i < N; ++i) {
                fx[ps.lin_index(i, j, k)] = func(ps.coo[i], ps.coo[j], ps.coo[k]);
            }

    // compute forward transform
    std::vector<double> fk = fourier3(ps, fx);
    CHECK_THAT(fk[5 + 5 * N + 5 * N * N], WithinRel(0.0005079819819736192, eps));
    CHECK_THAT(fk[6 + 5 * N + 5 * N * N], WithinRel(0.00036019289899672674, eps));
    CHECK_THAT(fk[7 + 5 * N + 5 * N * N], WithinRel(0.000216620583015753, eps));
    CHECK_THAT(fk[0], WithinRel(1.7043779074776304e-06, eps));
    CHECK_THAT(fk[4 + 4 * N + 4 * N * N], WithinRel(0.00011965752189560347, eps));
    CHECK_THAT(fk[4 + 5 * N + 7 * N * N], WithinRel(0.0001205415104388087, eps));

    // compute inverse transform
    std::vector<double> fx2 = inverse_fourier3(fs, fk);
    for (size_t i = 0; i < N * N * N; ++i) {
        CHECK_THAT(fx2[i], WithinRel(fx[i], eps));
    }

    std::cout << std::endl;
}

TEST_CASE("Test adopted fourier routine") {
    auto func = [](double x, double y, double z) {
        double v1 = 1.0 / (std::abs(x) + 1);
        double v2 = 1.0 / (std::abs(y) + 1);
        double v3 = 1.0 / (std::abs(z) + 1);

        return std::pow(v1, 6) * std::pow(v2, 4) * std::pow(v3, 5);
    };
    double eps = 1.e-4;

    // 1d partition
    size_t N = 11;

    // physical space boundaries
    double Lx = 2;

    CubeMeshBuilder<double> builder{Lx, N};

    const auto phys_mesh = builder.build();
    const auto fourier_mesh = phys_mesh->relation_table()->make_fourier_space<CubeMeshBuilder<double>>();

    // fill 3d function array
    std::vector<double> fx(N * N * N);

    for (size_t k = 0; k < N; ++k)
        for (size_t j = 0; j < N; ++j)
            for (size_t i = 0; i < N; ++i) {
                const auto lin_index = phys_mesh->relation_table()->lin_index(i, j, k);
                const auto point = phys_mesh->relation_table()->vertex(i, j, k);
                fx[lin_index] = func(point.get<0>(), point.get<1>(), point.get<2>());
            }

    // compute forward transform
    std::vector<double> fk = fourier3(*phys_mesh->relation_table(), fx);
    CHECK_THAT(fk[5 + 5 * N + 5 * N * N], WithinRel(0.0005079819819736192, eps));
    CHECK_THAT(fk[6 + 5 * N + 5 * N * N], WithinRel(0.00036019289899672674, eps));
    CHECK_THAT(fk[7 + 5 * N + 5 * N * N], WithinRel(0.000216620583015753, eps));
    CHECK_THAT(fk[0], WithinRel(1.7043779074776304e-06, eps));
    CHECK_THAT(fk[4 + 4 * N + 4 * N * N], WithinRel(0.00011965752189560347, eps));
    CHECK_THAT(fk[4 + 5 * N + 7 * N * N], WithinRel(0.0001205415104388087, eps));

    // compute inverse transform
    std::vector<double> fx2 = inverse_fourier3(*fourier_mesh, fk);
    for (size_t i = 0; i < N * N * N; ++i) {
        CHECK_THAT(fx2[i], WithinRel(fx[i], eps));
    }

    std::cout << std::endl;
}