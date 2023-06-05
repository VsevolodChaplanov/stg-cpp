#include "common.hpp"
#include "geometry/geometry.hpp"
#include "rtable/cube_relation_table.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <ranges>

struct IntegrationFixture {
    const double eps = 1.e-4;
    const std::size_t n = 100;
    const double edge = 3.;
    CubeMeshBuilder<double> builder{edge, n};
};


SCENARIO_METHOD(IntegrationFixture, "Integrate function over mesh") {
    GIVEN("Cube FE mesh") {
        const auto fe_mesh_ptr_ = builder.build();
        const auto nvert = fe_mesh_ptr_->n_vertices();

        THEN("Integrate constant function and get volume of domain") {
            const auto func_value = 1.5;
            const double exact_volume = func_value * edge * edge * edge;
            auto values = std::views::iota(0ul, nvert) | std::views::transform([func_value](auto index) {
                              return func_value;
                          });

            const auto result = fe_mesh_ptr_->integrate(values);
            CHECK_THAT(exact_volume, WithinRel(result, eps));
        }

        THEN("Integrate linear function") {
            auto func = [](const Point<double>& point) {
                const auto x = point.get<0>();
                const auto y = point.get<1>();
                const auto z = point.get<2>();
                return x + y + z;
            };
            const auto function_values = std::views::iota(0ul, nvert) | std::views::transform([&](auto index) {
                                             const auto point = fe_mesh_ptr_->relation_table()->vertex(index);
                                             return func(point);
                                         });

            const double exact_integratin_result = 0.;
            const auto result = fe_mesh_ptr_->integrate(function_values);

            CHECK_THAT(exact_integratin_result, WithinAbs(result, eps));
        }

        THEN("Integrate squared function") {
            auto func = [](const Point<double>& point) {
                const auto x = point.get<0>();
                const auto y = point.get<1>();
                const auto z = point.get<2>();
                return x * x + y * y + z * z;
            };
            const auto function_values = std::views::iota(0ul, nvert) | std::views::transform([&](auto index) {
                                             const auto point = fe_mesh_ptr_->relation_table()->vertex(index);
                                             return func(point);
                                         });

            const double exact_integratin_result = 243. / 4.;
            const auto result = fe_mesh_ptr_->integrate(function_values);

            CHECK_THAT(exact_integratin_result, WithinRel(result, eps));
        }
    }
}