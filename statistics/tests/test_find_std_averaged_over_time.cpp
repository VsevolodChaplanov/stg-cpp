#include "common.hpp"
#include <numbers>

struct AveragedOverTimeStdFixture {
    const std::size_t points_amount = 1000;
    const double period = std::numbers::pi / 2;
};

SCENARIO_METHOD(AveragedOverTimeStdFixture, "Test average over time algorithms") {
    {
        auto function = [](double x) { return 1.; };
        const double dh = period / (points_amount - 1);
        auto function_view = rv::iota(0ull, points_amount) | rv::transform([dh](auto val) { return val * dh; }) | rv::transform([function](auto point) { return function(point); });
        const auto result = StdIntegralAveraged::std_integral_averaged(function_view, dh, period);
        CHECK_THAT(result, WithinRel(1., 1.e-4));
    }

    {
        auto function = [](double x) {
            return std::sqrt(x);
        };
        const double dh = period / (points_amount - 1);
        auto function_view = rv::iota(0ull, points_amount) | rv::transform([dh](auto val) { return val * dh; }) | rv::transform([function](auto point) { return function(point); });
        const auto result = StdIntegralAveraged::std_integral_averaged(function_view, dh, period);
        CHECK_THAT(result, WithinRel(0.78540, 1.e-4));
    }

    {
        auto function = [](double x) {
            static double epi = std::pow(std::numbers::e, std::numbers::pi);
            return 5 / (epi - 2) * std::pow(std::numbers::e, 2 * x) * std::cos(x);
        };
        const double dh = period / (points_amount - 1);
        auto function_view = rv::iota(0ull, points_amount) | rv::transform([dh](auto val) { return val * dh; }) | rv::transform([function](auto point) { return function(point); });
        const auto result = StdIntegralAveraged::std_integral_averaged(function_view, dh, period);
        CHECK_THAT(result, WithinRel(0.46872, 1.e-4));
    }

    {
        auto function = [](double x) {
            return std::pow(std::numbers::e, -2 * x) * std::cos(18 * std::numbers::pi * x);
        };
        const double dh = period / (points_amount - 1);
        auto function_view = rv::iota(0ull, points_amount) | rv::transform([dh](auto val) { return val * dh; }) | rv::transform([function](auto point) { return function(point); });
        const auto result = StdIntegralAveraged::std_integral_averaged(function_view, dh, period);
        CHECK_THAT(result, WithinRel(0.079533, 1.e-4));
    }

    {
        auto function = [](double x) { return x; };
        const double dh = period / (points_amount - 1);
        auto function_view = rv::iota(0ull, points_amount) | rv::transform([dh](auto val) { return val * dh; }) | rv::transform([function](auto point) { return function(point); });
        const auto result = StdIntegralAveraged::std_integral_averaged(function_view, dh, period);
        CHECK_THAT(result, WithinRel(0.82247, 1.e-4));
    }
}

SCENARIO_METHOD(AveragedOverTimeStdFixture, "Test 1d integration by trapezoidal method routine") {
    const std::size_t points_amount = 1000;

    {
        auto function = [](double x) {
            static double epi = std::pow(std::numbers::e, std::numbers::pi);
            return 5 / (epi - 2) * std::pow(std::numbers::e, 2 * x) * std::cos(x);
        };
        const double dh = period / (points_amount - 1);
        std::vector<double> points;
        points.resize(points_amount);
        ranges::generate(points, [dh]() {
            static std::size_t init = 0;
            return dh * init++;
        });
        std::vector<double> values;
        std::transform(points.cbegin(), points.cend(), std::back_inserter(values), function);
        const auto result = detail::Integrator::integrate_trapezoidal_method(values | rv::all, points | rv::all);
        CHECK_THAT(result, WithinRel(1.0, 1.e-4));
    }

    {
        auto function = [](double x) {
            static double epi = std::pow(std::numbers::e, std::numbers::pi);
            return 5 / (epi - 2) * std::pow(std::numbers::e, 2 * x) * std::cos(x);
        };
        const double dh = period / (points_amount - 1);
        auto points_view = rv::iota(0ull, static_cast<std::size_t>(points_amount)) | rv::transform([dh](auto val) { return val * dh; });
        auto function_view = points_view | rv::transform([function](auto point) {
                                 return function(point);
                             });
        const auto result_by_views = detail::Integrator::integrate_trapezoidal_method(function_view, points_view);
        CHECK_THAT(result_by_views, WithinRel(1.0, 1.e-4));
        ;
    }

    {
        auto function = [](double x) {
            return 1.;
        };
        const double dh = period / (points_amount - 1);
        auto points_view = rv::iota(0ull, static_cast<std::size_t>(points_amount)) | rv::transform([dh](auto val) { return val * dh; });
        auto function_view = points_view | rv::transform([function](auto point) {
                                 return function(point);
                             });
        const auto result_by_views = detail::Integrator::integrate_trapezoidal_method(function_view, points_view);
        CHECK_THAT(result_by_views, WithinRel(std::numbers::pi / 2, 1.e-4));
        ;
    }

    {
        auto function = [](double x) {
            return 2.;
        };
        const double dh = period / (points_amount - 1);
        auto points_view = rv::iota(0ull, static_cast<std::size_t>(points_amount)) | rv::transform([dh](auto val) { return val * dh; });
        auto function_view = points_view | rv::transform([function](auto point) {
                                 return function(point);
                             });
        const auto result_by_views = detail::Integrator::integrate_trapezoidal_method(function_view, points_view);
        CHECK_THAT(result_by_views, WithinRel(std::numbers::pi, 1.e-4));
        ;
    }

    {
        auto function = [](double x) {
            return x;
        };
        const double dh = period / (points_amount - 1);
        auto points_view = rv::iota(0ull, static_cast<std::size_t>(points_amount)) | rv::transform([dh](auto val) { return val * dh; });
        auto function_view = points_view | rv::transform([function](auto point) {
                                 return function(point);
                             });
        const auto result_by_views = detail::Integrator::integrate_trapezoidal_method(function_view, points_view);
        CHECK_THAT(result_by_views, WithinRel(std::numbers::pi * std::numbers::pi / 8, 1.e-4));
        ;
    }
}