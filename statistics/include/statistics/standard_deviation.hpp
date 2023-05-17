#ifndef STG_STANDARD_DEVIATION_HPP
#define STG_STANDARD_DEVIATION_HPP

#include <execution>
#include <ranges>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <future>
#include <fmt/format.h>
#include "concepts.hpp"
#include "mean.hpp"


namespace stg::statistics {
namespace rv = ranges::views;

  class StandardDeviation {
    template<typename Iter>
    using IteratorValueType = typename std::iterator_traits<Iter>::value_type;

  public:
    template<NumericViewable Range>
    static auto std(Range&& range, RangeValueType<Range> mean) {
      const std::size_t size = ranges::distance(range);
      const auto xx_product = ranges::inner_product(range, range, 0.);
      const auto xx_mean = xx_product / size;
      auto result = std::sqrt(xx_mean - mean * mean);
      return result;
    }

    template<NumericViewable Range>
    static auto std(Range&& range) {
      const auto mean = Mean::mean(range);
      return std(std::forward<Range>(range), mean);
    }

    template<std::forward_iterator Iter>
    static auto std(Iter begin, Iter end) {
      const auto mean = Mean::mean(begin, end);
      const auto result = std(begin, end, mean);
      return result;
    }

    template<std::forward_iterator Iter>
    static auto std(Iter begin, Iter end, IteratorValueType<Iter> mean) {
      const size_t size = std::distance(begin, end);
      const auto sqr_sum = std::inner_product(begin, end, begin, 0.);
      const auto result = std::sqrt(sqr_sum / size - mean * mean);
      return result;
    }


    template<NumericViewable FirstRange, NumericViewable SecondRange, std::floating_point MeanValueType>
    static auto std(FirstRange&& f_range, SecondRange&& s_range, MeanValueType f_mean, MeanValueType s_mean) {
      const std::size_t f_size = ranges::distance(f_range);
      const std::size_t s_size = ranges::distance(s_range);
      if (f_size != s_size) throw std::invalid_argument("Ranges have different lengths");

      MeanValueType sum = 0.;
      for (const auto [f, s] : rv::zip(f_range, s_range)) {
        sum += (f - f_mean) * (s - s_mean);
      }

      return std::sqrt(sum / f_size);
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange>
    static auto std(FirstRange&& f_range, SecondRange&& s_range) {
      auto f_mean = std::async(std::launch::async, [f_range] {return Mean::mean(f_range); });
      auto s_mean = std::async(std::launch::async, [s_range] {return Mean::mean(s_range); });
      return std(std::forward<FirstRange>(f_range),
                 std::forward<SecondRange>(s_range),
                 f_mean.get(), s_mean.get());
    }
  };

  namespace detail {
    class Integrator final {
    public:
      template<NumericViewable Function, NumericViewable Points>
      static auto integrate_trapezoidal_method(Function&& function, Points&& x)  {
        using function_value_t = typename ranges::range_value_t<Function>;
        using point_value_t = typename ranges::range_value_t<Function>;
        using common_type = std::decay_t<typename std::common_type<function_value_t, point_value_t>::type>;
        static_assert(std::is_convertible_v<function_value_t, point_value_t>
                      && std::is_convertible_v<point_value_t, function_value_t>,
                      "Values under ranges are not convertible to each other");
        if (ranges::distance(function) != ranges::distance(x)) {
          throw std::runtime_error("Distances of ranges are not equal");
        }
        common_type result{0};
        for (const auto index : rv::iota(0ull, ranges::distance(function) - 1ull)) {
          const auto dx = x[index + 1] - x[index];
          const auto funcs_sum = function[index + 1] + function[index];
          result += dx * funcs_sum / 2;
        }

        return result;
      }
    };
  }

  class StdIntegralAveraged final {
  public:

    /*
     * Vert-based algorithm
     */
    template<NumericViewable Function, std::floating_point TimeDelta = ranges::range_value_t<Function>>
    static auto std_integral_averaged(Function&& function, TimeDelta dtime, TimeDelta time_period) {
      using function_value_t = typename ranges::range_value_t<Function>;
      static_assert(std::is_convertible_v<function_value_t, TimeDelta>
                    && std::is_convertible_v<TimeDelta, function_value_t>,
                    "Values under ranges are not convertible to each other");
      ranges::viewable_range auto function_squared = function | rv::transform([](auto value) { return value * value; });
      auto points_view = rv::iota(0ull, static_cast<std::size_t>(ranges::distance(function)))
          | rv::transform([dtime](auto val) { return val * dtime; });
      auto result = detail::Integrator::integrate_trapezoidal_method(function_squared, points_view) / time_period;
      return result;
    }

    /*
     * Vert-based algorithm
     */
    template<NumericViewable Function,
             NumericViewable TimePoints,
             std::floating_point MeanType = ranges::range_value_t<Function>,
             std::floating_point TimeType = ranges::range_value_t<TimePoints>>
    static auto std_integral_averaged(Function&& function, MeanType mean, TimePoints time_points, TimeType period) {
      auto function_squared = function | rv::transform([mean](auto value) { return (value - mean) * (value - mean); });
      auto result = detail::Integrator::integrate_trapezoidal_method(
          function_squared, std::forward<TimePoints>(time_points)
      );
      return result / period;
    }
  };
}

#endif //STG_STANDARD_DEVIATION_HPP
