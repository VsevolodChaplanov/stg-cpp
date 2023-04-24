#ifndef STG_STANDARD_DEVIATION_HPP
#define STG_STANDARD_DEVIATION_HPP

#include <execution>
#include <ranges>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <future>
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
    static auto std(FirstRange f_range, SecondRange s_range, MeanValueType f_mean, MeanValueType s_mean) {
      const std::size_t f_size = ranges::distance(f_range);
      const std::size_t s_size = ranges::distance(s_range);
      if (f_size != s_size) throw std::invalid_argument("Ranges have different lengths");

      MeanValueType sum = 0.;
      for (const auto [f, s] : rv::zip(f_range, s_range)) {
        sum += (f - f_mean) * (s - s_mean);
      };

      return std::sqrt(sum / f_size);
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange>
    static auto std(FirstRange f_range, SecondRange s_range) {
      auto f_mean = std::async(std::launch::async, [f_range] {return Mean::mean(f_range); });
      auto s_mean = std::async(std::launch::async, [s_range] {return Mean::mean(s_range); });
      return std(std::forward<FirstRange>(f_range),
                 std::forward<SecondRange>(s_range),
                 f_mean.get(), s_mean.get());
    }
  };
}

#endif //STG_STANDARD_DEVIATION_HPP
