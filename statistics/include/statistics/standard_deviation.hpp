#ifndef STG_STANDARD_DEVIATION_HPP
#define STG_STANDARD_DEVIATION_HPP

#include <execution>
#include <ranges>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "concepts.hpp"
#include "mean.hpp"

namespace stg::statistics {

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
  private:
  };
}

#endif //STG_STANDARD_DEVIATION_HPP
