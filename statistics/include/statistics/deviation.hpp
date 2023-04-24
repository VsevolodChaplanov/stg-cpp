#ifndef STG_DEVIATION_HPP
#define STG_DEVIATION_HPP

#include <algorithm>
#include <execution>
#include "concepts.hpp"

namespace stg::statistics {
namespace rv = ranges::views;

  class Deviation final {
  public:

    template<NumericViewable FirstRange, NumericViewable SecondRange>
    inline static NumericViewable auto deviation_range(FirstRange&& f_range, SecondRange&& s_range) {
      return rv::zip(std::forward<FirstRange>(f_range), std::forward<SecondRange>(s_range))
             | rv::transform([] (const auto& tuple_pair) {
                 return std::get<0>(tuple_pair) - std::get<1>(tuple_pair);
               });
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange>
    inline static NumericViewable auto deviation_range_abs(FirstRange&& f_range, SecondRange&& s_range) {
      return deviation_range(std::forward<FirstRange>(f_range), std::forward<SecondRange>(s_range))
             | rv::transform([] (const auto& value) { return std::fabs(value); });
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange>
    inline static NumericViewable auto deviation_range_sqr(FirstRange&& f_range, SecondRange&& s_range) {
      return deviation_range(std::forward<FirstRange>(f_range), std::forward<SecondRange>(s_range))
             | rv::transform([] (const auto& value) { return value * value; });
    }
  };
}

#endif //STG_DEVIATION_HPP
