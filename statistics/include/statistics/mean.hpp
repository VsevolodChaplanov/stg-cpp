#ifndef STG_MEAN_HPP
#define STG_MEAN_HPP

#include <execution>
#include <ranges>
#include <algorithm>
#include <numeric>
#include "concepts.hpp"

namespace stg::statistics {
  /*
   * Calculate statistics of given range
   * Has concepts overloading for
   * - range
   * - iterator pair
   */
  class Mean {
  public:
    template<AdditiveRange Range>
    static auto mean(Range&& range) {
      const std::size_t size = ranges::distance(range);
      const auto accum = ranges::accumulate(range, 0.);
      return accum / size;
    }

    /*
     * Requires forward iterator due using par execution policy
     */
    template<AdditiveIterator Iter>
    static auto mean(Iter begin, Iter end) {
      const size_t size = std::distance(begin, end);
      const auto accumulate = std::reduce(std::execution::par_unseq,
                                          begin, end, 0., std::plus{});
      const auto result = accumulate / size;
      return result;
    }
  };
}

#endif //STG_MEAN_HPP
