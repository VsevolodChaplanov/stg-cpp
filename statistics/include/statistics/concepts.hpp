#ifndef STG_CONCEPTS_STATISTICS_HPP
#define STG_CONCEPTS_STATISTICS_HPP

#include <concepts>
#include <iterator>
#include <range/v3/all.hpp>

namespace stg::statistics {
  template<typename Iter>
  using IterValueType = typename std::iterator_traits<Iter>::value_type;

  template<typename Range>
  using RangeValueType = typename ranges::range_value_t<Range>;

  template<typename T>
  concept Additive = requires (T a, T b) {
    { a + b } -> std::convertible_to<T>;
  };

  template<typename T>
  concept Deductible = requires (T a, T b) {
    { a - b } -> std::convertible_to<T>;
  };

  template<typename T>
  concept Multiplicable = requires (T a, T b) {
    { a * b } -> std::convertible_to<T>;
  };

  template<typename T, typename U>
  concept Divideable = requires (T a, U b) {
    { a / b } -> std::convertible_to<T>;
  };

  template<typename Range>
  concept AdditiveRange = ranges::viewable_range<Range> &&
    Additive<RangeValueType<Range>> &&
    Divideable<RangeValueType<Range>, std::size_t>;

  template<typename Iter>
  concept AdditiveIterator = std::forward_iterator<Iter> &&
    Additive<typename std::iterator_traits<Iter>::value_type> &&
    Divideable<typename std::iterator_traits<Iter>::value_type, std::size_t>;

  template<typename Range>
  concept NumericViewable = ranges::viewable_range<Range>
                            && Additive<RangeValueType<Range>>
                            && Divideable<RangeValueType<Range>, std::size_t>
                            && Deductible<RangeValueType<Range>>
                            && Multiplicable<RangeValueType<Range>>
                            && std::is_arithmetic_v<RangeValueType<Range>>;
}

#endif //STG_CONCEPTS_HPP
