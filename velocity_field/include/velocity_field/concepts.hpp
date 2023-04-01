#ifndef STG_CONCEPTS_HPP
#define STG_CONCEPTS_HPP

#include <concepts>
#include <ranges>
#include <range/v3/range_concepts.hpp>

namespace stg::field::concepts {

  template<typename T>
  concept ValueOrSampleConcept = std::is_floating_point_v<T> ||
                                 std::ranges::forward_range<T>;
}

#endif //STG_CONCEPTS_HPP
