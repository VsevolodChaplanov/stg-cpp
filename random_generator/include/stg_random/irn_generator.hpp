#ifndef STG_IRN_GENERATOR_HPP
#define STG_IRN_GENERATOR_HPP

#include <concepts>
#include <type_traits>
#include "generator_concept.hpp"

namespace stg {
  template<std::floating_point T>
  class IRNGenerator {
  public:
      using result_type = T;
      using engine_type = std::nullptr_t;
      using distribution = std::nullptr_t;
      virtual result_type operator()() = 0;
      virtual ~IRNGenerator() = default;
  };
}

#endif