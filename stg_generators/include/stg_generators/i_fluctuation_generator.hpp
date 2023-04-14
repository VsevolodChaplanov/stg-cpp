#ifndef STG_I_FLUCTUATION_GENERATOR_HPP
#define STG_I_FLUCTUATION_GENERATOR_HPP

#include <concepts>
#include <geometry/geometry.hpp>


namespace stg::generators {

  template<std::floating_point T>
  class IVelocityFluctuationGenerator {
  public:
    using value_type = T;

    virtual Vector<T> operator()(const Point<T>& real_space_point, T time) const = 0;
    virtual ~IVelocityFluctuationGenerator() = default;
  };
}

#endif //STG_I_FLUCTUATION_GENERATOR_HPP
