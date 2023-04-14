#ifndef STG_IAPPROXIMATOR_HPP
#define STG_IAPPROXIMATOR_HPP

#include <concepts>
#include <vector>
#include <functional>
#include "geometry/geometry.hpp"

namespace stg::mesh {

  template<std::floating_point T>
  class IApproximator {
  public:
    using value_type = T;
    virtual ~IApproximator() = default;

    virtual std::vector<value_type> approximate(std::function<double(Point<T>)> function) const = 0;
  };
}

#endif //STG_IAPPROXIMATOR_HPP
