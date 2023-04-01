#ifndef STG_IVELOCITY_FIELD_HPP
#define STG_IVELOCITY_FIELD_HPP

#include <vector>
#include <concepts>
#include <range/v3/range.hpp>
#include <range/v3/view.hpp>
#include <geometry/geometry.hpp>
#include "concepts.hpp"

namespace stg::field {

  template<std::floating_point T>
  class IVelocityField {
  public:
    virtual std::size_t size() const = 0;

    virtual ~IVelocityField() = default;
  };
}

#endif //STG_IVELOCITY_FIELD_HPP
