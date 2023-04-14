#ifndef STG_VELOCITY_FIELD_1D_HPP
#define STG_VELOCITY_FIELD_1D_HPP

#include "concepts.hpp"
#include "ivelocity_field.hpp"

namespace stg::field {
namespace rv = ranges::views;

  template<std::floating_point T>
  class VelocityField1D final : IVelocityField<T> {
  public:
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    using value_type = T;

    VelocityField1D(std::size_t n) : vx_(n) { }
    VelocityField1D(std::vector<value_type> values) : vx_{std::move(values)} { }
    VelocityField1D(std::vector<value_type>&& values) : vx_{std::move(values)} { }

    template<std::input_iterator Iter>
    VelocityField1D(Iter begin, Iter end) : vx_{begin, end} { }

    iterator begin() { return vx_.begin(); }
    iterator end() { return vx_.end(); }
    const_iterator cbegin() { return vx_.cbegin(); }
    const_iterator cend() { return vx_.cend(); }

    auto vx_view() const { return vx_ | rv::all; }

    std::size_t size() const { return vx_.size(); }

  private:
    std::vector<value_type> vx_;
  };
}

#endif //STG_VELOCITY_FIELD_1D_HPP
