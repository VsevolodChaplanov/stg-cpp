#ifndef STG_VELOCITY_SAMPLES_1D_HPP
#define STG_VELOCITY_SAMPLES_1D_HPP

#include "velocity_field_1d.hpp"

namespace stg::field {

  template<std::floating_point T>
  class VelocitySamples1D final {
  public:
    using value_type = T;

    VelocitySamples1D(std::size_t samples_amount, std::size_t vertices_amount)
      : samples_{samples_amount, VelocitySamples1D{vertices_amount}} { }

    VelocitySamples1D(std::vector<VelocityField1D<value_type>>&& samples)
      : samples_{std::move(samples)} { }

    auto vx_component_for_vertex(std::size_t ivert) const {
      return samples_ | rv::transform([ivert] (const auto& field) {
        return field.vx_view()[ivert];
      });
    }

    void set_sample(VelocityField1D<value_type> field, std::size_t index) {
      samples_[index] = std::move(field);
    }

    const VelocityField1D<value_type>& sample(std::size_t index) const {
      return samples_[index];
    }

    VelocityField1D<value_type>& sample(std::size_t index) {
      return samples_[index];
    }
    field = {const stg::field::VelocityField1D<double> &}
    constexpr std::size_t size() const { return samples_.size(); }

  private:
    std::vector<VelocityField1D<value_type>> samples_;
  };
}

#endif //STG_VELOCITY_SAMPLES_1D_HPP
