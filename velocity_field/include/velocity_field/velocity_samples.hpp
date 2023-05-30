#ifndef STG_VELOCITY_SAMPLES_HPP
#define STG_VELOCITY_SAMPLES_HPP

#include "velocity_field.hpp"

/*
 * [
 *  Field{
 *    vector vx
 *    vector vy
 *    vector vz
 *  },
 *  Field{
 *    vector vx
 *    vector vy
 *    vector vz
 *  }
 * ]
 *
 * Can take views on { Field_1{vx[ivert], Field_2{vx[ivert], ... }
 */

namespace stg::field {

    template<std::floating_point T>
    class VelocitySamples final {
    public:
        using value_type = T;

        VelocitySamples(std::size_t samples_amount, std::size_t vertices_amount)
            : velocity_samples_(samples_amount) {
            for (auto& field: velocity_samples_) {
                field.resize(vertices_amount);
            }
        }

        VelocitySamples(std::vector<VelocityField<T>> samples)
            : velocity_samples_{std::move(samples)} {}

        auto vx_component_for_vertex(std::size_t ivert) const {
            return velocity_samples_ | ranges::views::transform(
                                               [ivert](const VelocityField<value_type>& field) {
                                                   return field.vx_view()[ivert] /* | ranges::views::drop(ivert) | ranges::views::take(1)*/;
                                               }) /* | ranges::views::join*/;
        }

        auto vy_component_for_vertex(std::size_t ivert) const {
            return velocity_samples_ | ranges::views::transform(
                                               [ivert](const VelocityField<value_type>& field) {
                                                   return field.vy_view()[ivert] /* | ranges::views::drop(ivert) | ranges::views::take(1)*/;
                                               }) /* | ranges::views::join*/;
        }

        auto vz_component_for_vertex(std::size_t ivert) const {
            return velocity_samples_ | ranges::views::transform(
                                               [ivert](const VelocityField<value_type>& field) {
                                                   return field.vz_view()[ivert] /* | ranges::views::drop(ivert) | ranges::views::take(1)*/;
                                               }) /* | ranges::views::join*/;
        }

        auto velocity_for_vertex(std::size_t ivert) const {
            return velocity_samples_ | ranges::views::transform(
                                               [ivert](const VelocityField<value_type>& field) {
                                                   return field.values_view()[ivert];
                                               });
        }

        std::size_t size() const { return velocity_samples_.size(); }

        void resize(std::size_t new_size) { velocity_samples_.resize(new_size); }

        const VelocityField<value_type>& sample(std::size_t isample) const {
            return velocity_samples_[isample];
        }

        VelocityField<value_type>& sample(std::size_t isample) {
            return velocity_samples_[isample];
        }

        void set_sample(VelocityField<value_type> sample, std::size_t isample) {
            velocity_samples_[isample] = std::move(sample);
        }

    private:
        std::vector<VelocityField<value_type>> velocity_samples_;
    };
}// namespace stg::field

#endif//STG_VELOCITY_SAMPLES_HPP
