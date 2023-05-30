#ifndef STG_VELOCITY_FIELD_HPP
#define STG_VELOCITY_FIELD_HPP

#include "concepts.hpp"
#include "ivelocity_field.hpp"
#include <geometry/geometry.hpp>
#include <mutex>
#include <range/v3/all.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <vector>

namespace stg::field {
    namespace rv = ranges::views;

    template<std::floating_point T>
    class VelocityField final : IVelocityField<T> {
    public:
        using value_type = T;

        VelocityField() = default;

        explicit VelocityField(std::size_t nvert)
            : vx_(nvert), vy_(nvert), vz_(nvert) {}

        VelocityField(std::vector<value_type> vx,
                      std::vector<value_type> vy,
                      std::vector<value_type> vz) noexcept
            : vx_{std::move(vx)}, vy_{std::move(vy)}, vz_{std::move(vz)} {}

        explicit VelocityField(const std::vector<Vector<value_type>>& velocities)
            : vx_{velocities | ranges::views::transform([](const Vector<value_type>& vel) {
                      return vel.template get<0>();
                  }) |
                  ranges::to<std::vector<value_type>>()},
              vy_{velocities | ranges::views::transform([](const Vector<value_type>& vel) {
                      return vel.template get<1>();
                  }) |
                  ranges::to<std::vector<value_type>>()},
              vz_{velocities | ranges::views::transform([](const Vector<value_type>& vel) {
                      return vel.template get<2>();
                  }) |
                  ranges::to<std::vector<value_type>>()} {}

        VelocityField(const VelocityField& other) noexcept
            : vx_{other.vx_}, vy_{other.vy_}, vz_{other.vz_} {}

        VelocityField(VelocityField&& other) noexcept
            : vx_{std::move(other.vx_)}, vy_{std::move(other.vy_)}, vz_{std::move(other.vz_)} {}

        VelocityField& operator=(VelocityField&& other) noexcept {
            if (this == &other) { return *this; }
            vx_ = std::move(other.vx_);
            vy_ = std::move(other.vy_);
            vz_ = std::move(other.vz_);

            return *this;
        }

        [[nodiscard]] size_t size() const override { return vx_.size(); }

        ~VelocityField() override = default;

        auto vx_view() const { return ranges::views::all(vx_); }

        auto vy_view() const { return ranges::views::all(vx_); }

        auto vz_view() const { return ranges::views::all(vx_); }

        auto values_view() const { return ranges::views::zip(vx_, vy_, vz_); }

        Vector<value_type> value(std::size_t index) const {
            return {vx_[index], vy_[index], vz_[index]};
        }

        void set_value(Vector<value_type> value, std::size_t index) {
            std::lock_guard lg{mutex_};
            const value_type x = value.template get<0>();
            const value_type y = value.template get<1>();
            const value_type z = value.template get<2>();

            vx_[index] = x;
            vy_[index] = y;
            vz_[index] = z;
        }

        void set_value(value_type vx, value_type vy, value_type vz, std::size_t index) {
            std::lock_guard lg{mutex_};

            vx_[index] = vx;
            vy_[index] = vy;
            vz_[index] = vz;
        }

        [[nodiscard]] bool empty() const {
            return vx_.empty() && vy_.empty() && vz_.empty();
        }

        void resize(std::size_t new_size) {
            vx_.resize(new_size);
            vy_.resize(new_size);
            vz_.resize(new_size);
        }

    protected:
        std::vector<value_type> vx_;
        std::vector<value_type> vy_;
        std::vector<value_type> vz_;

        std::mutex mutex_;
    };
}// namespace stg::field

#endif//STG_VELOCITY_FIELD_HPP
