#ifndef STG_SPECTRAS_BASE_HPP
#define STG_SPECTRAS_BASE_HPP

#include "detail.hpp"
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>

namespace stg::generators {

    template<std::floating_point T>
    class ISpectra {
    public:
        using value_type = T;

        virtual value_type operator()(value_type k_x, value_type k_y, value_type k_z) const noexcept = 0;
        virtual value_type operator()(value_type k_mod) const noexcept = 0;
        virtual ~ISpectra() = default;
    };

    template<std::floating_point T>
    class DeltaSpectra final : public ISpectra<T> {
    public:
        using value_type = T;

        explicit DeltaSpectra(value_type v_0, value_type k_0) noexcept
            : v_0_{v_0}, k_0_{k_0} {}

        value_type operator()(value_type k_x, value_type k_y, value_type k_z) const noexcept override {
            const auto length = std::sqrt(k_x * k_x + k_y * k_y + k_z * k_z);
            return std::fabs(length - k_0_) <= std::numeric_limits<value_type>::epsilon() ? 3. / 2. * v_0_ * v_0_ : 0.;
        }

        value_type operator()(value_type k_mod) const noexcept override {
            return std::fabs(k_mod - k_0_) <= std::numeric_limits<value_type>::epsilon() ? 3. / 2. * v_0_ * v_0_ : 0.;
        }

        ~DeltaSpectra() override = default;

    private:
        const value_type v_0_;
        const value_type k_0_;
    };

    template<std::floating_point T>
    class GaussianSpectra final : public ISpectra<T> {
    public:
        using value_type = T;

        value_type operator()(value_type k_x, value_type k_y, value_type k_z) const noexcept override;
        value_type operator()(value_type k_mod) const noexcept override;
        ~GaussianSpectra() override = default;
    };
}// namespace stg::generators


#endif