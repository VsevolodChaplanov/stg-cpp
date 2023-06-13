#ifndef STG_STG_GENERATORS_DIRECT_ENERDY_SPECTRAL_METHOD_HPP
#define STG_STG_GENERATORS_DIRECT_ENERDY_SPECTRAL_METHOD_HPP

#include "geometry/geometry.hpp"
#include "stg_generators/fluctuation_fourier_sum.hpp"
#include "stg_generators/i_spectral_generator.hpp"
#include "stg_random/rn_generator_impl.hpp"
#include <cmath>
#include <concepts>
#include <coroutine>
#include <cstddef>
#include <functional>
#include <limits>
#include <numbers>
#include <random>
#include <ranges>
#include <stg_rangom.hpp>
#include <vector>

namespace stg::generators {
    using namespace std::numbers;

    template<std::floating_point T>
    class DirectEnergySpectralGenerator final : public ISpectralGenerator<T> {
    public:
        using value_type = T;
        Vector<T> operator()(const Point<value_type>& space_point, value_type time_point) const override {
            return EvenFluctuationFouirerSum<value_type>::calculate_sum(q_amplitudes_,
                                                                        k_wavenumbers_,
                                                                        phases_,
                                                                        frequencies_,
                                                                        k_unit_vectors_,
                                                                        sigma_unit_vectors_,
                                                                        space_point, time_point);
        }

        DirectEnergySpectralGenerator(std::size_t n_fourier,
                                      value_type l, value_type dh,
                                      value_type omega_0,
                                      std::function<value_type(value_type)> energy_function,
                                      std::size_t seed = std::mt19937_64::default_seed)
            : DirectEnergySpectralGenerator{n_fourier, l, l, l, dh, dh, dh, omega_0, energy_function, seed} {}

        DirectEnergySpectralGenerator(
                std::size_t n_fourier,
                value_type lx, value_type ly, value_type lz,
                value_type dx, value_type dy, value_type dz,
                value_type omega_0,
                std::function<value_type(value_type)> energy_function,
                std::size_t seed = std::mt19937_64::default_seed)
            : seed_{seed},
              n_fourier_{n_fourier},
              lx_{lx}, ly_{ly}, lz_{lz}, dx_{dx}, dy_{dy}, dz_{dz},
              omega_0_{omega_0},
              energy_function_{std::move(energy_function)} {
            k_unit_vectors_.resize(n_fourier_);
            sigma_unit_vectors_.resize(n_fourier_);
            q_amplitudes_.resize(n_fourier_);
            k_wavenumbers_.resize(n_fourier_);
            frequencies_.resize(n_fourier_);
            phases_.resize(n_fourier_);

            for (const std::size_t m: std::views::iota(0ul, n_fourier_)) {
                const auto k_mod = k_m(m);
                const auto dk = dk_m(m);
                const auto wave_vector_unit = k_unit();
                const auto zeta_vector_unit = zeta_unit();
                k_unit_vectors_[m] = wave_vector_unit;
                sigma_unit_vectors_[m] = sigma_m(zeta_vector_unit, wave_vector_unit);
                q_amplitudes_[m] = fourier_amplitude(k_mod, dk);
                k_wavenumbers_[m] = k_mod;
                frequencies_[m] = omega_generator_();
                phases_[m] = psi_generator_();
            }
        }

        ~DirectEnergySpectralGenerator() override = default;

    private:
        std::size_t seed_;
        std::size_t n_fourier_;
        value_type lx_, ly_, lz_, dx_, dy_, dz_, omega_0_;
        std::function<value_type(value_type)> energy_function_;
        value_type k_0_ = std::max({2 * pi_v<value_type> / lx_, 2 * pi_v<value_type> / ly_, 2 * pi_v<value_type> / lz_});
        value_type k_max_ = std::max({pi_v<value_type> / dx_, pi_v<value_type> / dy_, pi_v<value_type> / dz_});
        RNGenerator<std::mt19937_64, std::uniform_real_distribution<>> theta_generator_{
                generator_engines::get_engine<std::mt19937_64>(seed_),
                std::uniform_real_distribution<value_type>{-1, 1}};
        RNGenerator<std::mt19937_64, std::uniform_real_distribution<>> phi_generator_{
                generator_engines::get_engine<std::mt19937_64>(seed_),
                std::uniform_real_distribution<value_type>{0, 2 * std::numbers::pi_v<value_type>}};
        RNGenerator<std::mt19937_64, std::uniform_real_distribution<>> psi_generator_{
                generator_engines::get_engine<std::mt19937_64>(seed_),
                std::uniform_real_distribution<value_type>{-std::numbers::pi_v<value_type> / 2., std::numbers::pi_v<value_type> / 2.}};
        RNGenerator<std::mt19937_64, std::normal_distribution<>> omega_generator_{
                generator_engines::get_engine<std::mt19937_64>(seed_),
                std::normal_distribution<>{0, omega_0_}};
        std::vector<Vector<value_type>> k_unit_vectors_;
        std::vector<Vector<value_type>> sigma_unit_vectors_;
        std::vector<value_type> q_amplitudes_;
        std::vector<value_type> k_wavenumbers_;
        std::vector<value_type> frequencies_;
        std::vector<value_type> phases_;

        Vector<value_type> sigma_m(const Vector<value_type>& zeta, const Vector<T>& k_unit) const {
            const auto zeta_c_k = cross_product(zeta, k_unit);
            const auto sigma_m = scale_to_length(zeta_c_k, 1.);
            return sigma_m;
        }

        Vector<value_type> k_unit() {
            const auto cos_theta = theta_generator_();
            const auto theta_m = std::acos(cos_theta);
            const auto phi_m = phi_generator_();
            const auto kx = std::sin(theta_m) * std::cos(phi_m);
            const auto ky = std::sin(theta_m) * std::sin(phi_m);
            const auto kz = std::cos(theta_m);
            return {kx, ky, kz};
        }

        Vector<value_type> zeta_unit() {
            return k_unit();
        }

        value_type fourier_amplitude(value_type k_m, value_type dk_m) const {
            const auto energy = energy_function_(k_m);
            const auto kin_energy = energy * dk_m;
            const auto amplitude = std::sqrt(kin_energy);
            return amplitude;
        }

        value_type k_m(std::size_t m) const {
            const value_type dk = dk_m(m);
            const auto k_m = k_0_ + dk * m;
            return k_m;
        }

        value_type dk_m(std::size_t m) const {
            return (k_max_ - k_0_) / static_cast<value_type>(n_fourier_);
        }
    };
}// namespace stg::generators

#endif