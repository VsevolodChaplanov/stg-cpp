#ifndef STG_FLUCTUATION_FOURIER_SUM_HPP
#define STG_FLUCTUATION_FOURIER_SUM_HPP

#include <algorithm>
#include <bits/ranges_base.h>
#include <concepts>
#include <numeric>
#include <range/v3/view/transform.hpp>
#include <ranges>

#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/zip.hpp>

#include <geometry/geometry.hpp>
#include <stdexcept>
#include <utility>


namespace stg::generators {

    /*  Generic formulae 
        \vec{u} (\vec{x}, t) = 
            \sum_{n=1}^{N} 
                \left( 
                    \vec{p} \cos{(\vec{k_n} \cdot \vec{x} + \omega_n t)} 
                + 
                    \vec{q} \sin{(\vec{k_n} \cdot \vec{x} + \omega_n t)} 
                \right) 
        Only do summuation without tranformations
    */
    template<std::floating_point T>
    class FluctuationFouirerSum final {
    public:
        using value_type = T;

        template<std::ranges::viewable_range P, std::ranges::viewable_range Q,
                 std::ranges::viewable_range K, std::ranges::viewable_range Omega>
        static Vector<value_type> calculate_sum(P&& p_vectors, Q&& q_vectors, K&& wave_vectors_k, Omega&& frequencies, Point<value_type> space_point, value_type time) {
            return calculate_sum(ranges::views::zip(p_vectors, q_vectors, wave_vectors_k, frequencies),
                                 std::move(space_point), time);
        }

        // zipped as p, q, k, omega
        template<std::ranges::viewable_range ValuesZip>
        static Vector<value_type> calculate_sum(ValuesZip&& ranges_zip, Point<value_type> space_point, value_type time) {
            Vector<value_type> result_fluctuation{value_type{0}, value_type{0}, value_type{0}};

            result_fluctuation = ranges::accumulate(std::forward<ValuesZip>(ranges_zip), result_fluctuation,
                                                    [point = std::move(space_point), time](const Vector<value_type>& accumulator, const auto& vectors_tuple) {
                                                        const auto [p, q, k, omega] = vectors_tuple;
                                                        const auto phase = dot_product(k, point) + omega * time;
                                                        const auto u_m = p * std::cos(phase) + q * std::sin(phase);
                                                        return accumulator + u_m;
                                                    });

            return result_fluctuation;
        }
    };

    template<std::floating_point T>
    class EvenFluctuationFouirerSum final {
    public:
        using value_type = T;

        template<std::ranges::viewable_range Q, std::ranges::viewable_range K, std::ranges::viewable_range Psi,
                 std::ranges::viewable_range Omega,
                 std::ranges::viewable_range KV, std::ranges::viewable_range Sigma>
        static Vector<value_type> calculate_sum(Q&& q, K&& k, Psi&& psi, Omega&& omega, KV&& k_units, Sigma&& sigma_units,
                                                const Point<value_type>& space_point, value_type time) {
            auto zip = ranges::views::zip(std::forward<Q>(q),
                                          std::forward<K>(k),
                                          std::forward<Psi>(psi),
                                          std::forward<Omega>(omega),
                                          std::forward<KV>(k_units),
                                          std::forward<Sigma>(sigma_units));
            Vector<value_type> result{value_type{0}, value_type{0}, value_type{0}};

            result = ranges::accumulate(std::move(zip), result, [&space_point, time](const Vector<value_type>& accumulated, const auto& vectors_tuple) {
                const auto [q_m, k_m, psi_m, omega_m, wave_vector_unit, sigma_vector_unit] = vectors_tuple;
                const auto phase = dot_product(wave_vector_unit, space_point) * k_m + psi_m + omega_m * time;
                const auto amplitude_magnitude = 2 * q_m * std::cos(phase);
                const auto u_m = sigma_vector_unit * amplitude_magnitude;
                return accumulated + u_m;
            });

            return result;
        };
    };
}// namespace stg::generators

#endif