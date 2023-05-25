#ifndef STG_FLUCTUATION_FOURIER_SUM_HPP
#define STG_FLUCTUATION_FOURIER_SUM_HPP

#include <algorithm>
#include <concepts>
#include <numeric>
#include <ranges>

#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/zip.hpp>

#include <geometry/geometry.hpp>
#include <utility>


namespace stg::generator {

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
            return calculate_sum(ranges::views::zip(std::forward<P>(p_vectors), std::forward<P>(q_vectors),
                                                    std::forward<P>(wave_vectors_k), std::forward<P>(frequencies)),
                                 std::move(space_point), time);
        }

        // zipped as p, q, k, omega
        template<std::ranges::viewable_range ValuesZip>
        static Vector<value_type> calculate_sum(ValuesZip&& ranges_zip, Point<value_type> space_point, value_type time) {
            Vector<value_type> result_fluctuation{value_type{0}, value_type{0}, value_type{0}};

            ranges::accumulate(std::forward<ValuesZip>(ranges_zip), result_fluctuation,
                               [point = std::move(space_point), time](const Vector<value_type>& accumulator, const auto& vectors_tuple) {
                                   const auto [p, q, k, omega] = vectors_tuple;
                                   const auto phase = cross_product(k, point) + omega * time;
                                   return accumulator + p * std::cos(phase) + q * std::sin(q);
                               });

            return result_fluctuation;
        }
    };
}// namespace stg::generator

#endif