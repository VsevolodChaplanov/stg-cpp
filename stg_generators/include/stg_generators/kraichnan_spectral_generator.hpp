#ifndef STG_GENERATORS_KRAICHNAN_SPECTRAL_GENERATOR_HPP
#define STG_GENERATORS_KRAICHNAN_SPECTRAL_GENERATOR_HPP

#include "fluctuation_fourier_sum.hpp"
#include "geometry/geometry.hpp"
#include "i_spectral_generator.hpp"
#include "spectras_base.hpp"
#include "stg_coro_future/coroutine_future.hpp"
#include "stg_random/generator_engines.hpp"
#include "stg_random/irn_generator.hpp"
#include "stg_random/multidim_gaussian_genrator.hpp"
#include "stg_random/rn_generator_impl.hpp"
#include "stg_tensor/tensor.hpp"
#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <future>
#include <iterator>
#include <random>
#include <range/v3/range/conversion.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/generate.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/transform.hpp>
#include <ranges>
#include <stg_coro_future.hpp>
#include <stg_rangom.hpp>
#include <type_traits>
#include <utility>
#include <vector>

namespace stg::generators {

    template<std::floating_point T>
    class AKraichanGenerator : public ISpectralGenerator<T> {};

    /*
        Uses DeltaSpectra, generates velocity vectors with given amplitude
    */
    template<std::floating_point T>
    class KraichanGeneratorDeltaFunction final : public AKraichanGenerator<T> {
    public:
        using value_type = T;

        KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0, std::size_t seed = std::mt19937_64::default_seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value);
        KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0, std::vector<Vector<value_type>> k, std::size_t seed = std::mt19937_64::default_seed);
        KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0,
                                       std::vector<Vector<value_type>>&& v,
                                       std::vector<Vector<value_type>>&& w,
                                       std::vector<Vector<value_type>>&& k,
                                       std::vector<value_type>&& omega,
                                       std::size_t seed = std::mt19937_64::default_seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value);

        Vector<T> operator()(const Point<value_type>& space_point, value_type time_point) const override;

        ~KraichanGeneratorDeltaFunction() override = default;

    private:
        std::size_t n_;
        value_type k_0_, w_0_;
        std::vector<Vector<value_type>> v_, w_;
        std::vector<Vector<value_type>> k_;
        std::vector<value_type> omega_;


        static std::vector<Vector<value_type>> generate_wave_vectors(value_type k_0, std::size_t n, std::size_t seed);
        static std::vector<value_type> generate_frequencies(value_type w_0, std::size_t n, std::size_t seed);

        template<std::ranges::viewable_range WaveVectors>
        std::vector<Vector<T>> generate_ampltudes(const WaveVectors& wave_vectors, std::size_t n, std::size_t seed,
                                                  const Vector<value_type>& mean = Vector<value_type>{0, 0, 0},
                                                  const tensor::Tensor<value_type>& tensor = tensor::Tensor<value_type>{
                                                          1, 0, 0,
                                                          0, 1, 0,
                                                          0, 0, 1});
    };


    /*
        Uses Gaussian spectra, generates velocity vectors with given stddev
    */
    template<std::floating_point T>
    class KraichanGeneratorGaussian final : public AKraichanGenerator<T> {
    public:
        using value_type = T;

        KraichanGeneratorGaussian(std::size_t n, value_type k_0, value_type w_0, std::size_t seed = std::mt19937_64::default_seed) noexcept(std::is_nothrow_move_constructible_v<Vector<value_type>>&& std::is_nothrow_move_constructible_v<std::vector<Vector<value_type>>>);

        KraichanGeneratorGaussian(std::size_t n, value_type k_0, value_type w_0, std::vector<Vector<value_type>> k,
                                  std::size_t seed = std::mt19937_64::default_seed) noexcept(std::is_nothrow_move_constructible_v<Vector<value_type>>&&
                                                                                                     std::is_nothrow_move_constructible_v<std::vector<Vector<value_type>>>);

        KraichanGeneratorGaussian(std::size_t n, value_type k_0, value_type w_0,
                                  std::vector<Vector<value_type>>&& v,
                                  std::vector<Vector<value_type>>&& w,
                                  std::vector<Vector<value_type>>&& k,
                                  std::vector<value_type>&& omega) noexcept(std::is_nothrow_move_constructible_v<Vector<value_type>>&& std::is_nothrow_move_constructible_v<std::vector<Vector<value_type>>>);


        Vector<T> operator()(const Point<value_type>& space_point, value_type time_point) const override;
        ~KraichanGeneratorGaussian() override = default;

    private:
        std::size_t n_;
        value_type k_0_, w_0_;
        std::vector<Vector<value_type>> k_;
        std::vector<Vector<value_type>> v_, w_;
        std::vector<value_type> omega_;

        static std::vector<Vector<value_type>> generate_wave_vectors(value_type k_0, std::size_t n, std::size_t seed);
        static std::vector<value_type> generate_frequencies(value_type w_0, std::size_t n, std::size_t seed);

        template<std::ranges::viewable_range WaveVectors>
        std::vector<Vector<T>> generate_ampltudes(const WaveVectors& wave_vectors, std::size_t n, std::size_t seed,
                                                  const Vector<value_type>& mean = Vector<value_type>{0, 0, 0},
                                                  const tensor::Tensor<value_type>& tensor = tensor::Tensor<value_type>{
                                                          1, 0, 0,
                                                          0, 1, 0,
                                                          0, 0, 1});
    };


    /***********************************************\
    |   Kraichnan generator based of delta function |
    \***********************************************/


    template<std::floating_point T>
    KraichanGeneratorDeltaFunction<T>::KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0, std::size_t seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value)
        : KraichanGeneratorDeltaFunction{n, k_0, w_0, generate_wave_vectors(k_0, n, seed), seed} {
    }

    template<std::floating_point T>
    KraichanGeneratorDeltaFunction<T>::KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0, std::vector<Vector<value_type>> k, std::size_t seed)
        : KraichanGeneratorDeltaFunction{n, k_0, w_0,
                                         generate_ampltudes(k, n, seed),
                                         generate_ampltudes(k, n, seed),
                                         std::move(k),
                                         generate_frequencies(w_0, n, seed),
                                         seed} {}

    template<std::floating_point T>
    KraichanGeneratorDeltaFunction<T>::KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0,
                                                                      std::vector<Vector<value_type>>&& v,
                                                                      std::vector<Vector<value_type>>&& w,
                                                                      std::vector<Vector<value_type>>&& k,
                                                                      std::vector<value_type>&& omega,
                                                                      std::size_t seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value)
        : n_{n}, k_0_{k_0}, w_0_{w_0}, v_{std::move(v)}, w_{std::move(w)}, k_{std::move(k)}, omega_{std::move(omega)} {}


    template<std::floating_point T>
    Vector<T> KraichanGeneratorDeltaFunction<T>::operator()(const Point<value_type>& space_point, value_type time_point) const {
        return FluctuationFouirerSum<value_type>::calculate_sum(v_ | std::views::all,
                                                                w_ | std::views::all,
                                                                k_ | std::views::all,
                                                                omega_ | std::views::all,
                                                                space_point,
                                                                time_point);
    }

    template<std::floating_point T>
    template<std::ranges::viewable_range WaveVectors>
    std::vector<Vector<typename KraichanGeneratorDeltaFunction<T>::value_type>> KraichanGeneratorDeltaFunction<T>::generate_ampltudes(const WaveVectors& wave_vectors, std::size_t n, std::size_t seed,
                                                                                                                                      const Vector<value_type>& mean,
                                                                                                                                      const tensor::Tensor<value_type>& tensor) {
        random::MultidimenasionalGaussianGenerator<value_type> multidim_generator{mean, tensor};
        auto generator = [&multidim_generator](const Vector<T>& wave_vector) -> Vector<value_type> {
            const auto generated_xi = multidim_generator();
            return cross_product(generated_xi, wave_vector);
        };
        return wave_vectors |
               ranges::views::transform(generator) |
               ranges::to<std::vector<Vector<typename KraichanGeneratorDeltaFunction<T>::value_type>>>();
    }

    template<std::floating_point T>
    std::vector<Vector<typename KraichanGeneratorDeltaFunction<T>::value_type>> KraichanGeneratorDeltaFunction<T>::generate_wave_vectors(value_type k_0, std::size_t n, std::size_t seed) {
        RNGenerator generator{generator_engines::get_engine<std::mt19937_64>(seed), std::normal_distribution<>{0, k_0}};
        auto generate_random_vector = [&generator, k_0]() {
            auto x = generator();
            auto y = generator();
            auto z = generator();
            Vector<value_type> result{x, y, z};
            return scale_to_length(result, k_0);
        };

        return ranges::views::generate(generate_random_vector) |
               ranges::views::take(n) | ranges::to<std::vector<Vector<value_type>>>();
    }

    template<std::floating_point T>
    std::vector<typename KraichanGeneratorDeltaFunction<T>::value_type> KraichanGeneratorDeltaFunction<T>::generate_frequencies(value_type w_0, std::size_t n, std::size_t seed) {
        RNGenerator<std::mt19937_64, std::normal_distribution<>> gen{generator_engines::get_engine<std::mt19937_64>(seed),
                                                                     std::normal_distribution<value_type>{0, w_0}};
        return ranges::views::generate([&gen] {
                   return gen();
               }) |
               ranges::views::take(n) |
               ranges::to<std::vector>();
    }

    /**************************************************\
    |   Kraichnan generator based on gaussian function |
    \**************************************************/


    template<std::floating_point T>
    KraichanGeneratorGaussian<T>::KraichanGeneratorGaussian(std::size_t n, value_type k_0, value_type w_0, std::size_t seed) noexcept(std::is_nothrow_move_constructible_v<Vector<value_type>>&& std::is_nothrow_move_constructible_v<std::vector<Vector<value_type>>>)
        : KraichanGeneratorGaussian{n, k_0, w_0, generate_wave_vectors(k_0, n, seed), seed} {}

    template<std::floating_point T>
    KraichanGeneratorGaussian<T>::KraichanGeneratorGaussian(std::size_t n, value_type k_0, value_type w_0, std::vector<Vector<value_type>> k,
                                                            std::size_t seed) noexcept(std::is_nothrow_move_constructible_v<Vector<value_type>>&&
                                                                                               std::is_nothrow_move_constructible_v<std::vector<Vector<value_type>>>)
        : KraichanGeneratorGaussian{n, k_0, w_0, generate_ampltudes(k, n, seed), generate_ampltudes(k, n, seed), std::move(k), generate_frequencies(w_0, n, seed)} {}

    template<std::floating_point T>
    KraichanGeneratorGaussian<T>::KraichanGeneratorGaussian(std::size_t n, value_type k_0, value_type w_0,
                                                            std::vector<Vector<value_type>>&& v,
                                                            std::vector<Vector<value_type>>&& w,
                                                            std::vector<Vector<value_type>>&& k,
                                                            std::vector<value_type>&& omega) noexcept(std::is_nothrow_move_constructible_v<Vector<value_type>>&& std::is_nothrow_move_constructible_v<std::vector<Vector<value_type>>>)
        : n_{n}, k_0_{k_0}, w_0_{w_0}, v_{std::move(v)}, w_{std::move(w)}, k_{std::move(k)}, omega_{std::move(omega)} {}

    template<std::floating_point T>
    template<std::ranges::viewable_range WaveVectors>
    std::vector<Vector<T>> KraichanGeneratorGaussian<T>::generate_ampltudes(const WaveVectors& wave_vectors, std::size_t n, std::size_t seed,
                                                                            const Vector<value_type>& mean,
                                                                            const tensor::Tensor<value_type>& tensor) {

        random::MultidimenasionalGaussianGenerator<value_type> multidim_generator{mean, tensor};
        auto generate = [&multidim_generator](const Vector<T>& wave_vector) {
            const auto generated_xi = multidim_generator.template operator()();
            return cross_product(generated_xi, wave_vector);
        };
        return ranges::views::transform(generate) | ranges::to<std::vector<Vector<value_type>>>();
    }

    template<std::floating_point T>
    std::vector<Vector<typename KraichanGeneratorGaussian<T>::value_type>> KraichanGeneratorGaussian<T>::generate_wave_vectors(value_type k_0, std::size_t n, std::size_t seed) {
        RNGenerator generator{generator_engines::get_engine<std::mt19937_64>(seed), std::normal_distribution<>{0, k_0}};
        auto generate_random_vector = [&generator, k_0]() {
            auto x = generator();
            auto y = generator();
            auto z = generator();
            Vector<value_type> result{x, y, z};
            return scale_to_length(result, k_0);
        };

        return ranges::views::generate(generate_random_vector) |
               ranges::views::take(n) | ranges::to<std::vector<value_type>>();
    }

    template<std::floating_point T>
    std::vector<typename KraichanGeneratorGaussian<T>::value_type> KraichanGeneratorGaussian<T>::generate_frequencies(value_type w_0, std::size_t n, std::size_t seed) {
        RNGenerator<std::mt19937_64, std::normal_distribution<>> gen{generator_engines::get_engine<std::mt19937_64>(seed),
                                                                     std::normal_distribution<value_type>{0, w_0}};
        return ranges::views::generate([&gen] {
                   return gen();
               }) |
               ranges::views::take(n) |
               ranges::to<std::vector>();
    }

}// namespace stg::generators
#endif