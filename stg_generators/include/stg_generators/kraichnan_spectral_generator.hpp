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
        KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0,
                                       random::MultidimenasionalGaussianGenerator<value_type> multidim_generator,
                                       RNGenerator<std::mt19937_64, std::normal_distribution<>> one_dim_distr,
                                       std::size_t seed = std::mt19937_64::default_seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value);
        KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0,
                                       random::MultidimenasionalGaussianGenerator<value_type> multidim_generator,
                                       RNGenerator<std::mt19937_64, std::normal_distribution<>> one_dim_distr,
                                       std::vector<Vector<T>> wave_vectors, std::vector<value_type> frequencies,
                                       std::size_t seed = std::mt19937_64::default_seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value);
        KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0,
                                       std::vector<Vector<value_type>>&& v,
                                       std::vector<Vector<value_type>>&& w,
                                       std::vector<Vector<value_type>>&& k,
                                       std::vector<value_type>&& omega,
                                       std::size_t seed = std::mt19937_64::default_seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value);

        Vector<T> operator()(const Point<value_type>& space_point, value_type time_point) override;

        ~KraichanGeneratorDeltaFunction() override = default;

    private:
        std::size_t n_;
        value_type k_0_, w_0_;
        std::vector<Vector<value_type>> k_;
        std::vector<Vector<value_type>> v_, w_;
        std::vector<value_type> omega_;


        static std::vector<Vector<T>> generate_wave_vectors(ISpectralGenerator<value_type>& one_dim_generator, value_type k_0, std::size_t n);

        template<std::ranges::viewable_range WaveVectors>
        static std::vector<Vector<T>> generate_ampltudes(random::MultidimenasionalGaussianGenerator<value_type>& multidim_generator,
                                                         const WaveVectors& wave_vectors, std::size_t n);

        static std::vector<Vector<T>> generate_frequencies(value_type w_0, std::size_t n, std::size_t seed);
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


        Vector<T> operator()(const Point<value_type>& space_point, value_type time_point) override;
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
        : KraichanGeneratorDeltaFunction{n, k_0, w_0,
                                         random::MultidimenasionalGaussianGenerator<value_type>{Vector<value_type>{0, 0, 0},
                                                                                                tensor::Tensor<value_type>{1, 0, 0, 0, 1, 0, 0, 0, 1},
                                                                                                seed},
                                         RNGenerator<std::mt19937_64, std::normal_distribution<>>{
                                                 generator_engines::get_engine<std::mt19937_64>(seed),
                                                 std::normal_distribution<value_type>{0, 1}}} {
    }

    template<std::floating_point T>
    KraichanGeneratorDeltaFunction<T>::KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0,
                                                                      random::MultidimenasionalGaussianGenerator<value_type> multidim_generator,
                                                                      RNGenerator<std::mt19937_64, std::normal_distribution<>> one_dim_distr,
                                                                      std::size_t seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value)
        : KraichanGeneratorDeltaFunction{n, k_0, w_0, std::move(multidim_generator), one_dim_distr,
                                         generate_wave_vectors(one_dim_distr, k_0, n),
                                         generate_frequencies(w_0, n, seed)} {
    }

    template<std::floating_point T>
    KraichanGeneratorDeltaFunction<T>::KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0,
                                                                      random::MultidimenasionalGaussianGenerator<value_type> multidim_generator,
                                                                      RNGenerator<std::mt19937_64, std::normal_distribution<>> one_dim_distr,
                                                                      std::vector<Vector<T>> wave_vectors, std::vector<value_type> frequencies,
                                                                      std::size_t seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value)
        : KraichanGeneratorDeltaFunction{n, k_0, w_0,
                                         generate_ampltudes(multidim_generator, wave_vectors, n),
                                         generate_ampltudes(multidim_generator, wave_vectors, n),
                                         wave_vectors,
                                         frequencies,
                                         seed} {
    }

    template<std::floating_point T>
    KraichanGeneratorDeltaFunction<T>::KraichanGeneratorDeltaFunction(std::size_t n, value_type k_0, value_type w_0,
                                                                      std::vector<Vector<value_type>>&& v,
                                                                      std::vector<Vector<value_type>>&& w,
                                                                      std::vector<Vector<value_type>>&& k,
                                                                      std::vector<value_type>&& omega,
                                                                      std::size_t seed) noexcept(std::is_nothrow_copy_constructible<Vector<value_type>>::value)
        : n_{n}, k_0_{k_0}, w_0_{w_0}, v_{std::move(v)}, w_{std::move(w)}, k_{std::move(k)}, omega_{std::move(omega)} {}


    template<std::floating_point T>
    Vector<T> KraichanGeneratorDeltaFunction<T>::operator()(const Point<value_type>& space_point, value_type time_point) {
        return FluctuationFouirerSum<value_type>::calculate_sum(v_ | std::views::all,
                                                                w_ | std::views::all,
                                                                k_ | std::views::all,
                                                                omega_ | std::views::all,
                                                                std::forward<Point>(space_point),
                                                                time_point);
    }

    template<std::floating_point T>
    std::vector<Vector<T>> KraichanGeneratorDeltaFunction<T>::generate_wave_vectors(ISpectralGenerator<KraichanGeneratorDeltaFunction<T>::value_type>& one_dim_generator, value_type k_0, std::size_t n) {
        std::vector<Vector<T>> result{n};

        auto generate_random_vector = [&one_dim_generator, k_0]() {
            auto x = one_dim_generator.template operator()();
            auto y = one_dim_generator.template operator()();
            auto z = one_dim_generator.template operator()();
            Vector<value_type> result{x, y, z};
            return scale_to_length(result, k_0);
        };

        std::ranges::generate(result, std::move(generate_random_vector));

        return result;
    }

    template<std::floating_point T>
    template<std::ranges::viewable_range WaveVectors>
    std::vector<Vector<T>> KraichanGeneratorDeltaFunction<T>::generate_ampltudes(random::MultidimenasionalGaussianGenerator<value_type>& multidim_generator,
                                                                                 const WaveVectors& wave_vectors, std::size_t n) {
        auto generate = [&multidim_generator](const Vector<T>& wave_vector) {
            const auto generated_xi = multidim_generator.template operator()();
            return cross_product(generated_xi, wave_vector);
        };

        return std::ranges::views::transform(wave_vectors, std::move(generate)) | ranges::to<std::vector<value_type>>;
    }

    template<std::floating_point T>
    std::vector<Vector<T>> KraichanGeneratorDeltaFunction<T>::generate_frequencies(value_type w_0, std::size_t n, std::size_t seed) {
        RNGenerator<std::mt19937_64, std::normal_distribution<>> gen{generator_engines::get_engine<std::mt19937_64>(seed),
                                                                     std::normal_distribution<value_type>{0, w_0}};
        return ranges::views::generate([&gen] {
                   return gen();
               }) |
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
        co_return ranges::views::transform(co_await std::async(std::move(generate))) | ranges::to<std::vector<Vector<value_type>>>();
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

        co_return ranges::views::generate(co_await std::async(generate_random_vector)) | ranges::to<std::vector<value_type>>();
    }

    template<std::floating_point T>
    std::vector<typename KraichanGeneratorGaussian<T>::value_type> KraichanGeneratorGaussian<T>::generate_frequencies(value_type w_0, std::size_t n, std::size_t seed) {
        RNGenerator<std::mt19937_64, std::normal_distribution<>> gen{generator_engines::get_engine<std::mt19937_64>(seed),
                                                                     std::normal_distribution<value_type>{0, w_0}};
        return ranges::views::generate([&gen] {
                   return gen();
               }) |
               ranges::to<std::vector>();
    }

}// namespace stg::generators
#endif