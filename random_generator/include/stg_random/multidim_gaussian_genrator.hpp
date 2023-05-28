#ifndef STG_RANGOM_GENERATOR_MULTIDIM_GAUSSIAN_GENERATOR_HPP
#define STG_RANGOM_GENERATOR_MULTIDIM_GAUSSIAN_GENERATOR_HPP

#include "generator_engines.hpp"
#include "irn_generator.hpp"
#include <algorithm>
#include <armadillo>
#include <concepts>
#include <cstddef>
#include <geometry/geometry.hpp>
#include <random>
#include <stg_tensor/tensor.hpp>
#include <type_traits>

namespace stg::random {

    // Трёхмерное гауссово распределение с заданным средним и ковариационной матрицей
    template<std::floating_point T>
    class MultidimenasionalGaussianGenerator final {
    public:
        using distribution = std::normal_distribution<T>;
        using engine_type = std::mt19937_64;
        using value_type = T;
        using result_type = Vector<value_type>;

        MultidimenasionalGaussianGenerator(Vector<T> mean, tensor::Tensor<T> covariance, std::size_t seed = std::mt19937_64::default_seed) noexcept(
                std::is_nothrow_copy_constructible<Vector<T>>::value&&
                        std::is_nothrow_copy_constructible<tensor::Tensor<T>>::value);


        MultidimenasionalGaussianGenerator(const arma::vec3& mean, const arma::mat33& covariance, std::size_t seed = std::mt19937_64::default_seed) noexcept(
                std::is_nothrow_copy_constructible<arma::vec3>::value&&
                        std::is_nothrow_copy_constructible<arma::mat33>::value);

        result_type operator()();

        ~MultidimenasionalGaussianGenerator() = default;

    private:
        arma::vec3 mean_;
        arma::mat33 covariance_ = arma::mat::fixed<3, 3>();
    };


    template<std::floating_point T>
    MultidimenasionalGaussianGenerator<T>::MultidimenasionalGaussianGenerator(Vector<T> mean, tensor::Tensor<T> covariance, std::size_t seed) noexcept(
            std::is_nothrow_copy_constructible<Vector<T>>::value&&
                    std::is_nothrow_copy_constructible<tensor::Tensor<T>>::value)
        : MultidimenasionalGaussianGenerator<T>{arma::vec::fixed<3>({mean.template get<0>(), mean.template get<0>(), mean.template get<0>()}), arma::mat::fixed<3, 3>(covariance.cbegin())} {}

    template<std::floating_point T>
    MultidimenasionalGaussianGenerator<T>::MultidimenasionalGaussianGenerator(const arma::vec3& mean, const arma::mat33& covariance, std::size_t seed) noexcept(
            std::is_nothrow_copy_constructible<arma::vec3>::value&&
                    std::is_nothrow_copy_constructible<arma::mat33>::value)
        : mean_{mean}, covariance_{covariance} {
        arma::arma_rng::set_seed(seed);
    }

    template<std::floating_point T>
    MultidimenasionalGaussianGenerator<T>::result_type MultidimenasionalGaussianGenerator<T>::operator()() {
        arma::vec3 arma_result = arma::mvnrnd(mean_, covariance_);
        return {arma_result.at(0), arma_result.at(1), arma_result.at(2)};
    }
}// namespace stg::random

#endif