#ifndef STG_SPECTRAL_GENERATOR_HPP
#define STG_SPECTRAL_GENERATOR_HPP

#include <execution>
#include <vector>
#include <memory>
#include <range/v3/view/iota.hpp>
#include <stg_tensor/tensor.hpp>
#include <stg_random/rn_generator_impl.hpp>
#include <stg_random/generator_engines.hpp>
#include "i_spectral_generator.hpp"
#include "spectral_generator_config.hpp"

namespace stg::generators {
namespace engines = stg::generator_engines;

  template<std::floating_point T, std::size_t seed = 42>
  class SpectralGenerator final : ISpectralGenerator<T> {
  public:
    using value_type = typename ISpectralGenerator<T>::value_type;

    /*
     * Create generator with given rngenerators
     * values in generators can be not binded to
     * config values
     */
    template<concepts::GeneratorConcept AmplitudeGenerator,
      concepts::GeneratorConcept FrequenciesGenerator,
      concepts::GeneratorConcept WaveVectorGenerator>
    SpectralGenerator(std::shared_ptr<AmplitudeGenerator> amplitude_generator,
                      std::shared_ptr<FrequenciesGenerator> frequencies_generator,
                      std::shared_ptr<WaveVectorGenerator> wave_vector_generator,
                      SpectralGeneratorConfig<T> config)
      : config_{std::move(config)}
      , lower_triangular_{config_.reynolds_matrix_.cholesky()}
      , amplitude_generator_{std::move(amplitude_generator)}
      , wave_vector_generator_{std::move(wave_vector_generator)}
      , frequencies_generator_{std::move(frequencies_generator)}
      , wave_vectors_{wave_vectors(
          config_.wave_vectors_spectra_.cbegin(),
          config_.wave_vectors_spectra_.cend(),
          config_.fourier_nodes_number_
        )}
      , frequencies_{frequencies(config_.fourier_nodes_number_)}
    { }

    /*
     * Create generator using only config values
     */
    explicit SpectralGenerator(SpectralGeneratorConfig<T> config)
      : SpectralGenerator(std::make_shared<RNGenerator<std::mt19937_64, std::normal_distribution<>>>(
                            engines::get_engine<std::mt19937_64>(seed),
                            std::normal_distribution(config.amplitudes_generator_mean_, config.amplitudes_generator_std_)
                          ),
                          std::make_shared<RNGenerator<std::mt19937_64, std::normal_distribution<>>>(
                            engines::get_engine<std::mt19937_64>(seed),
                            std::normal_distribution(config.wave_vectors_generator_mean_, config.wave_vectors_generator_std_)
                          ),
                          std::make_shared<RNGenerator<std::mt19937_64, std::normal_distribution<>>>(
                            engines::get_engine<std::mt19937_64>(seed),
                            std::normal_distribution(config.frequencies_generator_mean_, config.frequencies_generator_std_)
                          ),
                          std::move(config))
    { }

    Vector<value_type> operator()(const Point<value_type>& real_space_point, value_type time) const override {
      const auto scaled_vertex = real_space_point / config_.turbulence_length_scale_;
      const auto scaled_time = time / config_.turbulence_time_scale_;

      Vector<value_type> result_fluctuation{0., 0., 0.};

      for (const size_t n: ranges::view::iota(0ul, config_.fourier_nodes_number_)) {
        const size_t phase = dot_product(wave_vectors_[n], scaled_vertex)
          + frequencies_[n] * scaled_time;
        const auto& wave_vector = wave_vectors_[n];
        const auto p_ampl = this->generate_amplitude(wave_vector);
        const auto q_ampl = this->generate_amplitude(wave_vector);
        result_fluctuation += p_ampl * std::cos(phase) + q_ampl * std::sin(phase);
      }

      result_fluctuation *= std::sqrt(2. / config_.fourier_nodes_number_);
      result_fluctuation = lower_triangular_ * result_fluctuation;

      return result_fluctuation;
    }

    ~SpectralGenerator() override = default;

  private:
    const SpectralGeneratorConfig<value_type> config_;
    const tensor::Tensor<value_type> lower_triangular_;
    const std::shared_ptr<IRNGenerator<value_type>> amplitude_generator_;
    const std::shared_ptr<IRNGenerator<value_type>> wave_vector_generator_;
    const std::shared_ptr<IRNGenerator<value_type>> frequencies_generator_;
    const std::vector<Vector<value_type>> wave_vectors_;
    const std::vector<value_type> frequencies_;

    std::vector<value_type> frequencies(std::size_t size) const {
      std::vector<value_type> result_frequencies(size);
      for (size_t i = 0; i < size; ++i) {
        result_frequencies[i] = frequencies_generator_->operator()();
      }
      return result_frequencies;
    }

    template<std::forward_iterator Iter>
    std::vector<Vector<value_type>> wave_vectors(Iter begin, Iter end, std::size_t size) const {
      const auto get_k_vector = [&] (const auto& spectra_value) {
        value_type x = wave_vector_generator_->operator()();
        value_type y = wave_vector_generator_->operator()();
        value_type z = wave_vector_generator_->operator()();

        const auto length = std::sqrt(x * x + y * y + z * z);
        const auto scaling = spectra_value / length;

        x *= scaling;
        y *= scaling;
        z *= scaling;

        return Vector<value_type>{x, y, z};
      };

      std::vector<Vector<value_type>> result(size);
      std::transform(std::execution::par_unseq, begin, end, result.begin(), get_k_vector);

      return result;
    }

    Vector<value_type> generate_amplitude(const Vector<value_type>& wave_vector) const {
      const value_type x = amplitude_generator_->operator()();
      const value_type y = amplitude_generator_->operator()();
      const value_type z = amplitude_generator_->operator()();
      const Vector<value_type> zeta{x, y, z};
      auto result = cross_product(zeta, wave_vector);
      return result;
    }
  };
}

#endif //STG_SPECTRAL_GENERATOR_HPP
