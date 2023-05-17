#ifndef STG_SPECTRAL_GENERATOR_HPP
#define STG_SPECTRAL_GENERATOR_HPP

#include <execution>
#include <vector>
#include <memory>
#include <range/v3/view/iota.hpp>
#include <stg_tensor/tensor.hpp>
#include <stg_random/rn_generator_impl.hpp>
#include <stg_random/generator_engines.hpp>
#include <geometry/geometry.hpp>
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

  template<std::floating_point T>
  class ISpectra {
  public:
    using value_type = T;
    virtual value_type operator()(value_type k) const = 0;
    virtual ~ISpectra() = default;
  };


  template<std::floating_point T>
  class VonKarmanSpectra final : public ISpectra<T> {
  public:
    using value_type = T;

    VonKarmanSpectra(value_type k_e, value_type k_eta, value_type k_cut) noexcept
      : k_e_{k_e}, k_eta_{k_eta}, k_cut_{k_cut} { }

    value_type operator()(value_type k) const override {
      const auto numerator = std::pow(k / k_e_, 4);
      const auto inner_denumerator_braces = std::pow(k / k_e_, 2);
      const auto denumerator = std::pow(1 + 2.4 * inner_denumerator_braces, 17. / 6.);
      return numerator / denumerator * f_cut(k) * f_eta(k);
    }

  protected:
    value_type k_e_, k_eta_, k_cut_;

    value_type f_eta(value_type k) const {
      return std::pow(std::numbers::e, - (12 * k / k_eta_) * (12 * k / k_eta_));
    }

    value_type f_cut(value_type k) const {
      const auto exp_arg = - std::pow(4 * std::max(k - k_cut_, 0) / k_cut_, 3);
      return std::pow(std::numbers::e, exp_arg);
    }
  };


  template<std::floating_point T>
  class SpectralGeneratorV2 final : public ISpectralGenerator<T> {
  private:
    struct NodesAmplitudes final {
      std::vector<T> p_;
      std::vector<T> q_;
    };
  public:
    using value_type = T;

    SpectralGeneratorV2(std::size_t seed_for_generators,
                        value_type ampl_mean, value_type ampl_std,
                        value_type wv_mean, value_type wv_std,
                        value_type freq_mean, value_type freq_std)
      : SpectralGeneratorV2(std::make_shared<RNGenerator<std::mt19937_64, std::normal_distribution<>>>(
                              engines::get_engine<std::mt19937_64>(seed_for_generators),
                              std::normal_distribution(ampl_mean, ampl_std)
                          ),
                          std::make_shared<RNGenerator<std::mt19937_64, std::normal_distribution<>>>(
                              engines::get_engine<std::mt19937_64>(seed_for_generators),
                              std::normal_distribution(wv_mean, wv_std)
                          ),
                          std::make_shared<RNGenerator<std::mt19937_64, std::normal_distribution<>>>(
                              engines::get_engine<std::mt19937_64>(seed_for_generators),
                              std::normal_distribution(freq_mean, freq_std)
                          ),
                          std::make_shared<RNGenerator<std::mt19937_64, std::uniform_real_distribution<>>>(
                              engines::get_engine<std::mt19937_64>(seed_for_generators),
                              std::uniform_real_distribution(0., 1.)
                          ))
    { }

    template<concepts::GeneratorConcept AmplitudeGenerator,
        concepts::GeneratorConcept FrequenciesGenerator,
        concepts::GeneratorConcept WaveVectorGenerator>
    SpectralGeneratorV2(std::shared_ptr<AmplitudeGenerator> amplitude_generator,
                      std::shared_ptr<FrequenciesGenerator> frequencies_generator,
                      std::shared_ptr<WaveVectorGenerator> wave_vector_generator,
                      std::shared_ptr<IRNGenerator<value_type>> rand_param_generator,
                      stg::tensor::Tensor<value_type> reynolds_tensor)
        : lower_triangular_{reynolds_tensor.cholesky()}
        , amplitude_generator_{std::move(amplitude_generator)}
        , wave_vector_generator_{std::move(wave_vector_generator)}
        , frequencies_generator_{std::move(frequencies_generator)}
        , rand_param_generator_{std::move(rand_param_generator)}
    { }

    void initialize_wave_vector_amplitudes(value_type k_min, value_type k_max, std::size_t n) {
      k_modules_.resize(n);
      const auto dk = (k_max - k_min) / (n - 1);
      std::ranges::generate(k_modules_, [k_m = k_min, dk] () {
        static std::size_t index = 0;
        return k_m + dk * index++;
      });
    }

    template<ranges::viewable_range Range>
    void initialize_wave_vector_amplitudes(Range&& values) {
      assert(k_modules_.empty());
      ranges::copy(values, ranges::back_inserter(k_modules_));
    }

    void initialize_random_coefficients() {
      ranges::generate(ranges::back_inserter(random_coeffs_), *rand_param_generator_);
    }

    auto initialize_spectra(std::shared_ptr<ISpectra<value_type>> spectra) {
      spectra_ = std::move(spectra);
    }

    void initialize_inner_generators() {
      assert(spectra_n_ == static_cast<std::size_t>(ranges::distance(k_modules_)));
      mode_fluctuations_generators_.resize(spectra_n_);

      for (const auto [index, k_module] : k_modules_) {
        const auto a_coeff = random_coeffs_[index];

        mode_fluctuations_generators_[index] = SpectralModeFluctuationGenerator::create(
            fourier_modes_n_,
            amplitude_generator_,
            frequencies_generator_,
            wave_vector_generator_,
            a_coeff, spectra_->operator()(k_module), k_module);
      }
    }

    Vector<T> operator()(const Point<T> &real_space_point, T time) const override {

    }

    ~SpectralGeneratorV2() override = default;

  private:

    /*
     * Generate u_m for spectra value E(k_m)
     * Main fluctuation is sum_m u_m;
     */
    struct SpectralModeFluctuationGenerator final {
      friend SpectralGeneratorV2;

      std::size_t fourier_modes_n_;
      std::vector<Vector<value_type>> wave_vectors_;
      std::vector<Vector<value_type>> p_vectors_;
      std::vector<Vector<value_type>> q_vectors_;
      std::vector<value_type> frequencies_;
      // to solve \sum_n \vec{p}^2 + \vec{q}^2 = 4E(k_m)

      template<stg::concepts::GeneratorConcept AmplitudeGenerator,
          stg::concepts::GeneratorConcept FrequenciesGenerator,
          stg::concepts::GeneratorConcept WaveVectorsGenarator,
          std::floating_point CoeffType,
          std::floating_point EnergyValueType,
          std::floating_point WaveVectorModuleType>
      static SpectralModeFluctuationGenerator create(
          std::size_t fourier_modes_n,
          std::shared_ptr<AmplitudeGenerator> ampl_generator,
          std::shared_ptr<FrequenciesGenerator> frequencies_generator,
          std::shared_ptr<WaveVectorsGenarator> wave_vectors_generator,
          CoeffType random_coeff, EnergyValueType energy, WaveVectorModuleType k_module) {
        SpectralModeFluctuationGenerator generator;
        auto freq_future = std::async(std::launch::async,
          &SpectralModeFluctuationGenerator::initialize_frequencies, generator, frequencies_generator);
        auto wave_vectors_future = std::async(
            &SpectralModeFluctuationGenerator::initialize_wave_vectors,
            generator, wave_vectors_generator, k_module);
        wave_vectors_future.wait();

        generator.initialize_amplitudes(ampl_generator, energy, random_coeff);

        freq_future.wait();

        return generator;
      }

      template<stg::concepts::GeneratorConcept FrequenciesGenerator>
      void initialize_frequencies(std::shared_ptr<FrequenciesGenerator> freq_generator) {
        frequencies_.resize(fourier_modes_n_);
        std::ranges::generate(frequencies_, *freq_generator);
      }

      template<stg::concepts::GeneratorConcept AmplitudesGenerator>
      void initialize_amplitudes(std::shared_ptr<AmplitudesGenerator> ampl_generator,
                                 value_type energy_value, value_type rand_coeff) {
        assert(!wave_vectors_.empty());
        p_vectors_.resize(fourier_modes_n_);
        q_vectors_.resize(fourier_modes_n_);

        const auto p_amplitude = std::sqrt(rand_coeff * energy_value * 4 / fourier_modes_n_);
        const auto q_amplitude = std::sqrt((1 - rand_coeff) * energy_value * 4 / fourier_modes_n_);

        std::ranges::transform(std::as_const(wave_vectors_), p_vectors_, [&] (const Vector<value_type>& wave_vector) {
          const value_type x = ampl_generator.operator()();
          const value_type y = ampl_generator.operator()();
          const value_type z = ampl_generator.operator()();
          const Vector<value_type> xi{x, y, z};
          const auto p_vector = cross_product(xi, wave_vector);
          return scale_to_length(p_vector, p_amplitude);
        });

        std::ranges::transform(std::as_const(wave_vectors_), q_vectors_, [&] (const Vector<value_type>& wave_vector) {
          const value_type x = ampl_generator.operator()();
          const value_type y = ampl_generator.operator()();
          const value_type z = ampl_generator.operator()();
          const Vector<value_type> xi{x, y, z};
          const auto q_vector = cross_product(xi, wave_vector);
          return scale_to_length(q_vector, q_amplitude);
        });
      }

      template<stg::concepts::GeneratorConcept WaveVectorsGenerator>
      void initialize_wave_vectors(std::shared_ptr<WaveVectorsGenerator> wave_vectors_generator, value_type k_module) {
        wave_vectors_.resize(fourier_modes_n_);
        std::ranges::generate(wave_vectors_, [gen = std::move(wave_vectors_generator), k_module, this] () {
          return generate_vector_with_length(gen, k_module);
        });
      }

      template<stg::concepts::GeneratorConcept Generator>
      Vector<value_type> generate_vector_with_length(std::shared_ptr<Generator> gen, value_type len) {
        value_type x = gen.operator()();
        value_type y = gen.operator()();
        value_type z = gen.operator()();
        Vector<value_type> vec{x, y, z};
        return scale_to_length(vec, len);
      };

      template<std::ranges::viewable_range SpectraFunction>
      NodesAmplitudes generate_fourier_modes_vectors_amplitudes(SpectraFunction&& function,
                                                                value_type a_m, std::size_t n) const {
        std::vector<value_type> p, q;
        assert(a_m >= 0);
        assert(a_m <= 1);
        assert(n == static_cast<std::size_t>(function));
        p = ranges::for_each(std::forward<SpectraFunction>(function),
                             [a_m, n] (const auto e_m) {
                               return std::sqrt(a_m * e_m * 4 / n);
                             });
        q = ranges::for_each(std::forward<SpectraFunction>(function),
                             [a_m, n] (const auto e_m) {
                               return std::sqrt((1 - a_m) * e_m * 4 / n);
                             });
        return { p ,q };
      }
    };

    const std::size_t spectra_n_;
    const std::size_t fourier_modes_n_;
    const tensor::Tensor<value_type> lower_triangular_;
    const std::shared_ptr<IRNGenerator<value_type>> amplitude_generator_;
    const std::shared_ptr<IRNGenerator<value_type>> wave_vector_generator_;
    const std::shared_ptr<IRNGenerator<value_type>> frequencies_generator_;
    const std::shared_ptr<IRNGenerator<value_type>> rand_param_generator_;
    std::shared_ptr<ISpectra<value_type>> spectra_;
    std::vector<value_type> k_modules_;
    std::vector<value_type> random_coeffs_;

    std::vector<SpectralModeFluctuationGenerator> mode_fluctuations_generators_;
  };
}

#endif //STG_SPECTRAL_GENERATOR_HPP
