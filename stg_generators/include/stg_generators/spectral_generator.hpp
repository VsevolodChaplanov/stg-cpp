#ifndef STG_SPECTRAL_GENERATOR_HPP
#define STG_SPECTRAL_GENERATOR_HPP

#include <execution>
#include <vector>
#include <memory>
#include <ranges>
#include <filesystem>
#include <algorithm>
#include <range/v3/view/iota.hpp>
#include <range/v3/all.hpp>
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
  class KolmogorovSpectra final : public ISpectra<T> {
  public:
    using value_type = T;
    value_type operator()(value_type k) const override {
      value_type logkappa = log10(k);
      value_type logE;
      if (logkappa < 0.0){
        logE = 2 * logkappa - 1;
      } else if (logkappa < 3.0){
        logE = -5.0/3.0 * logkappa - 1;
      } else {
        logE = -3 * logkappa + 3;
      }
      return std::pow(10, logE);
    }
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
      const auto exp_arg = - std::pow(4 * std::max(k - k_cut_, 0.) / k_cut_, 3);
      return std::pow(std::numbers::e, exp_arg);
    }
  };

  namespace fs = std::filesystem;
  template<std::floating_point T>
  struct SpectralParameters final {
    std::size_t seed = 42;
    T ampl_mean = 0., ampl_std = 1.;
    T wv_mean = 0., wv_std = 1.;
    T freq_mean = 0., freq_std = 1. / 2.;
    T k_min = 0.6, k_max = 32.;
    T length_scale = 1, time_scale = 1;
    std::size_t n_spectra = 100, n_fourier = 1000;
    std::size_t edge_points = 21;
    fs::path save_data_dir_path = fs::path{"./spectral_result/"};
    T cube_edge_len = 10.;
    stg::tensor::Tensor<T> reynolds_tensor_ {
        1., 0., 0.,
        0., 1., 0.,
        0., 0., 1.
    };
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

    explicit SpectralGeneratorV2(const SpectralParameters<value_type>& parameters)
      : SpectralGeneratorV2{parameters.seed, parameters.ampl_mean, parameters.ampl_std,
                            parameters.wv_mean, parameters.wv_std,
                            parameters.freq_mean, parameters.freq_std,
                            parameters.length_scale, parameters.time_scale,
                            parameters.reynolds_tensor_, parameters.n_spectra, parameters.n_fourier}
    { }

    SpectralGeneratorV2(std::size_t seed_for_generators,
                        value_type ampl_mean, value_type ampl_std,
                        value_type wv_mean, value_type wv_std,
                        value_type freq_mean, value_type freq_std,
                        value_type length_scale, value_type time_scale,
                        stg::tensor::Tensor<value_type> reynolds_tensor,
                        std::size_t spectra_n, std::size_t fourier_n)
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
                          ),
                          std::move(reynolds_tensor), length_scale, time_scale,
                          spectra_n, fourier_n)
    { }

    template<concepts::GeneratorConcept AmplitudeGenerator,
        concepts::GeneratorConcept FrequenciesGenerator,
        concepts::GeneratorConcept WaveVectorGenerator>
    SpectralGeneratorV2(std::shared_ptr<AmplitudeGenerator> amplitude_generator,
                      std::shared_ptr<FrequenciesGenerator> frequencies_generator,
                      std::shared_ptr<WaveVectorGenerator> wave_vector_generator,
                      std::shared_ptr<IRNGenerator<value_type>> rand_param_generator,
                      stg::tensor::Tensor<value_type> reynolds_tensor,
                      value_type length_scale, value_type time_scale,
                      std::size_t spectra_n, std::size_t fourier_n)
        : spectra_n_{spectra_n}, fourier_modes_n_{fourier_n}
        , lower_triangular_{reynolds_tensor.cholesky()}
        , length_scale_{length_scale}, time_scale_{time_scale}
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
    [[maybe_unused]] void initialize_wave_vector_amplitudes(Range&& values) {
      assert(k_modules_.empty());
      ranges::copy(values, ranges::back_inserter(k_modules_));
    }

    void initialize_random_coefficients() {
      random_coeffs_.resize(spectra_n_);
      std::generate(std::execution::par, random_coeffs_.begin(), random_coeffs_.end(),
                    [this] { return rand_param_generator_->operator()(); });
    }

    auto initialize_spectra(std::shared_ptr<ISpectra<value_type>> spectra) {
      spectra_ = std::move(spectra);
    }

    void initialize_inner_generators() {
      assert(spectra_n_ == static_cast<std::size_t>(ranges::distance(k_modules_)));
      mode_fluctuations_generators_.resize(spectra_n_);

      for (const auto [index, k_module] : k_modules_ | ranges::views::enumerate) {
        const auto a_coeff = random_coeffs_[index];

        mode_fluctuations_generators_[index] = SpectralModeFluctuationGenerator::create(
          fourier_modes_n_,
          amplitude_generator_,
          frequencies_generator_,
          wave_vector_generator_,
          a_coeff, spectra_->operator()(k_module), k_module,
          scale_factor_,
          { lower_triangular_.get(0, 0), lower_triangular_.get(1, 1), lower_triangular_.get(2, 2) }
        );
      }
    }

    Vector<T> operator()(const Point<T> &real_space_point, T time) const override {
      Vector<T> result{0., 0., 0.};
      Point<T> scaled_vert{real_space_point.template get<0>() / length_scale_,
                           real_space_point.template get<1>() / length_scale_,
                           real_space_point.template get<2>() / length_scale_};
      time /= time_scale_;
      for (const auto& inner_generator: mode_fluctuations_generators_) {
        result += inner_generator(scaled_vert, time);
      }

      return result;
    }

    value_type max_period() const {
      const auto min_omega = std::ranges::min(mode_fluctuations_generators_ | std::views::transform([](const auto& gen) {
        return std::ranges::min(gen.frequencies_ | std::views::transform([] (auto value) {
          return std::fabs(value);
        }));
      }));
      return 2 * std::numbers::pi / min_omega;
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

      Vector<value_type> operator()(const Point<value_type>& point, value_type time) const {
        Vector<value_type> result{0., 0., 0};

        for (const std::size_t node : std::views::iota(0ull, fourier_modes_n_)) {
          const auto phase = dot_product(wave_vectors_[node], point) + frequencies_[node] * time;
          result += p_vectors_[node] * std::cos(phase) + q_vectors_[node] * std::sin(phase);
        }

        return result * std::sqrt(2. / fourier_modes_n_);
      }

      template<stg::concepts::GeneratorConcept AmplitudeGenerator,
          stg::concepts::GeneratorConcept FrequenciesGenerator,
          stg::concepts::GeneratorConcept WaveVectorsGenerator,
          std::floating_point CoeffType,
          std::floating_point EnergyValueType,
          std::floating_point WaveVectorModuleType>
      static SpectralModeFluctuationGenerator create(
          std::size_t fourier_modes_n,
          std::shared_ptr<AmplitudeGenerator> ampl_generator,
          std::shared_ptr<FrequenciesGenerator> frequencies_generator,
          std::shared_ptr<WaveVectorsGenerator> wave_vectors_generator,
          CoeffType random_coeff, EnergyValueType energy, WaveVectorModuleType k_module,
          value_type scale_coeff, std::array<value_type, 3> diagonal) {
        SpectralModeFluctuationGenerator generator;
        generator.fourier_modes_n_ = fourier_modes_n;
        generator.initialize_frequencies(frequencies_generator);
        generator.initialize_wave_vectors(wave_vectors_generator, k_module, scale_coeff, diagonal);
        generator.initialize_amplitudes(ampl_generator, energy, random_coeff);
        return generator;
      }

      template<stg::concepts::GeneratorConcept FrequenciesGenerator>
      void initialize_frequencies(std::shared_ptr<FrequenciesGenerator> freq_generator) {
        frequencies_.resize(fourier_modes_n_);
        std::generate(std::execution::par, frequencies_.begin(), frequencies_.end(),
                      [gen = std::forward<std::shared_ptr<FrequenciesGenerator>>(freq_generator)]
                      { return gen->operator()(); });
      }

      template<stg::concepts::GeneratorConcept AmplitudesGenerator>
      void initialize_amplitudes(std::shared_ptr<AmplitudesGenerator> ampl_generator,
                                 value_type energy_value, value_type rand_coeff) {
        assert(!wave_vectors_.empty());
        p_vectors_.resize(fourier_modes_n_);
        q_vectors_.resize(fourier_modes_n_);

        const auto p_amplitude = std::sqrt(rand_coeff * energy_value * 4 / fourier_modes_n_);
        const auto q_amplitude = std::sqrt((1 - rand_coeff) * energy_value * 4 / fourier_modes_n_);

        auto createRandVector = [ampl_generator] {
          const value_type x = ampl_generator->operator()();
          const value_type y = ampl_generator->operator()();
          const value_type z = ampl_generator->operator()();
          return Vector<value_type>{x, y, z};
        };

        for (const std::size_t index : std::views::iota(0ull, fourier_modes_n_)) {
          auto xi = createRandVector();
          auto zeta = createRandVector();

          auto p_vector = cross_product(xi, wave_vectors_[index]);
          auto q_vector = cross_product(zeta, wave_vectors_[index]);

          p_vector = scale_to_length(p_vector, p_amplitude);
          q_vector = scale_to_length(q_vector, q_amplitude);

          p_vectors_[index] = std::move(p_vector);
          q_vectors_[index] = std::move(q_vector);
        }
      }

      template<stg::concepts::GeneratorConcept WaveVectorsGenerator>
      void initialize_wave_vectors(std::shared_ptr<WaveVectorsGenerator> wave_vectors_generator, value_type k_module,
                                   value_type scale_coeff, std::array<value_type, 3> diagonal) {
        wave_vectors_.resize(fourier_modes_n_);
        std::ranges::generate(wave_vectors_, [gen = std::move(wave_vectors_generator), k_module,
                                              scale_coeff, diagonal, this] () {
          auto generated = generate_vector_with_length(gen, k_module);
          generated = {generated.template get<0>() * scale_coeff / diagonal[0],
                       generated.template get<1>() * scale_coeff / diagonal[1],
                       generated.template get<2>() * scale_coeff / diagonal[2]};
          return generated;
        });
      }

      template<stg::concepts::GeneratorConcept Generator>
      Vector<value_type> generate_vector_with_length(std::shared_ptr<Generator> gen, value_type len) {
        assert(len > 0.);
        value_type x = gen->operator()();
        value_type y = gen->operator()();
        value_type z = gen->operator()();
        Vector<value_type> vec{x, y, z};
        return scale_to_length(vec, len);
      };
    };

    const std::size_t spectra_n_;
    const std::size_t fourier_modes_n_;
    const tensor::Tensor<value_type> lower_triangular_;
    const value_type length_scale_, time_scale_;
    const value_type scale_factor_ = length_scale_ / time_scale_;
    const std::shared_ptr<IRNGenerator<value_type>> amplitude_generator_;
    const std::shared_ptr<IRNGenerator<value_type>> wave_vector_generator_;
    const std::shared_ptr<IRNGenerator<value_type>> frequencies_generator_;
    const std::shared_ptr<IRNGenerator<value_type>> rand_param_generator_;
    std::shared_ptr<ISpectra<value_type>> spectra_
      = std::make_shared<VonKarmanSpectra<value_type>>(0.1, 10, 0.1);
    std::vector<value_type> k_modules_;
    std::vector<value_type> random_coeffs_;

    std::vector<SpectralModeFluctuationGenerator> mode_fluctuations_generators_;
  };
}

#endif //STG_SPECTRAL_GENERATOR_HPP
