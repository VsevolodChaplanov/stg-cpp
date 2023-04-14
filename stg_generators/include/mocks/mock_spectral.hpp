#ifndef STG_MOCK_SPECTRAL_HPP
#define STG_MOCK_SPECTRAL_HPP

#include <range/v3/view/generate.hpp>
#include <stg_generators/i_spectral_generator.hpp>
#include <stg_generators/spectral_generator_config.hpp>

namespace stg::generators::mock {

  class SpectralGeneratorMock final : public ISpectralGenerator<double> {
  public:
    Vector<double> operator()(const Point<double>& space, double time) const override {
      return {1, 1, 1};
    }

    ~SpectralGeneratorMock() override = default;
  };

  inline SpectralGeneratorConfig<double> MakeDefaultSpectralConfig(double k_init = 0.01, double k_end = 1.5, std::size_t nodes_n = 1000,
                                                            double length = 1.e-3, double time = 1.e-3) {
    std::vector<double> k_vec(nodes_n);
    std::generate(k_vec.begin(), k_vec.end(),
      [&] {
        static double init = k_init;
        static double dk = (k_end - k_init) / nodes_n;
        return init += dk;
    });
    std::vector<double> energy_array(nodes_n);
    const auto energy_func = [] (double k) {
      return 16 * sqrt(2. / std::numbers::pi) * k * k * k * k * exp( - 2. * k * k );
    };

    std::transform(k_vec.cbegin(), k_vec.cend(), energy_array.begin(), energy_func);

    return SpectralGeneratorConfig<double>{
      .energy_spectra_ = energy_array,
      .wave_vectors_spectra_ = k_vec,
      .turbulence_length_scale_ = length,
      .turbulence_time_scale_ = time,
      .fourier_nodes_number_ = nodes_n
    };
  }
}

#endif //STG_MOCK_SPECTRAL_HPP
