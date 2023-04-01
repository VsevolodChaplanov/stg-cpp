#include "common.hpp"

struct SpectralGeneratorFixture {


};

SCENARIO_METHOD(SpectralGeneratorFixture, "Smirnov generator tests with default config") {
  const auto get_k_spectra = [] {
    const double k_init = 0.01;
    const double k_last = 1.5;
    const size_t fourier_size = 1000;
    std::vector<double> k_spectra(fourier_size);

    const double kh = (k_last - k_init) / fourier_size;
    double k_val = k_init;
    k_spectra.reserve(fourier_size);
    for (const size_t i : ranges::view::iota(0ul, fourier_size)) {
      k_spectra[i] = k_val;
      k_val += kh;
    }
    return k_spectra;
  };


  GIVEN("Discretized energy spectra and config") {
    SpectralGeneratorConfig<double> config{
      .energy_spectra_ = {},
      .wave_vectors_spectra_ = get_k_spectra()
    };

    SpectralGenerator<double, seed> spectral_generator{ std::move(config) };
    const auto generated_value = spectral_generator({0, 0, 0}, 0);
  }
}