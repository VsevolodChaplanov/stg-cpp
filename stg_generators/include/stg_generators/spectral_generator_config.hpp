#ifndef STG_SPECTRAL_GENERATOR_CONFIG_HPP
#define STG_SPECTRAL_GENERATOR_CONFIG_HPP

namespace stg::generators {
  template<std::floating_point T>
  struct SpectralGeneratorConfig final {
    using value_type = T;
    const std::vector<T> energy_spectra_;
    const std::vector<T> wave_vectors_spectra_;
    const T turbulence_length_scale_ = 1.e-3;
    const T turbulence_time_scale_ = 1.e-3;
    const T scale_factor_ = turbulence_length_scale_ / turbulence_time_scale_;
    const std::size_t fourier_nodes_number_ = 1000;
    const stg::tensor::Tensor<T> reynolds_matrix_{1, 0, 0, 0, 1, 0, 0, 0, 1};
    const T amplitudes_generator_mean_ = 0.;
    const T amplitudes_generator_std_ = .5;
    const T frequencies_generator_mean_ = 0.;
    const T frequencies_generator_std_ = 1.;
    const T wave_vectors_generator_mean_ = 0.;
    const T wave_vectors_generator_std_ = 1.;
  };
}

#endif //STG_SPECTRAL_GENERATOR_CONFIG_HPP
