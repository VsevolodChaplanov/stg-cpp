#ifndef STG_GENERATOR_BUILDER_HPP
#define STG_GENERATOR_BUILDER_HPP

#include <memory>
#include <stg_tensor/tensor.hpp>
#include <stg_random/generator_concept.hpp>
#include <stg_random/irn_generator.hpp>

namespace stg::generators {

  template<std::floating_point T>
  class SpectralGeneratorBuilder final {
  public:
    using value_type = T;

    template<concepts::GeneratorConcept RNGenerator>
    SpectralGeneratorBuilder& set_amplitude_generator(RNGenerator amplitude_generator);

    template<concepts::GeneratorConcept RNGenerator>
    SpectralGeneratorBuilder& set_wave_vector_generator(RNGenerator amplitude_generator);

    template<concepts::GeneratorConcept RNGenerator>
    SpectralGeneratorBuilder& set_frequencies_generator(RNGenerator amplitude_generator);

    SpectralGeneratorBuilder& set_amplitude_generator(stg::tensor::Tensor<value_type> reynolds_tensor);


  private:
    std::shared_ptr<IRNGenerator<value_type>> amplitude_generator_;
    std::shared_ptr<IRNGenerator<value_type>> wave_vector_generator_;
    std::shared_ptr<IRNGenerator<value_type>> frequencies_generator_;
    stg::tensor::Tensor<value_type> reynolds_tensor_;
  };
}

#endif //STG_GENERATOR_BUILDER_HPP
