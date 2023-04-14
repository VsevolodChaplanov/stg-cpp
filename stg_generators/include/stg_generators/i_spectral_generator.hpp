#ifndef STG_I_SPECTRAL_GENERATOR_HPP
#define STG_I_SPECTRAL_GENERATOR_HPP

#include "i_fluctuation_generator.hpp"

namespace stg::generators {

  template<std::floating_point T>
  class ISpectralGenerator : public IVelocityFluctuationGenerator<T> {
  public:
    using typename IVelocityFluctuationGenerator<T>::value_type;
    using IVelocityFluctuationGenerator<T>::operator();
    virtual ~ISpectralGenerator() = default;
  };
}

#endif //STG_I_SPECTRAL_GENERATOR_HPP
