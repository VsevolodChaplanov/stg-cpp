#ifndef STG_RN_GENERATOR_IMPL_HPP
#define STG_RN_GENERATOR_IMPL_HPP

#include "generator_concept.hpp"
#include "irn_generator.hpp"

namespace stg {
  template<concepts::EngineConcept Engine, concepts::DistributionConcept<Engine> Distribution>
  class RNGenerator final : public IRNGenerator<typename Distribution::result_type> {
  public:
    using engine_type = Engine;
    using distribution_type = Distribution;
    using result_type = typename Distribution::result_type;

    RNGenerator(Engine& engine, Distribution distribution) noexcept
      : engine_{engine}, distribution_{std::move(distribution)}
    { }

    result_type operator()() override { return distribution_(engine_); }

    ~RNGenerator() override = default;

  private:
    Engine& engine_;
    Distribution distribution_;
  };
}
#endif //STG_RN_GENERATOR_IMPL_HPP
