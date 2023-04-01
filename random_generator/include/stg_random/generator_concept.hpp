#ifndef STG_RNGENERATOR_CONCEPTS_HPP
#define STG_RNGENERATOR_CONCEPTS_HPP

#include <concepts>
#include <utility>
#include <type_traits>

namespace stg::concepts {
  template<typename T>
  concept HasResultTypeConcept = requires (T t) { typename T::result_type; };
  template<typename T>
  concept HasEngineTypeConcept = requires (T t) { typename T::result_type; };

  template<typename Engine>
  concept EngineConcept = std::invocable<Engine> && HasEngineTypeConcept<Engine>;

  template<typename Func, typename Engine>
  concept DistributionConcept = EngineConcept<Engine>
                                && HasResultTypeConcept<Func>
                                && requires (Func f, Engine e) { f(e); };

  template<typename RNGenerator>
  concept GeneratorConcept = HasResultTypeConcept<RNGenerator>
                             && std::invocable<RNGenerator>
                             && std::is_floating_point_v<typename RNGenerator::result_type>;
}

#endif //STG_GENERATOR_CONCEPTS_HPP
