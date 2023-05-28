#ifndef STG_RNGENERATOR_CONCEPTS_HPP
#define STG_RNGENERATOR_CONCEPTS_HPP

#include <concepts>
#include <type_traits>
#include <utility>

namespace stg::concepts {

    template<typename Engine>
    concept EngineConcept = std::invocable<Engine>;
    template<typename Func, typename Engine>
    concept DistributionConcept = EngineConcept<Engine> && requires(Func f, Engine e) { f(e); };

    template<typename RNGenerator>
    concept GeneratorConcept = std::is_floating_point_v<typename RNGenerator::result_type>;
}// namespace stg::concepts

#endif//STG_GENERATOR_CONCEPTS_HPP
