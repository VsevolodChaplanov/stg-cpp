#ifndef STG_GENERATOR_ENGINES_HPP
#define STG_GENERATOR_ENGINES_HPP

#include "generator_concept.hpp"
#include <concepts>
#include <random>

namespace stg::generator_engines {
    /*
     * Template of common solution to have several refernces to
     * created by function engine type
     */
    template<concepts::GeneratorConcept Engine, std::invocable Device>
    static Engine& get_engine(Device device) {
        thread_local static Engine engine{device()};
        return engine;
    }

    template<typename Engine>
    static Engine& get_engine(size_t seed) {
        thread_local static Engine engine{seed};
        return engine;
    }

    inline static std::random_device& get_random_device() {
        thread_local static std::random_device rd;
        return rd;
    }
}// namespace stg::generator_engines


#endif