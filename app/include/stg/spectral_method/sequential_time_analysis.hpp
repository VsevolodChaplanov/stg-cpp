#ifndef STG_APP_SPECTRAL_METHOD_SEQUENTIAL_TIME_ANALYSIS_HPP
#define STG_APP_SPECTRAL_METHOD_SEQUENTIAL_TIME_ANALYSIS_HPP

#include "stg/spectral_method/data_loader.hpp"
#include <concepts>
#include <fmt/format.h>
namespace stg::spectral {

    template<std::floating_point T>
    class SquentialTimeAnalysis final {
    public:
        using value_type = T;

        SquentialTimeAnalysis(stg::spectral::DataLoader loader_)
            : loader_{std::move(loader_)} {}

    private:
        DataLoader loader_;
        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> velocity_mesh_;
    };
}// namespace stg::spectral

#endif