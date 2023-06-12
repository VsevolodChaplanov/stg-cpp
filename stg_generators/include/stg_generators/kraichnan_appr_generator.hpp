#ifndef STG_STG_GENERATORS_KRAICHNAN_APPR_GENERATOR_HHP
#define STG_STG_GENERATORS_KRAICHNAN_APPR_GENERATOR_HHP

#include "geometry/geometry.hpp"
#include "i_spectral_generator.hpp"
#include "kraichnan_spectral_generator.hpp"
#include "spectras_base.hpp"
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <random>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <utility>
#include <vector>

namespace stg::generators {

    template<std::floating_point T>
    struct NodeParameters final {
        T k_0;
        T v_0;
        T w_0 = k_0 * v_0;
        std::size_t nfourier;
        std::size_t seed;
    };

    // Uses delta function
    template<std::floating_point T>
    class SpectralApproximatorGenerator final : public ISpectralGenerator<T> {
    public:
        SpectralApproximatorGenerator(std::vector<NodeParameters<T>> nodes)
            : generator_nodes_params_{(std::move(nodes))},
              generators_(generator_nodes_params_ | ranges::views::transform([](const NodeParameters<T>& node) {
                              return KraichanGeneratorDeltaFunction<T>{node.nfourier, node.k_0, node.w_0, node.seed};
                          }) |
                          ranges::to<std::vector<KraichanGeneratorDeltaFunction<T>>>()) {
        }

        Vector<T> operator()(const Point<T>& space_point, T time_point) const override {
            Vector<T> result{0, 0, 0};
            for (const auto& generator: generators_) {
                result += generator(space_point, time_point);
            }
            return result;
        }

        ~SpectralApproximatorGenerator() override = default;

    private:
        std::vector<NodeParameters<T>> generator_nodes_params_;
        std::vector<KraichanGeneratorDeltaFunction<T>> generators_;
    };
}// namespace stg::generators

#endif