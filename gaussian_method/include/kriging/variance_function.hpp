#ifndef STG_VARIANCE_FUNCTION_HPP
#define STG_VARIANCE_FUNCTION_HPP

#include <concepts>
#include <fem/fe_mesh.hpp>
#include <memory>
#include <mesh_builders/cube_fe_mesh.hpp>
#include <range/v3/all.hpp>

namespace stg::gaussian {
    using namespace stg::mesh;

    namespace {
        template<std::floating_point T>
        constexpr std::size_t center_index(const std::shared_ptr<mesh::CubePrizmFEMesh<T>>& space) {
            const std::size_t nx = space->nx();
            const std::size_t ny = space->ny();
            const std::size_t nz = space->nz();
            const auto i = static_cast<std::size_t>((nx - 1) / 2);
            const auto j = static_cast<std::size_t>((ny - 1) / 2);
            const auto k = static_cast<std::size_t>((nz - 1) / 2);
            return space->lin_index(i, j, k);
        }
    }// namespace

    template<std::floating_point T>
    class VarianceFunction final {
    public:
        using value_type = T;

        VarianceFunction(std::shared_ptr<CubePrizmFEMesh<T>> space,
                         std::vector<T> val11,
                         std::vector<T> val12,
                         std::vector<T> val13,
                         std::vector<T> val22,
                         std::vector<T> val23,
                         std::vector<T> val33)
            : space_(space), val11_{std::move(val11)}, val12_{std::move(val12)}, val13_{std::move(val13)}, val22_{std::move(val22)}, val23_{std::move(val23)}, val33_{std::move(val33)}, center_{{val11_[center_index(space)],
                                                                                                                                                                                                  val12_[center_index(space)],
                                                                                                                                                                                                  val13_[center_index(space)],
                                                                                                                                                                                                  val22_[center_index(space)],
                                                                                                                                                                                                  val23_[center_index(space)],
                                                                                                                                                                                                  val33_[center_index(space)]}} {}

    private:
        const std::shared_ptr<mesh::CubePrizmFEMesh<value_type>> space_;
        const std::vector<value_type> val11_;
        const std::vector<value_type> val12_;
        const std::vector<value_type> val13_;
        const std::vector<value_type> val22_;
        const std::vector<value_type> val23_;
        const std::vector<value_type> val33_;
        const std::array<value_type, 6> center_;
    };
}// namespace stg::gaussian

#endif//STG_VARIANCE_FUNCTION_HPP
