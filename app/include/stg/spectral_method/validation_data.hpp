#ifndef STG_APP_SPECTRAL_METHOD_VALIDATION_DATA_HPP
#define STG_APP_SPECTRAL_METHOD_VALIDATION_DATA_HPP

#include "geometry/geometry.hpp"
#include "mesh_builders/cube_fe_mesh.hpp"
#include "mesh_builders/mesh_builders.hpp"
#include "rtable/cube_vtk_saver.hpp"
#include "stg/spectral_method/data_loader.hpp"
#include "stg_tensor/tensor.hpp"
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <ranges>
#include <string_view>
#include <utility>
#include <vector>

namespace stg {

    template<std::floating_point T>
    class EnergyFunction final {
    public:
        EnergyFunction(T cube_edge_len, std::size_t n, std::function<T(T)>&& function)
            : cube_mesh_{mesh::CubeMeshBuilder<T>{cube_edge_len, n}.build()},
              energy_function_{std::move(function)} {}

        EnergyFunction(spectral::DataLoader loader, std::string_view filewithmesh, std::function<T(T)>&& function)
            : cube_mesh_{loader.load_mesh<T>(filewithmesh)},
              energy_function_{std::move(function)} {}

        void fill() {
            const std::size_t nvert = cube_mesh_->n_vertices();
            energy_.resize(nvert);

            for (const std::size_t i: std::views::iota(0ul, nvert)) {
                const auto vertex = cube_mesh_->relation_table()->vertex(i);
                const auto k_mod = std::sqrt(dot_product(vertex, vertex));
                const auto energy_value = energy_function_(k_mod);
                energy_[i] = energy_value;
            }
        }

        void save_to_file(std::string_view filename, std::string_view table_name = "TableName") const {
            mesh::VtkRectilinearGridSaver saver(filename);
            saver.save_mesh(cube_mesh_->relation_table());
            saver.save_scalar_data(energy_.cbegin(), energy_.cend(), table_name);
        }

    private:
        std::shared_ptr<stg::mesh::CubeFiniteElementsMesh<T>> cube_mesh_;
        std::function<T(T)> energy_function_;
        std::vector<T> energy_;
    };

    // for delta function
    template<std::floating_point T>
    class KraichanCovarianceFunction final {
    public:
        KraichanCovarianceFunction(T cube_edge_len, std::size_t n)
            : cube_mesh_{mesh::CubeMeshBuilder<T>{cube_edge_len, n}.build()} {}

        void fill(T k_value) {
            const std::size_t nvert = cube_mesh_->n_vertices();
            covariance_.resize(nvert);

            for (const std::size_t i: std::views::iota(0ul, nvert)) {
                const auto vertex = cube_mesh_->relation_table()->vertex(i);
                const auto covariance = covariance_function(k_value, vertex);
                covariance_[i] = covariance;
            }
        }

        void save_to_file(std::string_view filename, std::string_view table_name = "TableName") const {
            mesh::VtkRectilinearGridSaver saver(filename);
            saver.save_mesh(cube_mesh_->relation_table());
            saver.save_tensor_data(covariance_.cbegin(), covariance_.cend(), table_name);
        }

    private:
        std::shared_ptr<stg::mesh::CubeFiniteElementsMesh<T>> cube_mesh_;
        std::vector<tensor::Tensor<T>> covariance_;

        T covariance_function(T k, const Point<T>& space_point) const {
            const auto r_mod = std::sqrt(dot_product(space_point, space_point));
            if (r_mod == 0) { return 1.; }
            const auto kr = k * r_mod;
            const auto sinkr = std::sin(kr);
            const auto result = 2 * sinkr / kr;
            return result;
        }
    };
}// namespace stg

#endif