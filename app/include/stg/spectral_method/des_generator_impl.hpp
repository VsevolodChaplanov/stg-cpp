#ifndef STG_APP_SPECTRAL_METHOD_DES_GENERATOR_HPP
#define STG_APP_SPECTRAL_METHOD_DES_GENERATOR_HPP

#include "data_loader.hpp"
#include "mesh_builders/cube_fe_mesh.hpp"
#include "mesh_builders/mesh_builders.hpp"
#include "rtable/cube_vtk_saver.hpp"
#include "velocity_field/velocity_field.hpp"
#include "velocity_field/velocity_samples.hpp"
#include <concepts>
#include <functional>
#include <stg_generators/des_method.hpp>

namespace stg::spectral {

    template<std::floating_point T>
    class DesGeneratorImpl final {
    public:
        using value_type = T;

        DesGeneratorImpl(value_type cube_edge_length, std::size_t cube_edge_n, std::size_t samples_amount,
                         std::size_t n_fourier, value_type omega_0,
                         std::function<value_type(value_type)> energy_function,
                         std::size_t seed = std::mt19937_64::default_seed)
            : fe_mesh_{CubeMeshBuilder<value_type>{cube_edge_length, cube_edge_n}.build()},
              generator_{n_fourier,
                         cube_edge_length,
                         cube_edge_length / static_cast<value_type>(cube_edge_n),
                         omega_0,
                         std::move(energy_function),
                         seed},
              velocity_field_{fe_mesh_->n_vertices()} {}

        void save_sample(const std::filesystem::path& directory, std::string_view filename, std::string_view table_name = "VelocityField") const {
            const std::string full_filename = fmt::format("{}/{}", directory.string(), filename);
            VtkRectilinearGridSaver saver{full_filename};
            saver.save_mesh(fe_mesh_->relation_table());
            saver.save_vector_data(velocity_field_.values_view(), table_name);
        }

        void generate_sample(value_type time) {
            const std::size_t nvert = fe_mesh_->n_vertices();

            for (const std::size_t ivert: std::views::iota(0ull, nvert)) {
                const auto vertex = fe_mesh_->relation_table()->vertex(ivert);
                auto fluct = generator_(vertex, time);
                velocity_field_.set_value(std::move(fluct), ivert);
            }
        }

    private:
        std::shared_ptr<CubeFiniteElementsMesh<value_type>> fe_mesh_;
        generators::DirectEnergySpectralGenerator<value_type> generator_;
        VelocityField<value_type> velocity_field_;
    };
}// namespace stg::spectral

#endif