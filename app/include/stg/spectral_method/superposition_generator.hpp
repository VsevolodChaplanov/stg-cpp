#ifndef STF_APP_SPECTRAL_SUPERPOSITION_GENERATOR_HPP
#define STF_APP_SPECTRAL_SUPERPOSITION_GENERATOR_HPP

#include "mesh_builders/cube_fe_mesh.hpp"
#include "mesh_builders/mesh_builders.hpp"
#include "rtable/cube_vtk_saver.hpp"
#include "stg_generators/kraichnan_appr_generator.hpp"
#include "velocity_field/velocity_field.hpp"
#include <concepts>
namespace stg::spectral {

    template<std::floating_point T>
    class SuperpositionGenerator final {
    public:
        using value_type = T;

        SuperpositionGenerator(value_type cube_edge_length, std::size_t cube_edge_n, generators::SpectralApproximatorGenerator<T> generator);

        void save_sample(const std::filesystem::path& directory, std::string_view filename, std::string_view table_name = "VelocityField") const;
        void generate_sample(value_type time);

    private:
        std::shared_ptr<mesh::CubeFiniteElementsMesh<value_type>> fe_mesh_;
        generators::SpectralApproximatorGenerator<value_type> generator_;
        field::VelocityField<value_type> velocity_field_;
    };

    template<std::floating_point T>
    SuperpositionGenerator<T>::SuperpositionGenerator(value_type cube_edge_length, std::size_t cube_edge_n, generators::SpectralApproximatorGenerator<T> generator)
        : fe_mesh_{mesh::CubeMeshBuilder<value_type>{cube_edge_length, cube_edge_n}.build()},
          generator_{std::move(generator)},
          velocity_field_{fe_mesh_->n_vertices()} {}

    template<std::floating_point T>
    void SuperpositionGenerator<T>::save_sample(const std::filesystem::path& directory, std::string_view filename, std::string_view table_name) const {
        const std::string full_filename = fmt::format("{}/{}", directory.string(), filename);
        mesh::VtkRectilinearGridSaver saver{full_filename};
        saver.save_mesh(fe_mesh_->relation_table());
        saver.save_vector_data(velocity_field_.values_view(), table_name);
    }

    template<std::floating_point T>
    void SuperpositionGenerator<T>::generate_sample(value_type time) {
        const std::size_t nvert = fe_mesh_->n_vertices();

        for (const std::size_t ivert: std::views::iota(0ull, nvert)) {
            const auto vertex = fe_mesh_->relation_table()->vertex(ivert);
            auto fluct = generator_(vertex, time);
            velocity_field_.set_value(std::move(fluct), ivert);
        }
    }


}// namespace stg::spectral

#endif