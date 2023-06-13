#ifndef STG_APP_KRAICHNAN_METHOD_IMPL_HPP
#define STG_APP_KRAICHNAN_METHOD_IMPL_HPP

#include "data_loader.hpp"
#include "mesh_builders/cube_fe_mesh.hpp"
#include "mesh_builders/mesh_builders.hpp"
#include "rtable/cube_vtk_saver.hpp"
#include "velocity_field/velocity_field.hpp"
#include "velocity_field/velocity_samples.hpp"
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <fmt/core.h>
#include <memory>
#include <random>
#include <ranges>
#include <stg_generators/kraichnan_spectral_generator.hpp>
#include <stg_thread_pool.hpp>
#include <string_view>

namespace stg::spectral {

    template<std::floating_point T>
    class KraichanMethodImpl final {
    public:
        using value_type = T;
        KraichanMethodImpl(value_type cube_edge_length, std::size_t cube_edge_n, std::size_t samples_amount,
                           std::size_t n_fourier, value_type k_0, value_type w_0, std::size_t seed = std::mt19937_64::default_seed);

        void save_samples(const std::filesystem::path& directory, std::string_view filename, std::string_view table_name = "VelocityField") const;
        void generate_sample(value_type time);
        void set_vo(value_type v0);

    private:
        std::shared_ptr<CubeFiniteElementsMesh<value_type>> fe_mesh_;
        generators::KraichanGeneratorDeltaFunction<value_type> generator_;
        VelocityField<value_type> velocity_field_;
    };

    template<std::floating_point T>
    class KraichanMethodGaussianSpectraImpl final {
    public:
        using value_type = T;
        KraichanMethodGaussianSpectraImpl(value_type cube_edge_length, std::size_t cube_edge_n, std::size_t samples_amount,
                                          std::size_t n_fourier, value_type k_0, value_type w_0, std::size_t seed = std::mt19937_64::default_seed);

        void save_samples(const std::filesystem::path& directory, std::string_view filename, std::string_view table_name = "VelocityField") const;
        void generate_sample(value_type time);

    private:
        std::shared_ptr<CubeFiniteElementsMesh<value_type>> fe_mesh_;
        generators::KraichanGeneratorGaussian<value_type> generator_;
        VelocityField<value_type> velocity_field_;
    };


    template<std::floating_point T>
    KraichanMethodImpl<T>::KraichanMethodImpl(value_type cube_edge_length, std::size_t cube_edge_n, std::size_t samples_amount,
                                              std::size_t n_fourier, value_type k_0, value_type w_0, std::size_t seed)
        : fe_mesh_{CubeMeshBuilder<value_type>{cube_edge_length, cube_edge_n}.build()},
          generator_{n_fourier, k_0, w_0, seed},
          velocity_field_{fe_mesh_->n_vertices()} {}

    template<std::floating_point T>
    void KraichanMethodImpl<T>::save_samples(const std::filesystem::path& directory, std::string_view filename, std::string_view table_name) const {
        const std::string full_filename = fmt::format("{}/{}", directory.string(), filename);
        VtkRectilinearGridSaver saver{full_filename};
        saver.save_mesh(fe_mesh_->relation_table());
        saver.save_vector_data(velocity_field_.values_view(), table_name);
    }

    template<std::floating_point T>
    void KraichanMethodImpl<T>::generate_sample(value_type time) {
        const std::size_t nvert = fe_mesh_->n_vertices();

        for (const std::size_t ivert: std::views::iota(0ull, nvert)) {
            const auto vertex = fe_mesh_->relation_table()->vertex(ivert);
            auto fluct = generator_(vertex, time);
            velocity_field_.set_value(std::move(fluct), ivert);
        }
    }

    template<std::floating_point T>
    void KraichanMethodImpl<T>::set_vo(value_type v0) {
        generator_.set_v0(v0);
    }


    template<std::floating_point T>
    KraichanMethodGaussianSpectraImpl<T>::KraichanMethodGaussianSpectraImpl(value_type cube_edge_length, std::size_t cube_edge_n, std::size_t samples_amount,
                                                                            std::size_t n_fourier, value_type k_0, value_type w_0, std::size_t seed)
        : fe_mesh_{CubeMeshBuilder<value_type>{cube_edge_length, cube_edge_n}.build()},
          generator_{n_fourier, k_0, w_0, seed},
          velocity_field_{fe_mesh_->n_vertices()} {}

    template<std::floating_point T>
    void KraichanMethodGaussianSpectraImpl<T>::save_samples(const std::filesystem::path& directory, std::string_view filename, std::string_view table_name) const {
        const std::string full_filename = fmt::format("{}/{}", directory.string(), filename);
        VtkRectilinearGridSaver saver{full_filename};
        saver.save_mesh(fe_mesh_->relation_table());
        saver.save_vector_data(velocity_field_.values_view(), table_name);
    }

    template<std::floating_point T>
    void KraichanMethodGaussianSpectraImpl<T>::generate_sample(value_type time) {
        const std::size_t nvert = fe_mesh_->n_vertices();

        for (const std::size_t ivert: std::views::iota(0ull, nvert)) {
            const auto vertex = fe_mesh_->relation_table()->vertex(ivert);
            auto fluct = generator_(vertex, time);
            velocity_field_.set_value(std::move(fluct), ivert);
        }
    }
}// namespace stg::spectral

#endif