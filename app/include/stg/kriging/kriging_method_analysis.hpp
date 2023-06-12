#ifndef STG_KRIGING_METHOD_ANALYSIS_HPP
#define STG_KRIGING_METHOD_ANALYSIS_HPP

#include "../spectral_method/data_loader.hpp"
#include <iostream>
#include <statistics/space_covariation.hpp>
#include <velocity_field/velocity_samples.hpp>

using namespace stg::spectral;

namespace stg::kriging {
    using namespace stg::field;
    using namespace stg::statistics;

    template<std::floating_point T>
    class KrigingAnalysis final {
    public:
        using value_type = T;
        using tensors_iterator = std::vector<Tensor<value_type>>::iterator;

        explicit KrigingAnalysis(DataLoader loader)
            : real_space_mesh_{loader.load_mesh<T>("r_11.vtk")}, fourier_space_mesh_{loader.load_mesh<T>("phi_11.vtk")}, covariations_{loader.load_covariation_data<T>()}, fert_values_{loader.load_fert<T>()}, velocity_samples_{loader.load_velocity_samples<T>()}, calc_covariations_{real_space_mesh_->n_vertices()}, calc_fert_values_{real_space_mesh_->n_vertices()} {}

        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> real_space_mesh() const { return real_space_mesh_; }

        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> fourier_space_mesh() const { return fourier_space_mesh_; }

        tensors_iterator covariations_begin() { return covariations_.begin(); }
        tensors_iterator covariations_end() { return covariations_.end(); }
        tensors_iterator fert_begin() { return fert_values_.begin(); }
        tensors_iterator fert_end() { return fert_values_.end(); }

        auto covariations_view() { return ranges::views::all(covariations_); }
        auto fert_view() { return ranges::views::all(fert_values_); }

        const std::vector<Tensor<value_type>> calculate_covariations() {
            const auto lin_center_index = real_space_mesh_->center_lin_index();

            try {
                for (const size_t ivert: rv::iota(0ull, real_space_mesh_->n_vertices())) {
                    auto vertex_velocity_sample_x = velocity_samples_.vx_component_for_vertex(ivert);
                    auto vertex_velocity_sample_y = velocity_samples_.vy_component_for_vertex(ivert);
                    auto vertex_velocity_sample_z = velocity_samples_.vz_component_for_vertex(ivert);

                    auto center_velocity_sample_x = velocity_samples_.vx_component_for_vertex(lin_center_index);
                    auto center_velocity_sample_y = velocity_samples_.vy_component_for_vertex(lin_center_index);
                    auto center_velocity_sample_z = velocity_samples_.vz_component_for_vertex(lin_center_index);

                    // Assumes that fluctuations has zero mean
                    auto&& covariance = SpaceCovariance::covariance_tensor(center_velocity_sample_x, vertex_velocity_sample_x,
                                                                           center_velocity_sample_y, vertex_velocity_sample_y,
                                                                           center_velocity_sample_z, vertex_velocity_sample_z,
                                                                           0., 0., 0., 0., 0., 0.);
                    calc_covariations_[ivert] = std::move(covariance);
                }
            } catch (const std::exception& ex) {
                std::cerr << ex.what() << std::endl;
            }

            return calc_covariations_;
        }

        const std::vector<Tensor<value_type>> covariations_from_data() const { return covariations_; }

    private:
        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> real_space_mesh_;
        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> fourier_space_mesh_;
        const std::vector<Tensor<value_type>> covariations_;
        const std::vector<Tensor<value_type>> fert_values_;
        VelocitySamples<value_type> velocity_samples_;

        std::vector<Tensor<value_type>> calc_covariations_;
        std::vector<Tensor<value_type>> calc_fert_values_;
    };
}// namespace stg::kriging

#endif//STG_KRIGING_METHOD_ANALYSIS_HPP
