#ifndef STG_KRIGING_1D_ANALYSIS_HPP
#define STG_KRIGING_1D_ANALYSIS_HPP

#include <velocity_field/velocity_field_1d.hpp>
#include <velocity_field/velocity_samples_1d.hpp>
#include "data_loader.hpp"

namespace stg::kriging {
using namespace stg::field;
using namespace stg::statistics;

  template<std::floating_point T>
  class KrigingAnalysis1D final {
  public:
    using value_type = T;
    using tensors_iterator = std::vector<value_type>::iterator;

    explicit KrigingAnalysis1D(DataLoader loader)
    : real_space_mesh_{loader.load_mesh<T>("r_11.vtk")}
    , fourier_space_mesh_{loader.load_mesh<T>("phi_11.vtk")}
    , covariations_{loader.load_scalar_data<T>("r_11.vtk")}
    , fert_values_{loader.load_scalar_data<T>("phi_11.vtk")}
    , velocity_samples_{loader.load_velocity_samples_1sg<T>()} {
      calc_covariations_.resize(real_space_mesh_->n_vertices());
      calc_fert_values_.resize(real_space_mesh_->n_vertices());
    }

    const std::shared_ptr<CubeFiniteElementsMesh<value_type>> real_space_mesh() const { return real_space_mesh_; }

    const std::shared_ptr<CubeFiniteElementsMesh<value_type>> fourier_space_mesh() const { return fourier_space_mesh_; }

    tensors_iterator covariations_begin() { return covariations_.begin(); }
    tensors_iterator covariations_end() { return covariations_.end(); }
    tensors_iterator fert_begin() { return fert_values_.begin(); }
    tensors_iterator fert_end() { return fert_values_.end(); }
    const std::vector<value_type> covariations_from_data() const { return covariations_; }

    const std::vector<value_type> calculate_covariations() {
      const auto center_vert_lin_index = real_space_mesh_->center_lin_index();
      const auto center_vert = real_space_mesh_->relation_table()->vertex(center_vert_lin_index);

      auto center_velocity_sample_x = velocity_samples_.vx_component_for_vertex(center_vert_lin_index);

      try {
        for (const std::size_t ivert: rv::iota(0ull, real_space_mesh_->n_vertices())) {
          auto vertex_vel_sample = velocity_samples_.vx_component_for_vertex(ivert);
          auto covariance = Covariance::covariance(center_velocity_sample_x, vertex_vel_sample, 0., 0.);
          calc_covariations_[ivert] = covariance;
        }
      } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
      }

      return calc_covariations_;
    }

    value_type std_btw_covariations() const {
      auto cal_view = calc_covariations_ | rv::all;
      auto data_view = covariations_ | rv::all;
      return StandardDeviation::std(cal_view, data_view);
    }

  private:
    const std::shared_ptr<CubeFiniteElementsMesh<value_type>> real_space_mesh_;
    const std::shared_ptr<CubeFiniteElementsMesh<value_type>> fourier_space_mesh_;
    const std::vector<value_type> covariations_;
    const std::vector<value_type> fert_values_;
    VelocitySamples1D<value_type> velocity_samples_;

    std::vector<value_type> calc_covariations_;
    std::vector<value_type> calc_fert_values_;
  };

}

#endif //STG_KRIGING_1D_ANALYSIS_HPP
