#ifndef STG_KRIGING_1D_ANALYSIS_HPP
#define STG_KRIGING_1D_ANALYSIS_HPP

#include <velocity_field/velocity_field_1d.hpp>
#include <velocity_field/velocity_samples_1d.hpp>
#include <rtable/cube_vtk_saver.hpp>
#include <statistics/deviation.hpp>
#include <fmt/format.h>
#include "data_loader.hpp"


namespace stg::kriging {
using namespace stg::field;
using namespace stg::statistics;

  struct Params {
    std::size_t eigen_cut, N, L, vel_files_amount;
  }

  template<std::floating_point T>
  class SequentialAnalysis final {
  public:
    using value_type = T;

    explicit SequentialAnalysis(DataLoader loader, Params params)
    : loader_{std::move(loader_)}
    , params_{std::move(params)}
    , files_amount_{params_.vel_files_amount}
    , velocity_mesh(loader_.load_mesh(fmt::format(vel_file_string_format, params_.L, params_.N, params_.eigen_cut, 0)))
    , covariance_mesh(loader_.load_mesh("r_11.vtk")) { 
      calc_covariations_.resize(velocity_mesh_.n_vertices());
      std::fill(calc_covariations_.begin(), calc_covariations_.end(), 0);
    }

    value_type peaks_deviation() const {
        const auto center_index = real_space_mesh_->relation_table()->center_lin_index();
        const auto result = calc_covariations_[center_index] - covariations_[center_index];
        return result;
    }

    value_type integrate_square_diff() const {
      auto calc_cov_v = calc_covariations_ | rv::all;
      auto data_cov_v = covariations_ | rv::all;
      auto sqr_diff_view = Deviation::deviation_range_sqr(calc_cov_v, data_cov_v);
      return real_space_mesh_->integrate(sqr_diff_view);
    }

  protected:
    const std::string vel_file_string_format{"sg1_L{}_N{}_eigcut{}_try{}.vtk"}
    DataLoader loader_;
    const Params params_;
    const std::size_t files_amount_;
    const std::shared_ptr<CubeFiniteElementsMesh<value_type>> velocity_mesh_;
    const std::shared_ptr<CubeFiniteElementsMesh<value_type>> covariance_mesh_;

    std::vector<value_type> calc_covariations_;

    template<std::floating_point T>
    void covariances_for_file(std::string_view first_filename) {
      auto first_vel = loader_.load_scalar_data<T>(first_filename);
      covariances_for_samples(std::move(first_vel));
    }

    /*
     * Ковариации для двух повторений
     */
    template<ranges::viewable_range Range>
    void covariances_for_two_samples(Range&& first) {
      const auto center_lin_index = velocity_mesh_->center_lin_index();
      const auto center_val = first[center_lin_index];
      auto products = first | rv::transform([center_val] (const auto value) { return center_val * value; });
      for (const auto [index, value] : products | rv::enumerate) { 
        calc_covariations_[index] += value;
      }
    }

    void run_along_files() {
      for (const std::size_t c_try : rv::iota(0ull, files_amount_)) {
        const std::string filename = fmt::format(vel_file_string_format, params_.L, params_.N, params_.eigen_cut, c_try);
        covariances_for_file(filename);
      }
    }
  }
}

#endif