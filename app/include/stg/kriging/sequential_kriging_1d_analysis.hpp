#ifndef STG_SEQUENTIAL_KRIGING_1D_ANALYSIS_HPP
#define STG_SEQUENTIAL_KRIGING_1D_ANALYSIS_HPP

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
  };

  template<std::floating_point T>
  class SequentialAnalysis final {
  public:
    using value_type = T;

    explicit SequentialAnalysis(DataLoader loader, Params params)
    : loader_{std::move(loader)}
    , params_{std::move(params)}
    , velocity_mesh_(loader_.load_mesh<value_type>(
        fmt::format("sg1_L{}_N{}_eigcut{}_try{}.vtk", params_.L, params_.N, params_.eigen_cut, 0)
      ))
    , covariance_mesh_(loader_.load_mesh<value_type>("r_11.vtk"))
    , covariations_{loader_.load_scalar_data<value_type>("r_11.vtk")} {
      calc_covariations_.resize(velocity_mesh_->n_vertices());
      std::fill(calc_covariations_.begin(), calc_covariations_.end(), 0);
    }

    value_type peaks_deviation() const {
        const auto center_index = covariance_mesh_->relation_table()->center_lin_index();
        const auto result = calc_covariations_[center_index] - covariations_[center_index];
        return result;
    }

    value_type integrate_square_diff() const {
      auto calc_cov_v = calc_covariations_ | rv::all;
      auto data_cov_v = covariations_ | rv::all;
      auto sqr_diff_view = Deviation::deviation_range_sqr(calc_cov_v, data_cov_v);
      return covariance_mesh_->integrate(sqr_diff_view);
    }

    void calc_covariations_for_amount(std::size_t take_n_fields) {
      run_along_files(take_n_fields);
    }

    void save_calculated_covariations(std::string_view filepath,
                                      std::string_view table_name = "TableName") const {
      VtkRectilinearGridSaver saver{filepath};
      saver.save_mesh(velocity_mesh_->relation_table());
      saver.save_scalar_data(calc_covariations_.cbegin(),
                             calc_covariations_.cend(),
                             table_name);
    }

  protected:
    DataLoader loader_;
    const Params params_;
    const std::shared_ptr<CubeFiniteElementsMesh<value_type>> velocity_mesh_;
    const std::shared_ptr<CubeFiniteElementsMesh<value_type>> covariance_mesh_;

    std::vector<value_type> calc_covariations_;
    std::vector<value_type> covariations_;

    void covariances_for_file(std::string_view first_filename) {
      auto first_vel = loader_.load_scalar_data<value_type>(first_filename);
      covariances_for_two_samples(first_vel | rv::all);
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

    void run_along_files(std::size_t amount) {for (const std::size_t c_try : rv::iota(0ull, amount)) {
        const std::string filename = fmt::format("sg1_L{}_N{}_eigcut{}_try{}.vtk", params_.L, params_.N, params_.eigen_cut, c_try);
        covariances_for_file(filename);
      }

      std::transform(calc_covariations_.begin(), calc_covariations_.end(),
                     calc_covariations_.begin(), [amount](auto val) { return val / amount; });
    }
  };
}

#endif