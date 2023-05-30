#ifndef STG_APP_SEQUENTIAL_ANALYSIS_HPP
#define STG_APP_SEQUENTIAL_ANALYSIS_HPP

#include "geometry/geometry.hpp"
#include "statistics/space_correlation.hpp"
#include "stg/spectral_method/data_loader.hpp"
#include "stg_tensor/tensor.hpp"
#include "velocity_field/velocity_field.hpp"
#include <bits/ranges_base.h>
#include <concepts>
#include <cstddef>
#include <functional>
#include <range/v3/view/transform.hpp>
#include <string_view>
#include <vector>
namespace stg::spectral {

    template<std::floating_point T>
    class SequentialCorrelations final {
    public:
        using value_type = T;
        SequentialCorrelations(DataLoader loader, std::string_view filename)
            : loader_{std::move(loader)},
              velocity_mesh_{loader.load_mesh<value_type>(filename)},
              correlations_{covariance_mesh_->n_vertices()} {}


        void calc_covariations_for_amount(std::size_t take_n_fields) {
            run_along_files(take_n_fields);
        }

        void save_calculated_covariations(std::string_view filepath,
                                          std::string_view table_name = "TableName") const {
            VtkRectilinearGridSaver saver{filepath};
            saver.save_mesh(covariance_mesh_->relation_table());
            saver.save_tensor_data(correlations_.cbegin(),
                                   correlations_.cend(),
                                   table_name);
        }

    private:
        DataLoader loader_;
        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> velocity_mesh_;
        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> covariance_mesh_{velocity_mesh_};
        std::vector<Tensor<value_type>> correlations_;

        void add_correlations_for_file(std::string_view filename) {
            const auto velocity_field = loader_.load_velocity_field<value_type>(filename);
            add_correlations(velocity_field);
        }

        void add_correlations(const VelocityField<value_type>& values) {
            const auto center_lin_index = velocity_mesh_->center_lin_index();
            const auto center_val = values.value(center_lin_index);

            auto tensors = values.values_view() | ranges::views::transform([&center_val, this](const auto& velocity_tuple) {
                               const auto [x, y, z] = velocity_tuple;
                               return calculate_add(Vector<value_type>{x, y, z}, center_val);
                           });
            for (std::size_t i = 0; const auto add_tensor: tensors) {
                correlations_[i] = correlations_[i] + add_tensor;
                i++;
            }
        }

        void run_along_files(std::size_t amount) {
            for (const std::size_t c_try: rv::iota(0ull, amount)) {
                const std::string filename = fmt::format("velocity_field_{}.vtk", c_try);
                add_correlations_for_file(filename);
            }

            std::transform(correlations_.begin(), correlations_.end(),
                           correlations_.begin(), [amount](auto val) { return val / amount; });
        }

        Tensor<value_type> calculate_add(const Vector<value_type>& vertex_value, Vector<value_type> center_value) const {
            return {vertex_value.template get<0>() * center_value.template get<0>(), vertex_value.template get<0>() * center_value.template get<1>(), vertex_value.template get<0>() * center_value.template get<2>(),
                    vertex_value.template get<1>() * center_value.template get<0>(), vertex_value.template get<1>() * center_value.template get<1>(), vertex_value.template get<1>() * center_value.template get<2>(),
                    vertex_value.template get<2>() * center_value.template get<0>(), vertex_value.template get<2>() * center_value.template get<1>(), vertex_value.template get<2>() * center_value.template get<2>()};
        }
    };
}// namespace stg::spectral

#endif