#ifndef STG_APP_SEQUENTIAL_ANALYSIS_HPP
#define STG_APP_SEQUENTIAL_ANALYSIS_HPP

#include "fourier/ft.hpp"
#include "geometry/geometry.hpp"
#include "mesh_builders/mesh_builders.hpp"
#include "spherical_mesh/sphere_mesh.hpp"
#include "statistics/space_correlation.hpp"
#include "stg/spectral_method/data_loader.hpp"
#include "stg_tensor/tensor.hpp"
#include "velocity_field/velocity_field.hpp"
#include <bits/ranges_base.h>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <ranges>
#include <string_view>
#include <vector>
namespace stg::spectral {

    template<std::floating_point T>
    class SequentialCorrelations final {
    public:
        struct WithoutGeneration {};

        using value_type = T;

        SequentialCorrelations(DataLoader loader, std::string_view filename)
            : loader_{std::move(loader)},
              velocity_mesh_{loader.load_mesh<value_type>(filename)},
              covariance_{covariance_mesh_->n_vertices()} {}

        SequentialCorrelations(DataLoader loader, std::string_view covariations_file, WithoutGeneration)
            : loader_{std::move(loader)},
              covariance_mesh_(loader_.load_mesh<value_type>(covariations_file)),
              covariance_{loader_.load_tensor_data<value_type>(covariations_file)} {}


        [[nodiscard]] std::size_t size() const { return covariance_mesh_->n_vertices(); }

        void calc_covariations_for_amount(std::size_t take_n_fields) {
            run_along_files(take_n_fields);
        }

        template<std::size_t i, std::size_t j>
        std::vector<value_type> calculate_energies(std::string_view filepath_for_fert,
                                                   std::string_view filepath_for_energ,
                                                   std::string_view table_name = "TableName") {
            const std::size_t nvert = covariance_mesh_->n_vertices();
            auto fourier_space = covariance_mesh_->relation_table()->template make_fourier_space<CubeMeshBuilder<value_type>>();
            std::vector<value_type> energies(nvert);

            std::vector<value_type> r_ii = covariance_ | std::views::transform([](const Tensor<value_type>& tensor) {
                                               return tensor.get(i, j);
                                           }) |
                                           ranges::to<std::vector<value_type>>();
            auto fert_values = kriging::fourier3(*covariance_mesh_->relation_table(), r_ii);
            VtkRectilinearGridSaver saver_fert{filepath_for_fert};
            saver_fert.save_mesh(fourier_space);
            saver_fert.save_scalar_data(fert_values.cbegin(),
                                        fert_values.cend(),
                                        table_name);

            for (const std::size_t index: std::views::iota(0ul, nvert)) {
                const auto k_vec = fourier_space->vertex(index);
                const auto fert = fert_values[index];
                const auto p_value = P<i, j>(k_vec);
                const auto k_vec_sqr = dot_product(k_vec, k_vec);
                const auto energy = fert * 4 * std::numbers::pi * k_vec_sqr / p_value;
                energies[index] = energy;
            }

            VtkRectilinearGridSaver saver{filepath_for_energ};
            saver.save_mesh(fourier_space);
            saver.save_scalar_data(energies.cbegin(),
                                   energies.cend(),
                                   table_name);

            return energies;
        }

        void save_calculated_covariations(std::string_view filepath,
                                          std::string_view table_name = "TableName") const {
            VtkRectilinearGridSaver saver{filepath};
            saver.save_mesh(covariance_mesh_->relation_table());
            saver.save_tensor_data(covariance_.cbegin(),
                                   covariance_.cend(),
                                   table_name);
        }

        void save_calculated_correlations(std::string_view filepath,
                                          std::string_view table_name = "TableName") const {
            VtkRectilinearGridSaver saver{filepath};
            saver.save_mesh(covariance_mesh_->relation_table());
            saver.save_tensor_data(correlations_.cbegin(),
                                   correlations_.cend(),
                                   table_name);
        }

        void convert_covariations_to_correlations() {
            correlations_.resize(covariance_mesh_->n_vertices());

            const std::size_t center_lin_index = covariance_mesh_->center_lin_index();
            const auto& center_tensor = covariance_[center_lin_index];
            const auto sigma_xx = std::sqrt(center_tensor.get(0, 0));
            const auto sigma_yy = std::sqrt(center_tensor.get(1, 1));
            const auto sigma_zz = std::sqrt(center_tensor.get(2, 2));

            for (const std::size_t ivert: std::views::iota(0ul, covariance_mesh_->n_vertices())) {
                const auto current_tensor = covariance_[ivert];
                correlations_[ivert] = Tensor<value_type>{
                        current_tensor.get(0, 0) / (sigma_xx * sigma_xx),
                        current_tensor.get(0, 1) / (sigma_xx * sigma_yy),
                        current_tensor.get(0, 2) / (sigma_xx * sigma_zz),
                        current_tensor.get(1, 0) / (sigma_yy * sigma_xx),
                        current_tensor.get(1, 1) / (sigma_yy * sigma_yy),
                        current_tensor.get(1, 2) / (sigma_yy * sigma_zz),
                        current_tensor.get(2, 0) / (sigma_zz * sigma_xx),
                        current_tensor.get(2, 1) / (sigma_zz * sigma_yy),
                        current_tensor.get(2, 2) / (sigma_zz * sigma_zz),
                };
            }
        }

        std::tuple<value_type, value_type, value_type> energy_for_wavenumber(value_type k_mod, std::size_t theta_n, std::size_t phi_n) const {
            SphericalMesh<value_type> sphere_mesh{k_mod, theta_n, phi_n};
            const auto k_points = sphere_mesh.elements_centers();
            auto r_xx = correlations_ | ranges::views::transform([](const Tensor<value_type>& corr) {
                            return corr.get(0, 0);
                        });
            auto r_yy = correlations_ | ranges::views::transform([](const Tensor<value_type>& corr) {
                            return corr.get(1, 1);
                        });
            auto r_zz = correlations_ | ranges::views::transform([](const Tensor<value_type>& corr) {
                            return corr.get(2, 2);
                        });

            std::vector<value_type> fert_xx;
            std::vector<value_type> fert_yy;
            std::vector<value_type> fert_zz;

            for (const auto& wave_vector: k_points) {
                const auto fert_integrated_xx = covariance_mesh_->integrate_fourier(r_xx, wave_vector);
                const auto fert_integrated_yy = covariance_mesh_->integrate_fourier(r_yy, wave_vector);
                const auto fert_integrated_zz = covariance_mesh_->integrate_fourier(r_zz, wave_vector);

                fert_xx.push_back(fert_integrated_xx);
                fert_yy.push_back(fert_integrated_yy);
                fert_zz.push_back(fert_integrated_zz);
            }

            const auto E_xx = sphere_mesh.integrate(fert_xx);
            const auto E_yy = sphere_mesh.integrate(fert_yy);
            const auto E_zz = sphere_mesh.integrate(fert_zz);

            return {E_xx, E_yy, E_zz};
        }

        template<std::size_t i, std::size_t j>
        std::vector<value_type> calculate_energies_by_symm_correlations(std::string_view filepath_for_fert,
                                                                        std::string_view filepath_for_energ,
                                                                        std::string_view table_name = "TableName") {
            const std::size_t nvert = covariance_mesh_->n_vertices();
            auto fourier_space = covariance_mesh_->relation_table()->template make_fourier_space<CubeMeshBuilder<value_type>>();
            std::vector<value_type> energies(nvert);

            std::vector<value_type> r_ii = symmetric_correlations_ | std::views::transform([](const Tensor<value_type>& tensor) {
                                               return tensor.get(i, j);
                                           }) |
                                           ranges::to<std::vector<value_type>>();
            auto fert_values = kriging::fourier3(*covariance_mesh_->relation_table(), r_ii);
            VtkRectilinearGridSaver saver_fert{filepath_for_fert};
            saver_fert.save_mesh(fourier_space);
            saver_fert.save_scalar_data(fert_values.cbegin(),
                                        fert_values.cend(),
                                        table_name);

            for (const std::size_t index: std::views::iota(0ul, nvert)) {
                const auto k_vec = fourier_space->vertex(index);
                const auto fert = fert_values[index];
                const auto p_value = P<i, j>(k_vec);
                const auto k_vec_sqr = dot_product(k_vec, k_vec);
                const auto energy = p_value < std::numeric_limits<double>::epsilon() ? 0 : fert * 4 * std::numbers::pi * k_vec_sqr / p_value;
                energies[index] = energy;
            }

            VtkRectilinearGridSaver saver{filepath_for_energ};
            saver.save_mesh(fourier_space);
            saver.save_scalar_data(energies.cbegin(),
                                   energies.cend(),
                                   table_name);

            return energies;
        }

        void make_correlations_symmetric() {
            const std::size_t nvert = covariance_mesh_->n_vertices();
            symmetric_correlations_.resize(nvert);
            for (std::size_t i = 0; i < nvert / 2 + 1; ++i) {
                const auto correlation = correlations_[i];
                symmetric_correlations_[i] = correlation;
                symmetric_correlations_[nvert - 1 - i] = correlation;
            }
        }

        void save_calculated_symm_correlations(std::string_view filepath,
                                               std::string_view table_name = "TableName") const {
            VtkRectilinearGridSaver saver{filepath};
            saver.save_mesh(covariance_mesh_->relation_table());
            saver.save_tensor_data(symmetric_correlations_.cbegin(),
                                   symmetric_correlations_.cend(),
                                   table_name);
        }

        template<std::size_t i, std::size_t j>
        std::vector<value_type> calculate_energies_by_symm_covariance(std::string_view filepath_for_fert,
                                                                      std::string_view filepath_for_energ,
                                                                      std::string_view table_name = "TableName") {
            const std::size_t nvert = covariance_mesh_->n_vertices();
            auto fourier_space = covariance_mesh_->relation_table()->template make_fourier_space<CubeMeshBuilder<value_type>>();
            std::vector<value_type> energies(nvert);

            std::vector<value_type> r_ii = symmetric_covariance_ | std::views::transform([](const Tensor<value_type>& tensor) {
                                               return tensor.get(i, j);
                                           }) |
                                           ranges::to<std::vector<value_type>>();
            auto fert_values = kriging::fourier3(*covariance_mesh_->relation_table(), r_ii);
            VtkRectilinearGridSaver saver_fert{filepath_for_fert};
            saver_fert.save_mesh(fourier_space);
            saver_fert.save_scalar_data(fert_values.cbegin(),
                                        fert_values.cend(),
                                        table_name);

            for (const std::size_t index: std::views::iota(0ul, nvert)) {
                const auto k_vec = fourier_space->vertex(index);
                const auto fert = fert_values[index];
                const auto p_value = P<i, j>(k_vec);
                const auto k_vec_sqr = dot_product(k_vec, k_vec);
                const auto energy = p_value < std::numeric_limits<double>::epsilon() ? 0 : fert * 4 * std::numbers::pi * k_vec_sqr / p_value;
                energies[index] = energy;
            }

            VtkRectilinearGridSaver saver{filepath_for_energ};
            saver.save_mesh(fourier_space);
            saver.save_scalar_data(energies.cbegin(),
                                   energies.cend(),
                                   table_name);

            return energies;
        }

        void make_covariance_symmetric() {
            const std::size_t nvert = covariance_mesh_->n_vertices();
            symmetric_covariance_.resize(nvert);
            for (std::size_t i = 0; i < nvert / 2 + 1; ++i) {
                const auto covariance = covariance_[i];
                symmetric_covariance_[i] = covariance;
                symmetric_covariance_[nvert - 1 - i] = covariance;
            }
        }

        void save_calculated_symm_covariance(std::string_view filepath,
                                             std::string_view table_name = "TableName") const {
            VtkRectilinearGridSaver saver{filepath};
            saver.save_mesh(covariance_mesh_->relation_table());
            saver.save_tensor_data(symmetric_covariance_.cbegin(),
                                   symmetric_covariance_.cend(),
                                   table_name);
        }

    private:
        DataLoader loader_;
        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> velocity_mesh_;
        const std::shared_ptr<CubeFiniteElementsMesh<value_type>> covariance_mesh_{velocity_mesh_};
        std::vector<Tensor<value_type>> covariance_;
        std::vector<Tensor<value_type>> correlations_;
        std::vector<Tensor<value_type>> symmetric_correlations_;
        std::vector<Tensor<value_type>> symmetric_covariance_;


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
                covariance_[i] = covariance_[i] + add_tensor;
                i++;
            }
        }

        void run_along_files(std::size_t amount) {
            for (const std::size_t c_try: rv::iota(0ull, amount)) {
                const std::string filename = fmt::format("velocity_field_{}.vtk", c_try);
                add_correlations_for_file(filename);
            }

            std::transform(covariance_.begin(), covariance_.end(),
                           covariance_.begin(), [amount](auto val) { return val / amount; });
        }

        Tensor<value_type> calculate_add(const Vector<value_type>& vertex_value, Vector<value_type> center_value) const {
            return {vertex_value.template get<0>() * center_value.template get<0>(), vertex_value.template get<0>() * center_value.template get<1>(), vertex_value.template get<0>() * center_value.template get<2>(),
                    vertex_value.template get<1>() * center_value.template get<0>(), vertex_value.template get<1>() * center_value.template get<1>(), vertex_value.template get<1>() * center_value.template get<2>(),
                    vertex_value.template get<2>() * center_value.template get<0>(), vertex_value.template get<2>() * center_value.template get<1>(), vertex_value.template get<2>() * center_value.template get<2>()};
        }

        template<std::size_t i, std::size_t j>
        constexpr T P(const Vector<T>& k_vec) const {
            const auto delta_ij = i == j ? 1 : 0;
            const auto k_i = k_vec.template get<i>();
            const auto k_j = k_vec.template get<j>();
            const auto k_sqr = dot_product(k_vec, k_vec);

            return delta_ij - (k_i * k_j) / k_sqr;
        }
    };
}// namespace stg::spectral

#endif