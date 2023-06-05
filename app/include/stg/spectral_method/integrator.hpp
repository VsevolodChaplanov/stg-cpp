#ifndef STG_APP_SPECTRAL_METHOD_INTEGRATOR_HPP
#define STG_APP_SPECTRAL_METHOD_INTEGRATOR_HPP

#include "fourier/ft.hpp"
#include "geometry/geometry.hpp"
#include "mesh_builders/cube_fe_mesh.hpp"
#include "stg/spectral_method/data_loader.hpp"
#include <catch2/internal/catch_clara.hpp>
#include <concepts>
#include <fourier.hpp>
#include <memory>
#include <ranges>
#include <spherical_mesh/sphere_mesh.hpp>
#include <string_view>
#include <utility>
#include <vector>
namespace stg::spectral {

    template<std::floating_point T>
    class IntegrateCovariations1D final {
    public:
        IntegrateCovariations1D(DataLoader loader, std::string_view cov_filename)
            : loader_(std::move(loader)), cov_mesh_(loader_.load_mesh<T>(cov_filename)), cov_values_{loader_.load_scalar_data<T>(cov_filename)} {
        }

        T integrate_for(T k_mod, std::size_t n_theta = 100, std::size_t n_phi = 100) const {
            stg::mesh::SphericalMesh<T> sphere_mesh_{k_mod, n_theta, n_phi};
            const auto k_wavevectors = sphere_mesh_.elements_centers();

            auto fert_values = k_wavevectors | std::views::transform([&](const Vector<T>& wave_vector) {
                                   return cov_mesh_->integrate_fourier(cov_values_ | std::views::all, wave_vector);
                               });

            auto result = sphere_mesh_.integrate(fert_values);

            return result;
        }

    private:
        DataLoader loader_;
        std::shared_ptr<mesh::CubeFiniteElementsMesh<T>> cov_mesh_;
        std::vector<T> cov_values_;
    };

    template<std::floating_point T>
    class IntegrateFert1D final {
    public:
        IntegrateFert1D(DataLoader loader, std::string_view cov_filename, std::string_view fert_mesh)
            : fert_values_filename_{fert_mesh},
              loader_(std::move(loader)),
              fert_mesh_{loader_.load_mesh<T>(fert_mesh)},
              cov_mesh_(loader_.load_mesh<T>(cov_filename)),
              cov_values_{loader_.load_scalar_data<T>(cov_filename)} {
        }

        void integrate_covariations() {
            const std::size_t n_vert = fert_mesh_->n_vertices();
            fert_values_.resize(n_vert);
            for (const std::size_t ivert: std::views::iota(0ul, n_vert)) {
                const auto wave_vector = fert_mesh_->relation_table()->vertex(ivert);
                const auto fert_value = cov_mesh_->integrate_fourier(cov_values_ | std::views::all, wave_vector);
                fert_values_[ivert] = fert_value;
            }
        }

        void save_calculated_fert(std::string_view filepath,
                                  std::string_view table_name = "TableName") const {
            VtkRectilinearGridSaver saver{filepath};
            saver.save_mesh(fert_mesh_->relation_table());
            saver.save_scalar_data(fert_values_.cbegin(),
                                   fert_values_.cend(),
                                   table_name);
        }

        std::pair<T, T> integrate_pure() {
            const auto fert_values_file = loader_.load_scalar_data<T>(fert_values_filename_);
            const auto fert_0 = fert_mesh_->integrate(fert_values_file);
            const auto center_lin_index = cov_mesh_->center_lin_index();
            const auto r_0 = cov_values_[center_lin_index];
            return {fert_0, r_0};
        }

        std::pair<T, T> apply_fourier_to_covariations() {
            const auto center_lin_index = cov_mesh_->center_lin_index();
            fert_values_ = kriging::fourier3(*cov_mesh_->relation_table(), cov_values_);
            const auto r_0 = cov_values_[center_lin_index];
            const auto fert_0 = fert_values_[center_lin_index];
            return {fert_0, r_0};
        }

    private:
        const std::string fert_values_filename_;
        DataLoader loader_;
        std::shared_ptr<mesh::CubeFiniteElementsMesh<T>> fert_mesh_;
        std::shared_ptr<mesh::CubeFiniteElementsMesh<T>> cov_mesh_;
        std::vector<T> cov_values_;
        std::vector<T> fert_values_;
    };
}// namespace stg::spectral


#endif