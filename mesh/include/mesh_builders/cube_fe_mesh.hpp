#ifndef STG_CUBE_FE_MESH_HPP
#define STG_CUBE_FE_MESH_HPP

#include "fem/i_fe_mesh.hpp"
#include "mesh_builder_fwd.hpp"
#include "rtable/cube_relation_table.hpp"
#include <cmath>
#include <complex>
#include <concepts>
#include <fem.hpp>
#include <geometry/geometry.hpp>
#include <memory>
#include <numbers>
#include <vector>

namespace stg::mesh {

    template<std::floating_point T>
    class CubePrizmFEMesh final : IFEMesh<T> {
    public:
        using IFEMesh<T>::integrate;

        CubePrizmFEMesh(T x_l, T x_r, T y_l, T y_r,
                        T z_l, T z_r, std::size_t nx,
                        std::size_t ny, std::size_t nz,
                        std::shared_ptr<CubePrizmRelationTable<T>>&& mesh,
                        std::vector<std::shared_ptr<IFiniteElement<T>>>&& elements)
            : x_l_{x_l}, x_r_{x_r}, y_l_{y_l},
              y_r_{y_r}, z_l_{z_l}, z_r_{z_r},
              nx_{nx}, ny_{ny}, nz_{nz},
              cube_relation_table_{std::move(mesh)},
              fe_elements_{std::move(elements)} {}

        [[nodiscard]] std::size_t n_vertices() const override {
            return nx_ * ny_ * nz_;
        }

        [[nodiscard]] std::size_t n_elements() const override {
            return size_;
        }

        const std::shared_ptr<IFiniteElement<T>>& element(std::size_t ielem) const override {
            return fe_elements_[ielem];
        }

        const std::vector<std::shared_ptr<IFiniteElement<T>>>& elements() const override {
            return fe_elements_;
        }

        std::shared_ptr<IRelationTable<T>> relation_table() const override { return cube_relation_table_; }

        [[nodiscard]] constexpr std::size_t nx() const { return nx_; }

        [[nodiscard]] constexpr std::size_t ny() const { return ny_; }

        [[nodiscard]] constexpr std::size_t nz() const { return nz_; }

        [[nodiscard]] std::array<std::size_t, 3> tri_index(std::size_t ivert) const {
            const size_t k = ivert / (nx_ * ny_);
            const size_t ij = ivert % (nx_ * ny_);
            const size_t j = ij / nx_;
            const size_t i = ij % nx_;

            return {i, j, k};
        }

        [[nodiscard]] std::size_t lin_index(std::size_t ix, std::size_t jy, std::size_t kz) const {
            return ix + jy * nx_ + kz * nx_ * ny_;
        }

    private:
        static const inline std::size_t vtk_cell_type = 12;
        const T x_l_;
        const T x_r_;
        const T y_l_;
        const T y_r_;
        const T z_l_;
        const T z_r_;
        const std::size_t nx_;
        const std::size_t ny_;
        const std::size_t nz_;
        const std::size_t size_ = (nx_ - 1) * (ny_ - 1) * (nz_ - 1);
        const std::shared_ptr<CubePrizmRelationTable<T>> cube_relation_table_;
        const std::vector<std::shared_ptr<IFiniteElement<T>>> fe_elements_;
    };

    template<std::floating_point T>
    class CubeFiniteElementsMesh final : public IFiniteElementsMesh<CubeRelationTable<T>> {
    public:
        using value_type = T;
        using TableType = typename IFiniteElementsMesh<CubeRelationTable<value_type>>::TableType;
        using ElementsIterator = typename std::vector<std::shared_ptr<IFiniteElement<value_type>>>::const_iterator;

        CubeFiniteElementsMesh(
                std::shared_ptr<CubeRelationTable<T>>&& cube_relation_table,
                std::vector<std::shared_ptr<IFiniteElement<T>>>&& fe_elements)
            : cube_relation_table_{std::forward<std::shared_ptr<CubeRelationTable<value_type>>>(cube_relation_table)}, fe_elements_{std::forward<std::vector<std::shared_ptr<IFiniteElement<value_type>>>>(fe_elements)} {}

        [[nodiscard]] size_t n_vertices() const noexcept override {
            return cube_relation_table_->n_vertices();
        }

        [[nodiscard]] size_t n_elements() const noexcept override {
            return cube_relation_table_->n_elements();
        }

        const std::shared_ptr<IFiniteElement<value_type>>& element(size_t ielem) const override {
            return fe_elements_[ielem];
        }

        const std::shared_ptr<TableType>& relation_table() const noexcept override {
            return cube_relation_table_;
        }

        ElementsIterator cbegin_elements() const noexcept { return fe_elements_.cbegin(); }

        ElementsIterator cend_elements() const noexcept { return fe_elements_.cend(); }

        auto elements_view() const noexcept { return ranges::views::all(fe_elements_); }

        [[nodiscard]] std::array<std::size_t, 3> center_tri_index() const {
            return cube_relation_table_->center_tri_index();
        }

        [[nodiscard]] std::size_t center_lin_index() const {
            return cube_relation_table_->center_lin_index();
        }

        value_type interpolate_at(const Point<value_type>& point, const std::vector<value_type>& values) const override {
            const auto element_base_indices = tri_index_from_point(point);
            const auto param_point_in_element = real_point_to_param(point);
            const std::size_t element_index = cube_relation_table_->lin_index(
                    element_base_indices[0],
                    element_base_indices[1],
                    element_base_indices[2]);
            const auto& element = fe_elements_[element_index];
            if (static_cast<std::size_t>(ranges::distance(values)) != element->basis_functions_n()) {
                throw std::logic_error("Range values not same with amount of element vertices");
            }
            std::vector<value_type> values_converted{ranges::begin(values), ranges::end(values)};
            auto result = element->interpolate(param_point_in_element, values_converted);
            return result;
        }

        template<std::ranges::viewable_range Range>
        value_type integrate(Range&& values) const {
            value_type result = 0.;
            for (auto elem: fe_elements_) {
                const auto indices = elem->global_indices();
                std::vector<value_type> element_values(elem->basis_functions_n());
                for (const auto ind: rv::iota(0ull, elem->basis_functions_n())) {
                    element_values[ind] = values[indices[ind]];
                }
                const auto elem_integral = elem->integrate(element_values);
                result += elem_integral;
            }

            return result;
        }

        template<std::ranges::viewable_range Range>
        value_type integrate_fourier(Range&& values, const Vector<value_type>& wave_vector) const {
            value_type fourier_coeff = 1 / (2 * std::numbers::pi_v<value_type>);
            fourier_coeff = fourier_coeff * fourier_coeff * fourier_coeff;
            value_type result = 0.;
            for (auto elem: fe_elements_) {
                const auto indices = elem->global_indices();
                std::vector<value_type> element_values(elem->basis_functions_n());
                for (const auto ind: rv::iota(0ull, elem->basis_functions_n())) {
                    const auto value = values[indices[ind]];
                    // e ^ (- i k r) = e ^ (phase)
                    const auto phase = dot_product(wave_vector, cube_relation_table_->vertex(ind));
                    const auto to_integrate = value * std::cos(phase) - value * std::sin(phase);//
                    element_values[ind] = to_integrate / fourier_coeff;
                }
                const auto elem_integral = elem->integrate(element_values);
                result += elem_integral;
            }

            return result;
        }

        template<std::ranges::viewable_range Range>
        value_type integrate_inverse_fourier(Range&& values, const Vector<value_type>& wave_vector) const {
            value_type result = 0.;
            for (auto elem: fe_elements_) {
                const auto indices = elem->global_indices();
                std::vector<value_type> element_values(elem->basis_functions_n());
                for (const auto ind: rv::iota(0ull, elem->basis_functions_n())) {
                    const auto value = values[indices[ind]];
                    // e ^ (- i k r) = e ^ (phase)
                    const auto phase = dot_product(wave_vector, cube_relation_table_->vertex(ind));
                    const auto to_integrate = value * (std::cos(phase) + std::sin(phase));//
                    element_values[ind] = to_integrate;
                }
                const auto elem_integral = elem->integrate(element_values);
                result += elem_integral;
            }

            return result;
        }

        template<std::ranges::viewable_range Range>
        std::complex<value_type> integrate_fourier2(Range&& values, const Vector<value_type>& wave_vector) const {
            value_type fourier_coeff = 1 / (2 * std::numbers::pi_v<value_type>);
            fourier_coeff = 1;//fourier_coeff * fourier_coeff * fourier_coeff;
            std::complex<value_type> result = 0.;
            for (auto elem: fe_elements_) {
                const auto indices = elem->global_indices();
                std::vector<std::complex<value_type>> element_values(elem->basis_functions_n());
                for (const auto ind: rv::iota(0ull, elem->basis_functions_n())) {
                    const auto value = values[indices[ind]];
                    // e ^ (- i k r) = e ^ (phase)
                    const auto phase = dot_product(wave_vector, cube_relation_table_->vertex(ind));
                    const auto to_integrate = std::complex<value_type>{value * std::cos(phase), -value * std::sin(phase)};
                    element_values[ind] = to_integrate / fourier_coeff;
                }
                const auto elem_integral = elem->integrate(element_values);
                result += elem_integral;
            }

            return result;
        }

        template<std::ranges::viewable_range Range>
        std::complex<value_type> integrate_inverse_fourier2(Range&& values, const Vector<value_type>& wave_vector) const {
            std::complex<value_type> result = 0.;
            for (auto elem: fe_elements_) {
                const auto indices = elem->global_indices();
                std::vector<std::complex<value_type>> element_values(elem->basis_functions_n());
                for (const auto ind: rv::iota(0ull, elem->basis_functions_n())) {
                    const auto value = values[indices[ind]];
                    // e ^ (- i k r) = e ^ (phase)
                    const auto phase = dot_product(wave_vector, cube_relation_table_->vertex(ind));
                    const auto value_to_integrate = value * std::complex<value_type>{std::cos(phase), std::sin(phase)};
                    // const auto to_integrate = std::complex<value_type>{value * std::cos(phase), +value * std::sin(phase)};
                    element_values[ind] = value_to_integrate;
                }
                const auto elem_integral = elem->integrate(element_values);
                result += elem_integral;
            }

            return result;
        }

    private:
        const std::shared_ptr<CubeRelationTable<value_type>> cube_relation_table_;
        const std::vector<std::shared_ptr<IFiniteElement<value_type>>> fe_elements_;

        std::array<std::size_t, 3> tri_index_from_point(const Point<value_type>& point) const {
            const auto base_point = -cube_relation_table_->l() / 2;

            const auto i_float_index = (point.template get<0>() - base_point) / cube_relation_table_->h();
            const auto j_float_index = (point.template get<1>() - base_point) / cube_relation_table_->h();
            const auto k_float_index = (point.template get<2>() - base_point) / cube_relation_table_->h();

            const auto i = static_cast<std::size_t>(std::floor(i_float_index));
            const auto j = static_cast<std::size_t>(std::floor(j_float_index));
            const auto k = static_cast<std::size_t>(std::floor(k_float_index));

            return {i, j, k};
        }

        Point<value_type> real_point_to_param(const Point<value_type>& point) const {
            const auto base_point = -cube_relation_table_->l() / 2;

            const auto i_float_index = (point.template get<0>() - base_point) / cube_relation_table_->h();
            const auto j_float_index = (point.template get<1>() - base_point) / cube_relation_table_->h();
            const auto k_float_index = (point.template get<2>() - base_point) / cube_relation_table_->h();

            const auto i = static_cast<std::size_t>(std::floor(i_float_index));
            const auto j = static_cast<std::size_t>(std::floor(j_float_index));
            const auto k = static_cast<std::size_t>(std::floor(k_float_index));

            return {i_float_index - i, j_float_index - j, k_float_index - k};
        }
    };
}// namespace stg::mesh

#endif//STG_CUBE_FE_MESH_HPP
