#ifndef STG_CUBE_RELATION_TABLE_HPP
#define STG_CUBE_RELATION_TABLE_HPP

#include "i_relation_table.hpp"
#include <array>
#include <range/v3/all.hpp>
#include <vector>

namespace stg::mesh {
    namespace rv = ranges::views;

    // template<std::floating_point T>
    // class CubeMeshBuilder;

    template<std::floating_point T>
    class CubePrizmRelationTable final : public IRelationTable<T> {
    public:
        using value_type = T;

        CubePrizmRelationTable(T x_l, T x_r, T y_l, T y_r, T z_l, T z_r,
                               std::size_t nx, std::size_t ny, std::size_t nz,
                               std::vector<value_type> vertices,
                               std::vector<std::vector<std::size_t>> bound_indices,
                               std::vector<std::size_t> elements_types)
            : x_l_{x_l}, x_r_{x_r}, y_l_{y_l},
              y_r_{y_r}, z_l_{z_l}, z_r_{z_r},
              nx_{nx}, ny_{ny}, nz_{nz},
              bound_indices_{std::move(bound_indices)},
              vertices_{std::move(vertices)}, element_types_{std::move(elements_types)} {}

        bool point_within(const Point<value_type>& point) const noexcept override {
            value_type x = point.template get<0>();
            value_type y = point.template get<1>();
            value_type z = point.template get<2>();

            if (x < x_l_ || x > x_r_) {
                return false;
            }

            if (y < y_l_ || y > y_r_) {
                return false;
            }

            if (z < z_l_ || z > z_r_) {
                return false;
            }

            return true;
        }

        [[nodiscard]] std::size_t lin_index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
            return i + j * nx_ + k * nx_ * ny_;
        }

        [[nodiscard]] std::array<std::size_t, 3> tri_index(std::size_t ivert) const noexcept {
            const size_t k = ivert / (nx_ * ny_);
            const size_t ij = ivert % (nx_ * ny_);
            const size_t j = ij / nx_;
            const size_t i = ij % nx_;

            return {i, j, k};
        }

        Point<T> vertex(std::size_t i, std::size_t j, std::size_t k) const {
            const auto lin_ind = lin_index(i, j, k);
            return vertex(lin_ind);
        }

        Point<T> vertex(std::size_t ivert) const override {
            return {vertices_[ivert * 3],
                    vertices_[ivert * 3 + 1],
                    vertices_[ivert * 3 + 2]};
        }

        [[nodiscard]] const std::vector<size_t>& element_vertices_indices(std::size_t ielem) const override {
            return bound_indices_[ielem];
        }

        [[nodiscard]] constexpr size_t element_type(std::size_t ielem) const override {
            return 12;
        }

        std::vector<Point<T>> element_vertices(std::size_t ielem) const override {
            const auto& bound_indices = bound_indices_[ielem];
            std::vector<Point<T>> result;
            std::for_each(bound_indices.cbegin(), bound_indices.cend(), [&result, this](auto index) {
                auto&& bound_vertex = this->vertex(index);
                result.push_back(bound_vertex);
            });
            return result;
        }

        const std::vector<value_type>& vertices() const noexcept override { return vertices_; }

        [[nodiscard]] const std::vector<std::vector<std::size_t>>& elements_global_indices() const noexcept override { return bound_indices_; }

        [[nodiscard]] const std::vector<std::size_t>& elements_types() const noexcept override { return element_types_; }

        [[nodiscard]] size_t n_vertices() const noexcept override { return n_vert_; }

        [[nodiscard]] size_t n_elements() const noexcept override { return n_elem_; }

        [[nodiscard]] size_t nx() const noexcept { return nx_; }

        [[nodiscard]] size_t ny() const noexcept { return ny_; }

        [[nodiscard]] size_t nz() const noexcept { return nz_; }

        value_type hx() const noexcept { return h_x_; }

        value_type hy() const noexcept { return h_y_; }

        value_type hz() const noexcept { return h_z_; }

        std::shared_ptr<CubePrizmRelationTable<value_type>> inverse_space() const noexcept {
            const value_type f_x_l_ = 2 * std::numbers::pi_v<value_type> / (x_r_ - x_l_) * nx_;
            const value_type f_y_l_ = 2 * std::numbers::pi_v<value_type> / (y_r_ - y_l_) * ny_;
            const value_type f_z_l_ = 2 * std::numbers::pi_v<value_type> / (z_r_ - z_l_) * nz_;

            CubePrizmRelationTable<value_type> cube_mesh_builder{f_x_l_, f_y_l_, f_z_l_,
                                                                 nx_, ny_, nz_};
            return cube_mesh_builder.build_relation_table();
        }

    private:
        const value_type x_l_;
        const value_type x_r_;
        const value_type y_l_;
        const value_type y_r_;
        const value_type z_l_;
        const value_type z_r_;
        const std::size_t nx_;
        const std::size_t ny_;
        const std::size_t nz_;
        const std::vector<std::vector<std::size_t>> bound_indices_;
        const std::vector<value_type> vertices_;
        const std::vector<std::size_t> element_types_;
        const std::size_t n_vert_ = nx_ * ny_ * nz_;
        const std::size_t n_elem_ = (nx_ - 1) * (ny_ - 1) * (nz_ - 1);
        const value_type l_x_ = x_r_ - x_l_;
        const value_type l_y_ = y_r_ - y_l_;
        const value_type l_z_ = z_r_ - z_l_;
        const value_type h_x_ = l_x_ / (nx_ - 1);
        const value_type h_y_ = l_y_ / (ny_ - 1);
        const value_type h_z_ = l_z_ / (nz_ - 1);
    };


    template<std::floating_point T>
    class CubeRelationTable final : public IRelationTable<T> {
    public:
        using value_type = T;

        CubeRelationTable(T L, std::size_t N,
                          std::vector<T> vertices,
                          std::vector<std::vector<std::size_t>> bound_indices,
                          std::vector<std::size_t> elements_types)
            : l_{L}, n_{N}, vertices_{std::move(vertices)}, bound_indices_{std::move(bound_indices)}, element_types_{std::move(elements_types)} {}

        bool point_within(const Point<value_type>& point) const noexcept override {
            value_type x = point.template get<0>();
            value_type y = point.template get<1>();
            value_type z = point.template get<2>();

            if (x < -l_ / 2 || x > l_ / 2) {
                return false;
            }

            if (y < -l_ / 2 || y > l_ / 2) {
                return false;
            }

            if (z < -l_ / 2 || z > l_ / 2) {
                return false;
            }

            return true;
        }

        Point<T> vertex(std::size_t ivert) const override {
            const auto [i, j, k] = tri_index(ivert);
            return {vertices_[i], vertices_[j], vertices_[k]};
        }

        [[nodiscard]] const std::vector<size_t>& element_vertices_indices(std::size_t ielem) const override {
            return bound_indices_[ielem];
        }

        [[nodiscard]] constexpr size_t element_type(std::size_t ielem) const override {
            return element_types_[ielem];
        }

        std::vector<Point<T>> element_vertices(std::size_t ielem) const override {
            const auto& bound_indices = bound_indices_[ielem];
            std::vector<Point<T>> result;
            std::for_each(bound_indices.cbegin(), bound_indices.cend(), [&result, this](auto index) {
                auto&& bound_vertex = this->vertex(index);
                result.push_back(bound_vertex);
            });
            return result;
        }

        const std::vector<value_type>& vertices() const noexcept override {
            return vertices_;
        }

        [[nodiscard]] const std::vector<std::vector<std::size_t>>& elements_global_indices() const noexcept override {
            return bound_indices_;
        }

        [[nodiscard]] const std::vector<size_t>& elements_types() const noexcept override {
            return element_types_;
        }

        [[nodiscard]] constexpr size_t n_vertices() const noexcept override {
            return n_ * n_ * n_;
        }

        [[nodiscard]] constexpr size_t n_elements() const noexcept override {
            return (n_ - 1) * (n_ - 1) * (n_ - 1);
        }

        // Return vertex of triplet index
        Point<T> vertex(std::size_t ix, std::size_t jy, std::size_t kz) const {
            return {vertices_[ix], vertices_[jy], vertices_[kz]};
        }

        // Convert triplet index on consecutive index
        [[nodiscard]] std::size_t lin_index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
            return i + j * n_ + k * n_ * n_;
        }

        // Convert consecutive index to tri index (ix, iy, iz)
        [[nodiscard]] std::array<std::size_t, 3> tri_index(std::size_t ivert) const noexcept {
            const size_t k = ivert / (n_ * n_);
            const size_t ij = ivert % (n_ * n_);
            const size_t j = ij / n_;
            const size_t i = ij % n_;

            return {i, j, k};
        }

        // Return discretization step h = l / (n - 1)
        constexpr value_type h() const noexcept { return h_; }

        // Return length of the edge, left coord is - l /2, right l / 2
        constexpr value_type l() const noexcept { return l_; }

        // Return discretization number along edge
        constexpr value_type n() const noexcept { return n_; }

        template<typename BuilderType>
        std::shared_ptr<CubeRelationTable<value_type>> make_fourier_space() const {
            const value_type hk = 2 * std::numbers::pi_v<value_type> / l_;
            const value_type lk = hk * n_;
            const BuilderType builder{lk, n_};
            return builder.build_relation_table();
        }

        template<typename BuilderType>
        std::shared_ptr<CubeRelationTable<value_type>> make_real_space() const {
            const value_type ls = 2 * std::numbers::pi_v<value_type> / h_;
            const BuilderType builder{ls, n_};
            return builder.build_relation_table();
        }

        [[nodiscard]] std::array<std::size_t, 3> center_tri_index() const {
            const auto center_n = static_cast<std::size_t>((n_ - 1) / 2);
            return {center_n, center_n, center_n};
        }

        std::size_t center_lin_index() const {
            const std::array<std::size_t, 3> tri_ind = center_tri_index();
            return lin_index(tri_ind[0], tri_ind[1], tri_ind[2]);
        }

        const std::vector<value_type>& linear_coordinates() const {
            return vertices_;
        }

    private:
        const value_type l_;                                       // length of cube's edges
        const std::size_t n_;                                      // number of points on edges
        const value_type h_ = l_ / (n_ - 1);                       // step between point on edges
        const std::vector<value_type> vertices_;                   // only one edge points {x0, x1, x2, ... }
        const std::vector<std::vector<std::size_t>> bound_indices_;// bound indices {{0, 1, 10, 11}, ... }
        const std::vector<std::size_t> element_types_;             // element types {12, 12, ... }
    };
}// namespace stg::mesh

#endif//STG_CUBE_RELATION_TABLE_HPP
