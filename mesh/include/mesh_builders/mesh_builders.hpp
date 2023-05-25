#ifndef STG_MESH_BUILDERS_HPP
#define STG_MESH_BUILDERS_HPP

#include "cube_fe_mesh.hpp"
#include "fem/fe_factory.hpp"
#include "rtable/cube_relation_table.hpp"
#include "rtable/i_relation_table.hpp"
#include <concepts>
#include <execution>
#include <future>
#include <iostream>
#include <memory>
#include <range/v3/all.hpp>
#include <ranges>
#include <syncstream>
#include <thread>
#include <vector>

namespace stg::mesh {
    namespace rv = ranges::views;

    template<std::floating_point T>
    class CubePrizmMeshBuilder final {
    public:
        using value_type = T;

        CubePrizmMeshBuilder(T x_l, T x_r, T y_l, T y_r,
                             T z_l, T z_r, std::size_t nx,
                             std::size_t ny, std::size_t nz)
            : x_l_{x_l}, x_r_{x_r}, y_l_{y_l}, y_r_{y_r}, z_l_{z_l}, z_r_{z_r},
              nx_{nx}, ny_{ny}, nz_{nz} {}

        /*
     * Symmetric relative to (0, 0, 0) point
     */
        CubePrizmMeshBuilder(T Lx, T Ly, T Lz, std::size_t nx,
                             std::size_t ny, std::size_t nz)
            : x_l_{-Lx / 2.}, x_r_{Lx / 2.}, y_l_{-Ly / 2.},
              y_r_{Ly / 2.}, z_l_{-Lz / 2.}, z_r_{Lz / 2.},
              nx_{nx}, ny_{ny}, nz_{nz} {}

        [[nodiscard("Heavy object construction")]] std::shared_ptr<CubePrizmRelationTable<value_type>> build_relation_table() const {
            auto vert_future = std::async(std::launch::async, &CubePrizmMeshBuilder::assemble_vertices, this);
            auto bound_ind_future = std::async(std::launch::async, &CubePrizmMeshBuilder::assemble_bounds_indices, this);
            auto element_types_future = std::async(std::launch::async, &CubePrizmMeshBuilder::assemble_element_types, this);
            auto&& vertices = vert_future.get();
            auto&& bound_indices = bound_ind_future.get();
            auto&& element_types = element_types_future.get();

            return std::make_shared<CubePrizmRelationTable<value_type>>(
                    x_l_, x_r_, y_l_, y_r_, z_l_, z_r_, nx_, ny_, nz_,
                    std::move(vertices), std::move(bound_indices), std::move(element_types));
        }

        [[nodiscard("Heavy object construction")]] std::shared_ptr<CubePrizmFEMesh<value_type>> build() const {
            std::shared_ptr<CubePrizmRelationTable<value_type>> rtable = build_relation_table();

            auto&& fe_elements = assemble_voxel_elements(rtable);

            return std::make_shared<CubePrizmFEMesh<value_type>>(
                    x_l_, x_r_,
                    y_l_, y_r_,
                    z_l_, z_r_,
                    nx_, ny_, nz_,
                    std::move(rtable),
                    std::move(fe_elements));
        }

    protected:
        std::vector<value_type> assemble_vertices() const {
            std::vector<value_type> vertices{};
            vertices.reserve(size_ * 3);

            double z = z_l_;

            for (const size_t iz: rv::iota(0ul, nz_)) {
                double y = y_l_;
                for (const size_t iy: rv::iota(0ul, ny_)) {
                    double x = x_l_;
                    for (const size_t ix: rv::iota(0ul, nx_)) {
                        vertices.push_back(x);
                        vertices.push_back(y);
                        vertices.push_back(z);
                        x += dx_;
                    }
                    y += dy_;
                }
                z += dz_;
            }

            return vertices;
        }

        std::vector<std::vector<size_t>> assemble_bounds_indices() const {
            std::vector<std::vector<size_t>> bounds;
            bounds.reserve(size_);
            for (const size_t iz: rv::iota(0ul, nz_ - 1)) {
                for (const size_t iy: rv::iota(0ul, ny_ - 1)) {
                    for (const size_t ix: rv::iota(0ul, nx_ - 1)) {
                        const size_t base_index = iy * nx_ + ix + iz * (nx_ * ny_);
                        const size_t next_layer_index = base_index + nx_ * ny_;

                        bounds.push_back({base_index,
                                          base_index + 1,
                                          base_index + nx_ + 1,
                                          base_index + nx_,
                                          next_layer_index,
                                          next_layer_index + 1,
                                          next_layer_index + nx_ + 1,
                                          next_layer_index + nx_});
                    }
                }
            }

            return bounds;
        }

        std::vector<size_t> assemble_element_types() const {
            std::vector<size_t> elements(size_);
            std::fill(std::execution::par_unseq, elements.begin(), elements.end(), 12);
            return elements;
        }

        std::vector<std::shared_ptr<IFiniteElement<value_type>>>
        assemble_voxel_elements(const std::shared_ptr<CubePrizmRelationTable<value_type>>& rtable) const {
            const size_t nelem = rtable->n_elements();
            std::vector<std::shared_ptr<IFiniteElement<value_type>>> fe_elements(nelem);

            ranges::for_each(ranges::view::iota(0ul, nelem), [&fe_elements, &rtable, this](size_t index) {
                const auto& element_indices = rtable->element_vertices_indices(index);
                const auto vertices_2d = rtable->element_vertices(index);
                fe_elements[index] = fe_factory_(vertices_2d, element_indices, FiniteElementsFactory::VtkTypes::Voxel);
            });

            return fe_elements;
        }

        const T x_l_;
        const T x_r_;
        const T y_l_;
        const T y_r_;
        const T z_l_;
        const T z_r_;
        const std::size_t nx_;
        const std::size_t ny_;
        const std::size_t nz_;
        const T dx_ = (x_r_ - x_l_) / (nx_ - 1);
        const T dy_ = (y_r_ - y_l_) / (ny_ - 1);
        const T dz_ = (z_r_ - z_l_) / (nz_ - 1);
        const std::size_t size_ = (nx_ - 1) * (ny_ - 1) * (nz_ - 1);
        static const inline size_t vtk_cell_type_ = 12;
        const FiniteElementsFactory fe_factory_;
    };


    template<std::floating_point T>
    class CubeMeshBuilder final {
    public:
        using value_type = T;

        CubeMeshBuilder(T l, std::size_t n) : l_{l}, n_{n} {}

        [[nodiscard("Heavy object construction")]] std::shared_ptr<CubeRelationTable<value_type>> build_relation_table() const {
            auto vert_future = std::async(std::launch::async, &CubeMeshBuilder::assemble_vertices, this);
            auto bound_ind_future = std::async(std::launch::async, &CubeMeshBuilder::assemble_bounds_indices, this);
            auto element_types_future = std::async(std::launch::async, &CubeMeshBuilder::assemble_element_types, this);
            auto&& vertices = vert_future.get();
            auto&& bound_indices = bound_ind_future.get();
            auto&& element_types = element_types_future.get();

            return std::make_shared<CubeRelationTable<value_type>>(
                    l_, n_, std::move(vertices), std::move(bound_indices), std::move(element_types));
        }

        [[nodiscard("Heavy object construction")]] std::shared_ptr<CubeFiniteElementsMesh<value_type>> build() const {
            std::shared_ptr<CubeRelationTable<value_type>> rtable = build_relation_table();
            auto&& fe_elements = assemble_voxel_elements(rtable);

            return std::make_shared<CubeFiniteElementsMesh<value_type>>(
                    std::move(rtable),
                    std::move(fe_elements));
        }

    protected:
        std::vector<value_type> assemble_vertices() const {
            std::vector<value_type> vertices(n_);
            value_type left = -l_ / 2 - h_;
            std::generate(vertices.begin(), vertices.end(),
                          [&] { return left += h_; });
            return vertices;
        }

        std::vector<std::vector<size_t>> assemble_bounds_indices() const {
            std::vector<std::vector<size_t>> bounds;
            bounds.reserve(vertices_size_);
            for (const size_t iz: rv::iota(0ul, n_ - 1)) {
                for (const size_t iy: rv::iota(0ul, n_ - 1)) {
                    for (const size_t ix: rv::iota(0ul, n_ - 1)) {
                        const size_t base_index = iy * n_ + ix + iz * (n_ * n_);
                        const size_t next_layer_index = base_index + n_ * n_;

                        bounds.push_back({base_index, base_index + 1, base_index + n_ + 1, base_index + n_,
                                          next_layer_index, next_layer_index + 1, next_layer_index + n_ + 1, next_layer_index + n_});
                    }
                }
            }

            return bounds;
        }

        std::vector<size_t> assemble_element_types() const { return std::vector<size_t>(vertices_size_, 11); }

        std::vector<std::shared_ptr<IFiniteElement<value_type>>> assemble_voxel_elements(const std::shared_ptr<IRelationTable<value_type>>& rtable) const {
            const size_t nelem = rtable->n_elements();
            std::vector<std::shared_ptr<IFiniteElement<value_type>>> fe_elements(nelem);

            ranges::for_each(ranges::views::iota(0ul, nelem), [&fe_elements, &rtable, this](size_t index) {
                const auto& element_indices = rtable->element_vertices_indices(index);
                const auto vertices_2d = rtable->element_vertices(index);
                fe_elements[index] = fe_factory_(vertices_2d, element_indices, FiniteElementsFactory::VtkTypes::Voxel);
            });

            /*      ranges::views::iota(0ul, nelem) | ranges::for_each([&fe_elements, &rtable, this](size_t index) {
        const auto& element_indices = rtable->element_vertices_indices(index);
        const auto vertices_2d = rtable->element_vertices(index);
        fe_elements[index] = fe_factory_(vertices_2d, element_indices, FiniteElementsFactory::VtkTypes::Voxel);
      });*/

            return fe_elements;
        }

        const value_type l_;
        const std::size_t n_;
        const value_type h_ = (l_) / (n_ - 1);
        const std::size_t vertices_size_ = (n_ - 1) * (n_ - 1) * (n_ - 1);
        static const inline size_t vtk_cell_type_ = 12;
        const FiniteElementsFactory fe_factory_;
    };
}// namespace stg::mesh

#endif//STG_MESH_BUILDERS_HPP
