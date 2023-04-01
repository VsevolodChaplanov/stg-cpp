#ifndef STG_CUBE_FE_MESH_HPP
#define STG_CUBE_FE_MESH_HPP

#include <concepts>
#include <vector>
#include <fem.hpp>
#include <geometry/geometry.hpp>
#include "fem/i_fe_mesh.hpp"
#include "rtable/cube_relation_table.hpp"

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
        fe_elements_{std::move(elements)}
    { }

    std::size_t n_vertices() const override {
      return nx_ * ny_ * nz_;
    }

    std::size_t n_elements() const override {
      return size_;
    }

    const std::shared_ptr<IFiniteElement<T>>& element(size_t ielem) const override {
      return fe_elements_[ielem];
    }

    const std::vector<std::shared_ptr<IFiniteElement<T>>>& elements() const override {
      return fe_elements_;
    }

    std::shared_ptr<IRelationTable<T>> relation_table() const override { return cube_relation_table_; }

    constexpr std::size_t nx() const { return nx_; }

    constexpr std::size_t ny() const { return ny_; }

    constexpr std::size_t nz() const { return nz_; }

    std::array<std::size_t, 3> tri_index(std::size_t ivert) const {
      const size_t k = ivert / (nx_ * ny_);
      const size_t ij = ivert % (nx_ * ny_);
      const size_t j = ij / nx_;
      const size_t i = ij % nx_;

      return {i, j, k};
    }

    std::size_t lin_index(std::size_t ix, std::size_t jy, std::size_t kz) const {
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
  class CubeFiniteElementsMesh final : public  IFiniteElementsMesh<CubeRelationTable<T>> {
  public:
    using value_type = T;
    using TableType = typename IFiniteElementsMesh<CubeRelationTable<value_type>>::TableType;
    using ElementsIterator = typename std::vector<std::shared_ptr<IFiniteElement<value_type>>>::const_iterator;

    CubeFiniteElementsMesh(
      std::shared_ptr<CubeRelationTable<T>>&& cube_relation_table,
      std::vector<std::shared_ptr<IFiniteElement<T>>>&& fe_elements)
      : cube_relation_table_{std::forward<std::shared_ptr<CubeRelationTable<value_type>>>(cube_relation_table)}
      , fe_elements_{std::forward<std::vector<std::shared_ptr<IFiniteElement<value_type>>>>(fe_elements)}
    { }

    size_t n_vertices() const noexcept override {
      return cube_relation_table_->n_vertices();
    }

    size_t n_elements() const noexcept override {
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

    ranges::views::view_closure<ranges::views::all_fn>
    elements_view() const noexcept { return ranges::view::all(fe_elements_); }

    std::array<std::size_t, 3> center_tri_index() const {
      return cube_relation_table_->center_tri_index();
    }

    std::size_t center_lin_index() const {
      return cube_relation_table_->center_lin_index();
    }

  private:
    const std::shared_ptr<CubeRelationTable<value_type>> cube_relation_table_;
    const std::vector<std::shared_ptr<IFiniteElement<value_type>>> fe_elements_;
  };
}

#endif //STG_CUBE_FE_MESH_HPP
