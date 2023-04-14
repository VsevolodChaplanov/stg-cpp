#ifndef STG_FE_MESH_HPP
#define STG_FE_MESH_HPP

#include <iterator>
#include <memory>
#include <ranges>
#include <range/v3/range_concepts.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/all.hpp>
#include <rtable/relation_table.hpp>
#include "ifinite_element.hpp"
#include "i_fe_mesh.hpp"

namespace stg::mesh {

  template<std::floating_point T>
  class FEMesh final : public IFEMesh<T> {
  public:
    FEMesh(std::shared_ptr<IRelationTable<T>>&& table,
           std::vector<std::shared_ptr<IFiniteElement<T>>>&& elements)
           : rtable_{std::move(table)},
             elements_{std::move(elements)}
    { }

    std::size_t n_vertices() const override { return rtable_.n_vertices(); }

    std::size_t n_elements() const override { return elements_.size(); }

    const std::shared_ptr<IFiniteElement<T>>& element(size_t ielem) const override { return elements_[ielem]; }

    const std::vector<std::shared_ptr<IFiniteElement<T>>>& elements() const override { return elements_; }

    std::shared_ptr<IRelationTable<T>> relation_table() const override { return rtable_; }

  protected:
    const std::shared_ptr<IRelationTable<T>> rtable_;
    const std::vector<std::shared_ptr<IFiniteElement<T>>> elements_;
  };
}

#endif //STG_FE_MESH_HPP
