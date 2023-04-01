#ifndef STG_I_FE_MESH_HPP
#define STG_I_FE_MESH_HPP

#include <iterator>
#include <memory>
#include <ranges>
#include <range/v3/range_concepts.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/all.hpp>
#include <rtable/relation_table.hpp>
#include "ifinite_element.hpp"

namespace stg::mesh {

  template<std::floating_point T>
  class IFEMesh {
  public:
    virtual std::size_t n_vertices() const = 0;

    virtual std::size_t n_elements() const = 0;

    virtual const std::shared_ptr<IFiniteElement<T>>& element(size_t ielem) const = 0;

    virtual const std::vector<std::shared_ptr<IFiniteElement<T>>>& elements() const = 0;

    virtual std::shared_ptr<IRelationTable<T>> relation_table() const = 0;

    template<ranges::bidirectional_range Range>
    T integrate(const Range& range) const {
      const auto integrate_lambda = [this, &range](std::size_t index) {
        const auto& fe_element = element(index);
        const auto& element_g_indices = fe_element->global_indices();
        std::vector<T> values_to_integrate;
        values_to_integrate.reserve(element_g_indices.size());
        for (const auto global_index: element_g_indices) {
          values_to_integrate.push_back(range[global_index]);
        }
        auto element_integration_result = fe_element->integrate(values_to_integrate);
        return element_integration_result;
      };
      double result = 0.;

      for (const size_t index : ranges::views::iota(0ul, n_elements())) {
        result += integrate_lambda(index);
      }
      return result;
    }
  };

  template<typename RelationTable>
  class IFiniteElementsMesh {
  public:
    using TableType = RelationTable;
    using value_type = typename RelationTable::value_type;

    virtual std::size_t n_vertices() const noexcept = 0;

    virtual std::size_t n_elements() const noexcept = 0;

    virtual const std::shared_ptr<IFiniteElement<value_type>>& element(size_t ielem) const = 0;

    virtual const std::shared_ptr<TableType>& relation_table() const noexcept = 0;

    template<ranges::bidirectional_range Range>
    value_type integrate(const Range& range) const {
      const auto integrate_lambda = [this, &range](std::size_t index) {
        const auto& fe_element = element(index);
        const auto& element_g_indices = fe_element->global_indices();
        std::vector<value_type> values_to_integrate;
        values_to_integrate.reserve(element_g_indices.size());
        for (const auto global_index: element_g_indices) {
          values_to_integrate.push_back(range[global_index]);
        }
        auto element_integration_result = fe_element->integrate(values_to_integrate);
        return element_integration_result;
      };
      double result = 0.;

      for (const size_t index : ranges::views::iota(0ul, n_elements())) {
        result += integrate_lambda(index);
      }
      return result;
    }
  };
}

#endif //STG_I_FE_MESH_HPP
