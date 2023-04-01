#ifndef STG_RELATION_TABLE_HPP
#define STG_RELATION_TABLE_HPP

#include <concepts>
#include <vector>
#include <geometry/geometry.hpp>
#include "i_relation_table.hpp"

namespace stg::mesh {

  template<std::floating_point T>
  class RelationTable final : public IRelationTable<T> {
  public:
    using value_type = T;

    RelationTable(std::vector<T>&& vertices,
                  std::vector<std::vector<std::size_t>>&& bound_indices,
                  std::vector<std::size_t>&& elements_types)
      : RelationTable(vertices.size() / 3, elements_types.size(),
                      std::move(vertices), std::move(bound_indices),
                      std::move(elements_types))
    { }

    Point<T> vertex(std::size_t ivert) const override {
      return {vertices_[ivert * 3],
              vertices_[ivert * 3 + 1],
              vertices_[ivert * 3 + 2]};
    }

    virtual Point<T> vertex(std::size_t i, std::size_t j, std::size_t k) const override {
      return {vertices_[i],
              vertices_[j],
              vertices_[k]};
    }

    const std::vector<size_t>& element_vertices_indices(std::size_t ielem) const override {
      return bound_indices_[ielem];
    }

    size_t element_type(std::size_t ielem) const override {
      return element_types_[ielem];
    }

    std::vector<Point<T>> element_vertices(std::size_t ielem) const override {
      const auto& bound_indices = bound_indices_[ielem];
      std::vector<Point<T>> result;
      for (const auto& vertex_index: bound_indices) {
        auto&& bound_vertex = vertex(vertex_index);
        result.push_back(bound_vertex);
      }
      return result;
    }

    const std::vector<value_type>& vertices() const override { return vertices_; }

    const std::vector<std::vector<std::size_t>>& elements_global_indices() const override { return bound_indices_; }

    const std::vector<std::size_t>& elements_types() const override { return element_types_; }

    std::shared_ptr<IRelationTable<value_type>> inverse_space() const override {
      std::vector<value_type> f_vertices;
      f_vertices.reserve(vertices_.size());
      ranges::for_each(vertices_, [&f_vertices](value_type coord) {
        f_vertices.push_back(2 * std::numbers::pi_v<value_type> / coord);
      });

      return std::make_shared<RelationTable<value_type>>(
        n_vert_, n_elem_, std::move(f_vertices),
        bound_indices_, element_types_
      );
    }

    size_t n_vertices() const override { return n_vert_; }

    size_t n_elements() const override { return n_elem_; }

    size_t nx() const override { return n_vert_ / 3; }

    size_t ny() const override { return n_vert_ / 3; }

    size_t nz() const override { return n_vert_ / 3; }

    size_t hx() const override {
      value_type min, max = vertices_ | ranges::views::take_fn([]() {
        static size_t index = 0;
        return index % 3 == 0;
      }) | ranges::minmax;
      return (max - min) / (n_vertices() / 3);
    }

    size_t hy() const override {
      value_type min, max = vertices_ | ranges::views::take_fn([]() {
        static size_t index = 0;
        return index + 1 % 3 == 0;
      }) | ranges::minmax;
      return (max - min) / (n_vertices() / 3); }

    size_t hz() const override {
      value_type min, max = vertices_ | ranges::views::take_fn([]() {
        static size_t index = 0;
        return index + 2 % 3 == 0;
      }) | ranges::minmax;
      return (max - min) / (n_vertices() / 3); }

  private:
    const std::size_t n_vert_;
    const std::size_t n_elem_;
    const std::vector<std::vector<std::size_t>> bound_indices_;
    const std::vector<value_type> vertices_;
    const std::vector<std::size_t> element_types_;

    RelationTable(std::size_t n_vert, std::size_t n_elem,
                  std::vector<value_type>&& vertices,
                  std::vector<std::vector<std::size_t>>&& bound_indices,
                  std::vector<std::size_t>&& elements_types)
      : n_vert_{n_vert}, n_elem_{n_elem}, bound_indices_{std::move(bound_indices)},
        vertices_{std::move(vertices)}, element_types_{std::move(elements_types)}
    { }

    RelationTable(std::size_t n_vert, std::size_t n_elem,
                  std::vector<value_type>&& vertices,
                  std::vector<std::vector<std::size_t>> bound_indices,
                  std::vector<std::size_t> elements_types)
      : n_vert_{n_vert}, n_elem_{n_elem}, bound_indices_{std::move(bound_indices)},
        vertices_{std::move(vertices)}, element_types_{std::move(elements_types)}
    { }
  };
}

#endif //STG_RELATION_TABLE_HPP
