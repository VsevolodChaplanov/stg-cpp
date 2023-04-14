#ifndef STG_I_RELATION_TABLE_HPP
#define STG_I_RELATION_TABLE_HPP

#include <concepts>
#include <vector>
#include <geometry/geometry.hpp>

namespace stg::mesh {

  template<std::floating_point T>
  class IRelationTable {
  public:
    using value_type = T;

    // Check point belongs to domain
    virtual bool point_within(const Point<value_type>& point) const noexcept = 0;

    // Return vertex of the consecutive index
    virtual Point<T> vertex(std::size_t ivert) const = 0;

    // Return global indices to the given element
    virtual const std::vector<size_t>& element_vertices_indices(std::size_t ielem) const = 0;

    // Return type of given element
    virtual constexpr size_t element_type(std::size_t ielem) const = 0;

    // Returns nodes belonging to the given element
    virtual std::vector<Point<T>> element_vertices(std::size_t ielem) const = 0;

    // Return all vertices
    virtual const std::vector<value_type>& vertices() const noexcept = 0;

    // Return global indices of each element
    virtual const std::vector<std::vector<std::size_t>>& elements_global_indices() const noexcept = 0;

    // Return elements types in relation table
    virtual const std::vector<std::size_t>& elements_types() const noexcept = 0;

    // Return amount of all vertices in relation table
    virtual constexpr size_t n_vertices() const noexcept = 0;

    // Return amount of all element in relation table
    virtual constexpr size_t n_elements() const noexcept = 0;
  };
}


#endif //STG_I_RELATION_TABLE_HPP
