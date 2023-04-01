#ifndef STG_LAGRANGIAN_VOXEL_FE_HPP
#define STG_LAGRANGIAN_VOXEL_FE_HPP

#include <cmath>
#include "ifinite_element.hpp"

namespace stg::mesh {

  template<std::floating_point T>
  class VoxelFiniteElement final : public LagrangianFiniteElement<T, 8> {
    using BasisFunctions = std::array<std::function<double(Point<T>)>, 8>;
    using BasisFunctionsGradients = std::array<std::function<Vector<T>(Point<T>)>, 8>;
  public:
    using value_type = LagrangianFiniteElement<T, 8>::value_type;
//    using LagrangianFiniteElement<T, 8>::lumped;
    using LagrangianFiniteElement<T, 8>::integrate;
    using LagrangianFiniteElement<T, 8>::basis_functions_n;
    using LagrangianFiniteElement<T, 8>::global_indices;

    VoxelFiniteElement(std::vector<T> vertices, std::vector<size_t> global_indices)
       : VoxelFiniteElement(std::fabs(vertices[3] - vertices[0]),
                            std::fabs(vertices[10] - vertices[1]),
                            std::fabs(vertices[14] - vertices[2]),
                            {vertices[0], vertices[1], vertices[2]},
                            std::forward<decltype(global_indices)>(global_indices))
    { }

    /*
     * Converts param point to point of real space
     */
    Point<T> real_space_point(Point<T> param_point) const noexcept override {
      const auto semi_point = Point<T>{
        param_point.template get<0>() * dx_,
        param_point.template get<1>() * dy_,
        param_point.template get<2>() * dz_
      };
      auto result = base_point_ + semi_point;
      return result;
    }

    value_type jacobian(std::size_t i, std::size_t j) const override {
      if (i == 0 && j == 0) { return dx_; }
      if (i == 1 && j == 1) { return dy_; }
      if (i == 1 && j == 2) { return dz_; }
      return 0.;
    }

    value_type inverse_j(size_t i, size_t j) const override {
      if (i == 0 && j == 0) { return 1 / dx_; }
      if (i == 1 && j == 1) { return 1 / dy_; }
      if (i == 2 && j == 2) { return 1 / dz_; }
      return 0.;
    }

    value_type interpolate(const Point<T>& param_point, const std::vector<value_type>& values) const override {
      if (values.size() != 8) {
        throw std::logic_error("The array of values has the wrong size");
      }
      value_type result = 0;
      auto val_iter = values.cbegin();
      std::for_each(basis_functions_.cbegin(), basis_functions_.cend(), [&](const auto& basis_func) {
        value_type node_weight = basis_func(param_point);
        value_type value = *val_iter;
        result += node_weight * value;
        val_iter++;
      });
      return result;
    }

    ~VoxelFiniteElement() override = default;

  protected:
    VoxelFiniteElement(T dx, T dy, T dz,
                       Point<T> base_vertex,
                       std::vector<size_t> global_indices)
      : VoxelFiniteElement(dx, dy, dz,
                           dx * dy * dz,
                           std::forward<decltype(base_vertex)>(base_vertex),
                           std::forward<decltype(global_indices)>(global_indices))
    { }

    VoxelFiniteElement(T dx, T dy, T dz,
                       T det_j,
                       Point<T> base_vertex,
                       std::vector<size_t> global_indices)
      : VoxelFiniteElement(
          dx, dy, dz,
          det_j,
          {det_j / 8., det_j / 8., det_j / 8., det_j / 8.,
           det_j / 8., det_j / 8., det_j / 8., det_j / 8.},
          std::forward<decltype(base_vertex)>(base_vertex),
          std::forward<decltype(global_indices)>(global_indices))
    { }

    VoxelFiniteElement(T dx, T dy, T dz,
                       T det_j,
                       std::array<T, 8>&& lumped,
                       Point<T> base_vertex,
                       std::vector<size_t> global_indices)
      : LagrangianFiniteElement<T, 8>(
          std::forward<std::vector<size_t>>(global_indices),
          std::forward<std::array<T, 8>>(lumped),
          std::forward<Point<T>>(base_vertex)
        ),
        dx_{dx}, dy_{dy}, dz_{dz},
        det_j_{det_j}
    { }

    /*
     * Return values of gradients of basis function in given point
     */
    std::array<Vector<T>, 8> basis_gradients_at_point(Point<T> point) const override {
      std::array<Vector<T>, 8> result;
      for (size_t k = 0; const auto& grad_func : basis_functions_gradients_) {
        result[k] = grad_func(point);
        k++;
      }
      return result;
    }

    /*
     * Return array of array of gradients of basis functions
     * calculated in vertices
     */
    std::array<std::array<Vector<T>, 8>, 8> basis_gradients_at_vertices() const {
      return {
        basis_gradients_at_point({0, 0, 0}),
        basis_gradients_at_point({1, 0, 0}),
        basis_gradients_at_point({1, 1, 0}),
        basis_gradients_at_point({0, 1, 0}),
        basis_gradients_at_point({0, 0, 1}),
        basis_gradients_at_point({1, 0, 1}),
        basis_gradients_at_point({1, 1, 1}),
        basis_gradients_at_point({0, 1, 1}),
      };
    }

    const T dx_;
    const T dy_;
    const T dz_;
    const T det_j_;
    static inline const BasisFunctions basis_functions_{
      [](Point<T> point) { return (1 - point.template get<0>()) * (1 - point.template get<1>()) * (1 - point.template get<2>()); },
      [](Point<T> point) { return point.template get<0>() * (1 - point.template get<1>()) * (1 - point.template get<2>()); },
      [](Point<T> point) { return point.template get<0>() * point.template get<1>() * (1 - point.template get<2>()); },
      [](Point<T> point) { return (1 - point.template get<0>()) * point.template get<1>() * (1 - point.template get<2>()); },
      [](Point<T> point) { return (1 - point.template get<0>()) * (1 - point.template get<1>()) * point.template get<2>(); },
      [](Point<T> point) { return point.template get<0>() * (1 - point.template get<1>()) * point.template get<2>(); },
      [](Point<T> point) { return point.template get<0>() * point.template get<1>() * point.template get<2>(); },
      [](Point<T> point) { return (1 - point.template get<0>()) * point.template get<1>() * point.template get<2>(); },
    };
    static inline const BasisFunctionsGradients basis_functions_gradients_{
      [](Point<T> point) -> Vector<T> {
        return {
          - ((1 - point.template get<2>()) * (1 - point.template get<1>())),
          - ((1 - point.template get<2>()) * (1 - point.template get<0>())),
          - ((1 - point.template get<1>()) * (1 - point.template get<0>()))
        }; },
      [](Point<T> point) -> Vector<T> {
        return {
          ((1 - point.template get<2>()) * (1 - point.template get<1>())),
          - ((1 - point.template get<2>()) * point.template get<0>()),
          - ((1 - point.template get<1>()) * point.template get<0>())
        }; },
      [](Point<T> point) -> Vector<T> {
        return {
          (1 - point.template get<2>()) * point.template get<1>(),
          (1 - point.template get<2>()) * point.template get<0>(),
          (-point.template get<1>() * point.template get<0>())
        }; },
      [](Point<T> point) -> Vector<T> {
        return {
          - ((1 - point.template get<2>()) * point.template get<1>()),
          (1 - point.template get<2>()) * (1 - point.template get<0>  ()),
          - point.template get<1>() * (1 - point.template get<0>())
        }; },
      [](Point<T> point) -> Vector<T> {
        return {
          - point.template get<2>() * (1 - point.template get<1>()),
          - point.template get<2>() * (1 - point.template get<0>()),
          (1 - point.template get<1>()) * (1 - point.template get<0>())
        }; },
      [](Point<T> point) -> Vector<T> {
        return {
          point.template get<2>() * (1 - point.template get<1>()),
          -point.template get<0>() * point.template get<2>(),
          (1 - point.template get<1>()) * point.template get<0>()
        }; },
      [](Point<T> point) -> Vector<T> {
        return {
          point.template get<2>() * point.template get<1>(),
          point.template get<2>() * point.template get<0>(),
          point.template get<1>() * point.template get<0>()
        }; },
      [](Point<T> point) -> Vector<T> {
        return {
          - point.template get<2>() * point.template get<1>(),
          point.template get<2>() * (1 - point.template get<0>()),
          point.template get<1>() * (1 - point.template get<0>())
        }; },
    };
    using LagrangianFiniteElement<T, 8>::global_indices_;
    using LagrangianFiniteElement<T, 8>::lumped_;
    using LagrangianFiniteElement<T, 8>::base_point_;
  };
}

#endif //STG_LAGRANGIAN_VOXEL_FE_HPP
