#ifndef STG_FE_FACTORY_HPP
#define STG_FE_FACTORY_HPP

#include "lagrangian_voxel_fe.hpp"

namespace stg::mesh {

  class FiniteElementsFactory {
  public:
    enum struct VtkTypes {
      Voxel = 12,
    };

    FiniteElementsFactory() = default;

    template<std::floating_point T>
    std::shared_ptr<IFiniteElement<T>> operator()(
      const std::vector<Point<T>>& vertices,
      const std::vector<std::size_t>& global_indices,
      VtkTypes type) const {

      std::vector<T> vertices_1d{};
      vertices_1d.reserve(vertices.size() * 3);

      std::for_each(vertices.cbegin(), vertices.cend(),
        [&vertices_1d](const auto& point) {
          vertices_1d.push_back(point.template get<0>());
          vertices_1d.push_back(point.template get<1>());
          vertices_1d.push_back(point.template get<2>());
        }
      );

      switch (type) {
        case VtkTypes::Voxel:
          return std::make_shared<VoxelFiniteElement<T>>(vertices_1d, global_indices);
        default:
          throw std::logic_error("No such element type");
      }
    }

  };
}

#endif //STG_FE_FACTORY_HPP
