#ifndef STG_STOCHASTIC_GAUSSIAN_HPP
#define STG_STOCHASTIC_GAUSSIAN_HPP

#include <memory>
#include <armadillo>
#include <rtable/cube_relation_table.hpp>

namespace stg::gaussian {

  template<std::floating_point T>
  class StochasticGaussian final {
  public:
    struct Params {
      std::size_t eigen_cut = 1000;
      double variance_cut = 0.05;
    };

    StochasticGaussian(std::shared_ptr<mesh::CubeRelationTable<T>> space)
  };
}

#endif //STG_STOCHASTIC_GAUSSIAN_HPP
