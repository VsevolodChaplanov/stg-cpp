#ifndef STG_TRANSFORM_ACTIONS_HPP
#define STG_TRANSFORM_ACTIONS_HPP

#include <stg_tensor/tensor.hpp>

namespace stg::statistics {
using namespace stg::tensor;

  template<std::floating_point T>
  Tensor<T> normalize(const Tensor<T>& covariation) {
    const auto c_11 = covariation.template get(0, 0);
    const auto c_22 = covariation.template get(1, 1);
    const auto c_33 = covariation.template get(2, 2);

    const auto std_11 = std::sqrt(c_11);
    const auto std_22 = std::sqrt(c_22);
    const auto std_33 = std::sqrt(c_33);

    Tensor result {{
      covariation.template get(0, 0) / (std_11 * std_11),
          covariation.template get(0, 0) / (std_11 * std_22),
              covariation.template get(0, 0) / (std_11 * std_33),
      covariation.template get(0, 0) / (std_22 * std_11),
        covariation.template get(0, 0) / (std_22 * std_22),
          covariation.template get(0, 0) / (std_22 * std_33),
      covariation.template get(0, 0) / (std_33 * std_11),
        covariation.template get(0, 0) / (std_33 * std_22),
          covariation.template get(0, 0) / (std_33 * std_33),
    }};

    return result;
  }
}

#endif //STG_TRANSFORM_ACTIONS_HPP
