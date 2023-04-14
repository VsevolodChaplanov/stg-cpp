#include "common.hpp"

TEST_CASE("Tensor main usages tests") {
  GIVEN("Symmetric tensor") {
    constexpr double eps = 1.e-6;
    std::array<double, 9> init_array{1., 0.5, 0.2, 0.5, 1., 0.3, 0.2, 0.3, 1.};
    stg::tensor::Tensor<double> test_tensor{init_array};

    WHEN("Tensor is not changes it values are same as constructed from") {
      CHECK_THAT(test_tensor.get(0, 0), WithinRel(1., eps));
      CHECK_THAT(test_tensor.get(0, 2), WithinRel(0.2, eps));
      CHECK_THAT(test_tensor.get(1, 0), WithinRel(0.5, eps));
      CHECK_THAT(test_tensor.get(2, 1), WithinRel(0.3, eps));
    }

    WHEN("In tensor changed value its change") {
      test_tensor.set(0, 1, 0.4);
      test_tensor.set(1, 0, 0.4);
      CHECK_THAT(test_tensor.get(1, 0), WithinRel(0.4, eps));
      CHECK_THAT(test_tensor.get(0, 1), WithinRel(0.4, eps));
    }

    WHEN("Check cholesky decomposition") {
      auto lower_tensor = test_tensor.cholesky();
      CHECK_THAT(lower_tensor.get(0, 0), WithinRel(1., eps));
      CHECK_THAT(lower_tensor.get(1, 0), WithinRel(0.5, eps));
      CHECK_THAT(lower_tensor.get(1, 1), WithinRel(0.8660254, eps));
      CHECK_THAT(lower_tensor.get(2, 0), WithinRel(0.2, eps));
      CHECK_THAT(lower_tensor.get(2, 1), WithinRel(0.23094011, eps));
      CHECK_THAT(lower_tensor.get(2, 2), WithinRel(0.95219046, eps));

      CHECK_THAT(lower_tensor.get(0, 1), WithinAbs(0, eps));
    }
  }
}
