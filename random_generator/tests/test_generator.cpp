#include "common.hpp"

using namespace stg;

struct GeneratorFixture {
  constexpr static inline double eps = 1.e-6;
  std::mt19937_64& mt_engine = stg::generator_engines::get_engine<std::mt19937_64>(seed);
  std::normal_distribution<> normal_distribution{.0, 1.};
  std::shared_ptr<RNGenerator<std::mt19937_64, std::normal_distribution<>>> rn_generator
    = std::make_shared<RNGenerator<std::mt19937_64, std::normal_distribution<>>>(
        mt_engine,
        std::move(normal_distribution)
      );
};

SCENARIO_METHOD(GeneratorFixture, "Check that generator generate defined values "
                                  "for defined seed, values are not the same") {
  GIVEN("Generator with given seed") {
    const auto first_value = rn_generator->operator()();
    const auto second_value = rn_generator->operator()();

    CHECK_THAT(first_value, Catch::Matchers::WithinRel(0.704988, eps));
    CHECK_THAT(second_value, Catch::Matchers::WithinRel(1.29382, eps));
    CHECK(std::fabs(second_value - first_value) > eps);
  }
}