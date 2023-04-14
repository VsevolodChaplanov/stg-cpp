#include "common.hpp"

using namespace stg;

struct Fixture {
  constexpr static inline double eps = 1.e-6;
  std::mt19937_64& mt_engine = stg::generator_engines::get_engine<std::mt19937_64>(seed);
  std::normal_distribution<> normal_distribution{.0, 1.};
  RNGenerator<std::mt19937_64, std::normal_distribution<>> rn_generator{mt_engine, normal_distribution};

  constexpr static inline std::array<double, 10> first_sample
    = {0.34820076, 0.26857129, 0.10389639, 0.21398258, 0.94669873,
       0.0425099 , 0.31029466, 0.67020328, 0.04364531, 0.28159116};
  constexpr static inline std::array<double, 10> second_sample
    = {0.18361852, 0.65461521, 0.85717529, 0.42393114, 0.75886066,
       0.15420898, 0.70441173, 0.45965862, 0.91917156, 0.91907483};
  constexpr static inline std::array<double, 10> third_sample
    = {0.66376923, 0.07483858, 0.90539227, 0.09945824, 0.42822669,
       0.30501249, 0.35339393, 0.66808912, 0.03961761, 0.356688326};

  constexpr static inline double first_mean = 0.3229594060975379;
  constexpr static inline double second_mean = 0.6034726541690436;
  constexpr static inline double third_mean = 0.389448647674321;

  constexpr static inline double first_std = 0.27048638106834644;
  constexpr static inline double second_std = 0.2704844207066285;
  constexpr static inline double third_std = 0.2712106545412239;

  constexpr static inline Tensor<double> covariation_tensor{
    std::array{0.07316288, 0.00210752, 0.01867616, 0.00210752,
      0.07316182, -0.00688355, 0.01867616, -0.00688355, 0.07355522}
  };

  constexpr static inline Tensor<double> correlation_tensor{
    std::array{1., 0.02880604, 0.25458655,
               0.02880604, 1., -0.09383473,
               0.25458655, -0.09383473, 1.}
  };
};

SCENARIO_METHOD(Fixture, "Check statistical values equality") {
  const auto is_floating_equal = [](const auto& first_seq_elem, const auto& second_seq_elem){
    return std::fabs(first_seq_elem - second_seq_elem) < eps;
  };
  GIVEN("Values calculated by algorithms") {
    const auto test_first_mean = statistics::Mean::mean(first_sample.cbegin(), first_sample.cend());
    const auto test_second_mean = statistics::Mean::mean(second_sample.cbegin(), second_sample.cend());
    const auto test_third_mean = statistics::Mean::mean(third_sample.cbegin(), third_sample.cend());

    const auto test_first_std = statistics::StandardDeviation::std(first_sample.cbegin(), first_sample.cend());
    const auto test_second_std = statistics::StandardDeviation::std(second_sample.cbegin(), second_sample.cend());
    const auto test_third_std = statistics::StandardDeviation::std(third_sample.cbegin(), third_sample.cend());

    const auto test_cov_tensor = statistics::Covariance::covariance_tensor(first_sample.cbegin(), first_sample.cend(),
                                                                           second_sample.cbegin(), second_sample.cend(),
                                                                           third_sample.cbegin(), third_sample.cend());
    const auto test_cor_tensor = statistics::Correlation::correlation_tensor(first_sample.cbegin(), first_sample.cend(),
                                                                             second_sample.cbegin(), second_sample.cend(),
                                                                             third_sample.cbegin(), third_sample.cend());

    CHECK_THAT(test_first_mean, Catch::Matchers::WithinRel(first_mean, eps));
    CHECK_THAT(test_second_mean, Catch::Matchers::WithinRel(second_mean, eps));
    CHECK_THAT(test_third_mean, Catch::Matchers::WithinRel(third_mean, eps));

    CHECK_THAT(test_first_std, Catch::Matchers::WithinRel(first_std, eps));
    CHECK_THAT(test_second_std, Catch::Matchers::WithinRel(second_std, eps));
    CHECK_THAT(test_third_std, Catch::Matchers::WithinRel(third_std, eps));

    const bool is_covariations_equal = std::equal(test_cov_tensor.cbegin(), test_cov_tensor.cend(),
                                                  covariation_tensor.cbegin(),
                                                  is_floating_equal);
    const bool is_correlation_equal = std::equal(test_cor_tensor.cbegin(), test_cor_tensor.cend(), correlation_tensor.cbegin(),
                                                 is_floating_equal);

    CHECK(is_covariations_equal);
    CHECK(is_correlation_equal);
  }
}