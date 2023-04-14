#include "common.hpp"

struct RangesSpaceStatisticsFixture {
  constexpr static inline double eps = 1.e-5;

  std::vector<double> first_x_sample
    = {0.34820076, 0.26857129, 0.10389639, 0.21398258, 0.94669873,
       0.0425099, 0.31029466, 0.67020328, 0.04364531, 0.28159116};
  std::vector<double> first_y_sample
    = {0.18361852, 0.65461521, 0.85717529, 0.42393114, 0.75886066,
       0.15420898, 0.70441173, 0.45965862, 0.91917156, 0.91907483};
  std::vector<double> first_z_sample
    = {0.66376923, 0.07483858, 0.90539227, 0.09945824, 0.42822669,
       0.30501249, 0.35339393, 0.66808912, 0.03961761, 0.356688326};

  std::vector<double> second_x_sample
    = {-0.31329386,  0.09340065, -0.81094991, -2.00560014, -0.55417502,
       0.33007673,  0.03751994,  0.74763207,  0.68756813,  2.11817622};

  std::vector<double> second_y_sample
    = {0.62454724, -0.67122117, -2.00168978, -0.59064146,  0.76010282,
       -0.35949893,  0.31602925,  0.42544118, -0.34743041, -0.21730326};

  std::vector<double> second_z_sample
    = {-1.48821841, -0.583823  ,  0.61463187, -0.79433738, -0.89903332,
       0.54243189,  0.52758341,  0.86055713,  1.16767381,  1.39355251};

  Tensor<double> covariation_tensor{{-0.00755426,  0.13798217, -0.08584327,
                                      0.0838092 , -0.06262648,  0.12480679,
                                     -0.01044446, -0.01833379,  0.0008852}};

  Tensor<double> correlation_tensor{{-0.02698447,  0.66261494, -0.33809468,
                                      0.29937575, -0.3007457,   0.49155653,
                                     -0.03720888, -0.08780703,  0.00347708}};
  static inline const auto is_floating_equal = [](const auto& first_seq_elem, const auto& second_seq_elem){
    CHECK_THAT(first_seq_elem, WithinRel(second_seq_elem, eps));
    return true;
  };
};

SCENARIO_METHOD(RangesSpaceStatisticsFixture, "Calculate space covariation") {
  GIVEN("Calculated space covariation") {
    const auto test_cov = SpaceCovariance::covariance_tensor(
      first_x_sample, second_x_sample,
      first_y_sample, second_y_sample,
      first_z_sample, second_z_sample);

    CHECK(std::equal(covariation_tensor.cbegin(), covariation_tensor.cend(),
                     test_cov.cbegin(), is_floating_equal));
  }

  GIVEN("Calculated space correlation") {
    const auto test_corr = SpaceCorrelation::correlation_tensor(
      first_x_sample, second_x_sample,
      first_y_sample, second_y_sample,
      first_z_sample, second_z_sample);

    CHECK(std::equal(correlation_tensor.cbegin(), correlation_tensor.cend(),
                     test_corr.cbegin(), is_floating_equal));
  }
}
