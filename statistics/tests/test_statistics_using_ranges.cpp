#include "common.hpp"

struct RangesStatisticsFixture {
  constexpr static inline double eps = 1.e-6;

  std::vector<double> first_sample
  = {0.34820076, 0.26857129, 0.10389639, 0.21398258, 0.94669873,
     0.0425099 , 0.31029466, 0.67020328, 0.04364531, 0.28159116};
  std::vector<double> second_sample
  = {0.18361852, 0.65461521, 0.85717529, 0.42393114, 0.75886066,
     0.15420898, 0.70441173, 0.45965862, 0.91917156, 0.91907483};
  std::vector<double> third_sample
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

  inline static const auto check_ranges_equal = [](auto f_range, auto s_range) -> void {
    const std::size_t f_size = ranges::size(f_range);
    const std::size_t s_size = ranges::size(s_range);

    REQUIRE(f_size == s_size);

    const auto zip_view = ranges::zip_view(f_range, s_range);

    ranges::for_each(zip_view, [](const auto& pair) {
      const auto [x, y] = pair;
      CHECK_THAT(x, WithinRel(y, eps));
    });
  };

  inline static const auto is_floating_equal = [](const auto& first_seq_elem, const auto& second_seq_elem){
    return std::fabs(first_seq_elem - second_seq_elem) < eps;
  };
};


SCENARIO_METHOD(RangesStatisticsFixture, "Calculate statistics for ranges") {
  GIVEN("Copied samples to check to the ranges doesn't change values") {
    auto f_sample_cpy = first_sample;
    auto s_sample_cpy = second_sample;
    auto t_sample_cpy = third_sample;

    THEN("Calculate means and check") {
      const auto f_test_mean = Mean::mean(f_sample_cpy);
      const auto s_test_mean = Mean::mean(s_sample_cpy);
      const auto t_test_mean = Mean::mean(t_sample_cpy);

      CHECK_THAT(f_test_mean, WithinRel(first_mean, eps));
      CHECK_THAT(s_test_mean, WithinRel(second_mean, eps));
      CHECK_THAT(t_test_mean, WithinRel(third_mean, eps));

      AND_THEN("Check that the initial array doesn't changed") {
        CHECK_THAT(first_sample, Approx(f_sample_cpy).margin(eps));
        CHECK_THAT(second_sample, Approx(s_sample_cpy).margin(eps));
        CHECK_THAT(third_sample, Approx(t_sample_cpy).margin(eps));
      }

      AND_THEN("Take a view of samples and do the same") {
        auto f_sample_view = std::ranges::views::all(first_sample);
        auto s_sample_view = ranges::views::all(second_sample);
        auto t_sample_view = ranges::views::all(third_sample);

        const auto f_view_mean = Mean::mean(f_sample_view);
        const auto s_view_mean = Mean::mean(s_sample_view);
        const auto t_view_mean = Mean::mean(t_sample_view);

        CHECK_THAT(f_view_mean, WithinRel(first_mean, eps));
        CHECK_THAT(s_view_mean, WithinRel(second_mean, eps));
        CHECK_THAT(t_view_mean, WithinRel(third_mean, eps));

        AND_THEN("Check that the initial array doesn't changed") {
          check_ranges_equal(f_sample_view, first_sample);
          check_ranges_equal(s_sample_view, second_sample);
          check_ranges_equal(t_sample_view, third_sample);
        }
      }

      AND_THEN("Check that can calculate mean using initial array") {
        const auto f_s_mean = Mean::mean(first_sample);
        const auto s_s_mean = Mean::mean(second_sample);
        const auto t_s_mean = Mean::mean(third_sample);

        CHECK_THAT(f_s_mean, WithinRel(first_mean, eps));
        CHECK_THAT(s_s_mean, WithinRel(second_mean, eps));
        CHECK_THAT(t_s_mean, WithinRel(third_mean, eps));
      }
    }
  }

  GIVEN("Copied ranges to check that algos for std don't change values") {
    auto f_sample_cpy = first_sample;
    auto s_sample_cpy = second_sample;
    auto t_sample_cpy = third_sample;

    THEN("Calculate standard deviations and check") {
      const auto f_test_std = StandardDeviation::std(f_sample_cpy);
      const auto s_test_std = StandardDeviation::std(s_sample_cpy);
      const auto t_test_std = StandardDeviation::std(t_sample_cpy);

      CHECK_THAT(f_test_std, WithinRel(first_std, eps));
      CHECK_THAT(s_test_std, WithinRel(second_std, eps));
      CHECK_THAT(t_test_std, WithinRel(third_std, eps));

      AND_THEN("Check that the initial array doesn't changed") {
        CHECK_THAT(first_sample, Approx(f_sample_cpy).margin(eps));
        CHECK_THAT(second_sample, Approx(s_sample_cpy).margin(eps));
        CHECK_THAT(third_sample, Approx(t_sample_cpy).margin(eps));
      }

      AND_THEN("Take a view of samples and do the same") {
        auto f_sample_view = std::ranges::views::all(first_sample);
        auto s_sample_view = ranges::views::all(second_sample);
        auto t_sample_view = ranges::views::all(third_sample);

        const auto f_view_std = StandardDeviation::std(f_sample_view);
        const auto s_view_std = StandardDeviation::std(s_sample_view);
        const auto t_view_std = StandardDeviation::std(t_sample_view);

        CHECK_THAT(f_view_std, WithinRel(first_std, eps));
        CHECK_THAT(s_view_std, WithinRel(second_std, eps));
        CHECK_THAT(t_view_std, WithinRel(third_std, eps));

        AND_THEN("Check that the initial array doesn't changed") {
          check_ranges_equal(f_sample_view, first_sample);
          check_ranges_equal(s_sample_view, second_sample);
          check_ranges_equal(t_sample_view, third_sample);
        }
      }

      AND_THEN("Check that can calculate mean using initial array") {
        const auto f_s_std = StandardDeviation::std(first_sample);
        const auto s_s_std = StandardDeviation::std(second_sample);
        const auto t_s_std = StandardDeviation::std(third_sample);

        CHECK_THAT(f_s_std, WithinRel(first_std, eps));
        CHECK_THAT(s_s_std, WithinRel(second_std, eps));
        CHECK_THAT(t_s_std, WithinRel(third_std, eps));
      }
    }
  }
}

SCENARIO_METHOD(RangesStatisticsFixture, "Check ranges algos for covariance matrices") {
  GIVEN("Copied samples") {
    auto f_sample_cpy = first_sample;
    auto s_sample_cpy = second_sample;
    auto t_sample_cpy = third_sample;

    THEN("Check covariation tensor for copied samples and compare") {
      const auto test_covariation_tensor = Covariance::covariance_tensor(f_sample_cpy, s_sample_cpy, t_sample_cpy);

      const bool is_covariations_equal = std::equal(test_covariation_tensor.cbegin(), test_covariation_tensor.cend(),
                                                    covariation_tensor.cbegin(),
                                                    is_floating_equal);
      CHECK(is_covariations_equal);
      check_ranges_equal(f_sample_cpy, first_sample);
      check_ranges_equal(s_sample_cpy, second_sample);
      check_ranges_equal(t_sample_cpy, third_sample);

      AND_THEN("Create views and check") {
        auto f_sample_view = ranges::views::all(f_sample_cpy);
        auto s_sample_view = ranges::views::all(s_sample_cpy);
        auto t_sample_view = ranges::views::all(t_sample_cpy);

        const auto covariation_views = Covariance::covariance_tensor(f_sample_view, s_sample_view, t_sample_view);
        const bool is_covariations_equal = std::equal(covariation_views.cbegin(), covariation_views.cend(),
                                                      covariation_tensor.cbegin(),
                                                      is_floating_equal);
        CHECK(is_covariations_equal);
        check_ranges_equal(f_sample_cpy, first_sample);
        check_ranges_equal(s_sample_cpy, second_sample);
        check_ranges_equal(t_sample_cpy, third_sample);
      }
    }
  }
}


SCENARIO_METHOD(RangesStatisticsFixture, "Check ranges algos for correlation matrices") {
  GIVEN("Copied samples") {
    auto f_sample_cpy = first_sample;
    auto s_sample_cpy = second_sample;
    auto t_sample_cpy = third_sample;

    THEN("Check covariation tensor for copied samples and compare") {
      const auto corr_views = Correlation::correlation_tensor(f_sample_cpy, s_sample_cpy, t_sample_cpy);

      const bool is_corr_equal = std::equal(corr_views.cbegin(), corr_views.cend(),
                                            correlation_tensor.cbegin(), is_floating_equal);
      CHECK(is_corr_equal);
      check_ranges_equal(f_sample_cpy, first_sample);
      check_ranges_equal(s_sample_cpy, second_sample);
      check_ranges_equal(t_sample_cpy, third_sample);

      AND_THEN("Create views and check") {
        auto f_sample_view = ranges::views::all(f_sample_cpy);
        auto s_sample_view = ranges::views::all(s_sample_cpy);
        auto t_sample_view = ranges::views::all(t_sample_cpy);

        const auto corr_views = Correlation::correlation_tensor(f_sample_view, s_sample_view, t_sample_view);
        const bool is_corr_equal = std::equal(corr_views.cbegin(), corr_views.cend(),
                                              correlation_tensor.cbegin(), is_floating_equal);
        CHECK(is_corr_equal);
        check_ranges_equal(f_sample_cpy, first_sample);
        check_ranges_equal(s_sample_cpy, second_sample);
        check_ranges_equal(t_sample_cpy, third_sample);
      }
    }
  }
}