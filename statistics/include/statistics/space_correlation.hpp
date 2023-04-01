#ifndef STG_SPACE_CORRELATION_HPP
#define STG_SPACE_CORRELATION_HPP

#include "standard_deviation.hpp"
#include "correlation.hpp"
#include "space_covariation.hpp"

namespace stg::statistics {

  class SpaceCorrelation final {
  public:
    template<NumericViewable FirstRange, NumericViewable SecondRange, NumericViewable ThirdRange>
    static auto correlation_tensor(FirstRange&& f_x_range, FirstRange&& s_x_range,
                                   SecondRange&& f_y_range, SecondRange&& s_y_range,
                                   ThirdRange&& f_z_range, ThirdRange&& s_z_range) {
      auto f_x_mean = std::async(std::launch::async, [&] {
        return Mean::mean(f_x_range);
      });
      auto f_y_mean = std::async(std::launch::async, [&] {
        return Mean::mean(f_y_range);
      });
      auto f_z_mean = std::async(std::launch::async, [&] {
        return Mean::mean(f_z_range);
      });
      auto s_x_mean = std::async(std::launch::async, [&] {
        return Mean::mean(s_x_range);
      });
      auto s_y_mean = std::async(std::launch::async, [&] {
        return Mean::mean(s_y_range);
      });
      auto s_z_mean = std::async(std::launch::async, [&] {
        return Mean::mean(s_z_range);
      });

      return correlation_tensor(std::forward<FirstRange>(f_x_range), std::forward<FirstRange>(s_x_range),
                                std::forward<SecondRange>(f_y_range), std::forward<SecondRange>(s_y_range),
                                std::forward<ThirdRange>(f_z_range), std::forward<ThirdRange>(s_z_range),
                                f_x_mean.get(), s_x_mean.get(),
                                f_y_mean.get(), s_y_mean.get(),
                                f_z_mean.get(), s_z_mean.get());
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange,
             NumericViewable ThirdRange, std::floating_point MeanValueType>
    static auto correlation_tensor(FirstRange&& f_x_range, FirstRange&& s_x_range,
                                   SecondRange&& f_y_range, SecondRange&& s_y_range,
                                   ThirdRange&& f_z_range, ThirdRange&& s_z_range,
                                   MeanValueType first_x_mean, MeanValueType second_x_mean,
                                   MeanValueType first_y_mean, MeanValueType second_y_mean,
                                   MeanValueType first_z_mean, MeanValueType second_z_mean) {
      auto f_x_std = std::async(std::launch::async, [&] {
        return StandardDeviation::std(f_x_range, first_x_mean);
      });
      auto f_y_std = std::async(std::launch::async, [&] {
        return StandardDeviation::std(f_y_range, first_y_mean);
      });
      auto f_z_std = std::async(std::launch::async, [&] {
        return StandardDeviation::std(f_z_range, first_z_mean);
      });
      auto s_x_std = std::async(std::launch::async, [&] {
        return StandardDeviation::std(s_x_range, second_x_mean);
      });
      auto s_y_std = std::async(std::launch::async, [&] {
        return StandardDeviation::std(s_y_range, second_y_mean);
      });
      auto s_z_std = std::async(std::launch::async, [&] {
        return StandardDeviation::std(s_z_range, second_z_mean);
      });

      return correlation_tensor(std::forward<FirstRange>(f_x_range), std::forward<FirstRange>(s_x_range),
                                std::forward<SecondRange>(f_y_range), std::forward<SecondRange>(s_y_range),
                                std::forward<ThirdRange>(f_z_range), std::forward<ThirdRange>(s_z_range),
                                first_x_mean, second_x_mean,
                                first_y_mean, second_y_mean,
                                first_z_mean, second_z_mean,
                                f_x_std.get(), s_x_std.get(),
                                f_y_std.get(), s_y_std.get(),
                                f_z_std.get(), s_z_std.get());
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange,
             NumericViewable ThirdRange, std::floating_point MeanValueType>
    static auto correlation_tensor(FirstRange&& f_x_range, FirstRange&& s_x_range,
                                   SecondRange&& f_y_range, SecondRange&& s_y_range,
                                   ThirdRange&& f_z_range, ThirdRange&& s_z_range,
                                   MeanValueType first_x_mean, MeanValueType second_x_mean,
                                   MeanValueType first_y_mean, MeanValueType second_y_mean,
                                   MeanValueType first_z_mean, MeanValueType second_z_mean,
                                   MeanValueType first_x_std, MeanValueType second_x_std,
                                   MeanValueType first_y_std, MeanValueType second_y_std,
                                   MeanValueType first_z_std, MeanValueType second_z_std) {
      auto c_11_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_x_range, s_x_range,
                                        first_x_mean, second_x_mean,
                                        first_x_std, second_x_std); });
      auto c_12_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_x_range, s_y_range,
                                        first_x_mean, second_y_mean,
                                        first_x_std, second_y_std); });
      auto c_13_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_x_range, s_z_range,
                                        first_x_mean, second_z_mean,
                                        first_x_std, second_z_std); });
      auto c_21_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_y_range, s_x_range,
                                        first_y_mean, second_x_mean,
                                        first_y_std, second_x_std); });
      auto c_22_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_y_range, s_y_range,
                                        first_y_mean, second_y_mean,
                                        first_y_std, second_y_std); });
      auto c_23_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_y_range, s_z_range,
                                        first_y_mean, second_z_mean,
                                        first_y_std, second_z_std); });
      auto c_31_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_z_range, s_x_range,
                                        first_z_mean, second_x_mean,
                                        first_z_std, second_x_std); });
      auto c_32_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_z_range, s_y_range,
                                        first_z_mean, second_y_mean,
                                        first_z_std, second_y_std); });
      auto c_33_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_z_range, s_z_range,
                                        first_z_mean, second_z_mean,
                                        first_z_std, second_z_std); });

      auto c_11 = c_11_f.get();
      auto c_12 = c_12_f.get();
      auto c_13 = c_13_f.get();
      auto c_21 = c_21_f.get();
      auto c_22 = c_22_f.get();
      auto c_23 = c_23_f.get();
      auto c_31 = c_31_f.get();
      auto c_32 = c_32_f.get();
      auto c_33 = c_33_f.get();

      return Tensor{ c_11, c_12, c_13,
                     c_21, c_22, c_23,
                     c_31, c_32, c_33 };
    }
  };
}

#endif //STG_SPACE_CORRELATION_HPP
