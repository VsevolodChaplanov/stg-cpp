#ifndef STG_SPACE_COVARIATION_HPP
#define STG_SPACE_COVARIATION_HPP

#include "covariance.hpp"

namespace stg::statistics {

  class SpaceCovariance final {
  public:

    template<NumericViewable FirstRange, NumericViewable SecondRange, NumericViewable ThirdRange>
    static auto covariance_tensor(FirstRange&& f_x_range, FirstRange&& s_x_range,
                                  SecondRange&& f_y_range, SecondRange&& s_y_range,
                                  ThirdRange&& f_z_range, ThirdRange&& s_z_range) {
      auto f_x_mean_f = std::async(std::launch::async, [&] {
        return Mean::mean(f_x_range);
      });
      auto f_y_mean_f = std::async(std::launch::async, [&] {
        return Mean::mean(f_y_range);
      });
      auto f_z_mean_f = std::async(std::launch::async, [&] {
        return Mean::mean(f_z_range);
      });
      auto s_x_mean_f = std::async(std::launch::async, [&] {
        return Mean::mean(s_x_range);
      });
      auto s_y_mean_f = std::async(std::launch::async, [&] {
        return Mean::mean(s_y_range);
      });
      auto s_z_mean_f = std::async(std::launch::async, [&] {
        return Mean::mean(s_z_range);
      });

      auto f_x_mean = f_x_mean_f.get();
      auto f_y_mean = f_y_mean_f.get();
      auto f_z_mean = f_z_mean_f.get();
      auto s_x_mean = s_x_mean_f.get();
      auto s_y_mean = s_y_mean_f.get();
      auto s_z_mean = s_z_mean_f.get();

      return covariance_tensor(f_x_range, s_x_range,
                               f_y_range, s_y_range,
                               f_z_range, s_z_range,
                               f_x_mean, s_x_mean,
                               f_y_mean, s_y_mean,
                               f_z_mean, s_z_mean);
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange,
             NumericViewable ThirdRange, std::floating_point MeanValueType>
    static auto covariance_tensor(FirstRange&& f_x_range, FirstRange&& s_x_range,
                                  SecondRange&& f_y_range, SecondRange&& s_y_range,
                                  ThirdRange&& f_z_range, ThirdRange&& s_z_range,
                                  MeanValueType first_x_mean, MeanValueType second_x_mean,
                                  MeanValueType first_y_mean, MeanValueType second_y_mean,
                                  MeanValueType first_z_mean, MeanValueType second_z_mean) {
      auto c_11_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_x_range, s_x_range,
                                      first_x_mean, second_x_mean); });
      auto c_12_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_x_range, s_y_range,
                                      first_x_mean, second_y_mean); });
      auto c_13_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_x_range, s_z_range,
                                      first_x_mean, second_z_mean); });
      auto c_21_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_y_range, s_x_range,
                                      first_y_mean, second_x_mean); });
      auto c_22_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_y_range, s_y_range,
                                      first_y_mean, second_y_mean); });
      auto c_23_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_y_range, s_z_range,
                                      first_y_mean, second_z_mean); });
      auto c_31_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_z_range, s_x_range,
                                      first_z_mean, second_x_mean); });
      auto c_32_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_z_range, s_y_range,
                                      first_z_mean, second_y_mean); });
      auto c_33_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_z_range, s_z_range,
                                      first_z_mean, second_z_mean); });

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

#endif //STG_SPACE_COVARIATION_HPP
