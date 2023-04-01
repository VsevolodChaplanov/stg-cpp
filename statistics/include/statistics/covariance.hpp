#ifndef STG_STATISTICS_covariance_HPP
#define STG_STATISTICS_covariance_HPP

#include <execution>
#include <ranges>
#include <future>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <range/v3/numeric/inner_product.hpp>
#include <stg_tensor/tensor.hpp>
#include "concepts.hpp"
#include "mean.hpp"
#include "standard_deviation.hpp"


namespace stg::statistics {
using namespace stg::tensor;

  class Covariance {
  public:

    template<NumericViewable FirstRange, NumericViewable SecondRange>
    static auto covariance(FirstRange&& f_range, SecondRange&& s_range) {
      auto f_mean = std::async(std::launch::async, [&] {
        return Mean::mean(s_range);
      });
      auto s_mean = std::async(std::launch::async, [&] {
        return Mean::mean(s_range);
      });

      return covariance(f_range, s_range, f_mean.get(), s_mean.get());
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange,
             std::floating_point MeanValueType>
    static auto covariance(FirstRange&& f_range, SecondRange&& s_range,
                           MeanValueType f_mean, MeanValueType s_mean) {
      const std::size_t f_size = std::ranges::distance(f_range);
      const std::size_t s_size = std::ranges::distance(s_range);
      if (f_size != s_size) {
        throw std::logic_error("Ranges sizes are not equal");
      }

      const auto xy_product = ranges::inner_product(f_range, s_range, 0.);
      const auto mean_xy = xy_product / f_size;
      return mean_xy - f_mean * s_mean;
    }

    template<NumericViewable Range, std::floating_point MeanValueType>
    static auto covariance_tensor(Range&& f_range, Range&& s_range, Range&& t_range,
                                  MeanValueType first_mean, MeanValueType second_mean,
                                  MeanValueType third_mean) {
      auto c_11_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_range, f_range,
                                      first_mean, first_mean); });
      auto c_12_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_range, s_range,
                                      first_mean, second_mean); });
      auto c_13_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(f_range, t_range,
                                      first_mean, third_mean); });
      auto c_22_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(s_range, s_range,
                                      second_mean, second_mean); });
      auto c_23_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(s_range, t_range,
                                      second_mean, third_mean); });
      auto c_33_f = std::async(std::launch::async, [&] {
        return Covariance::covariance(t_range, t_range,
                                      third_mean, third_mean); });

      auto c_11 = c_11_f.get();
      auto c_12 = c_12_f.get();
      auto c_13 = c_13_f.get();
      auto c_22 = c_22_f.get();
      auto c_23 = c_23_f.get();
      auto c_33 = c_33_f.get();

      return Tensor{ c_11, c_12, c_13,
                     c_12, c_22, c_23,
                     c_13, c_23, c_33 };
    }

    template<NumericViewable Range>
    static auto covariance_tensor(Range&& f_range, Range&& s_range, Range&& t_range) {
      const auto f_mean = Mean::mean(std::forward<Range>(f_range));
      const auto s_mean = Mean::mean(std::forward<Range>(s_range));
      const auto t_mean = Mean::mean(std::forward<Range>(t_range));

      return covariance_tensor(std::forward<Range>(f_range),
                               std::forward<Range>(s_range),
                               std::forward<Range>(t_range),
                               f_mean, s_mean, t_mean);
    }

    template<std::forward_iterator Iter>
    static auto covariance(Iter first_begin, Iter first_end,
                           Iter second_begin, Iter second_end) {
      const auto first_mean = Mean::mean(first_begin, first_end);
      const auto second_mean = Mean::mean(second_begin, second_end);

      const auto result = covariance(first_begin, first_end,
                                     second_begin, second_end,
                                     first_mean, second_mean);
      return result;
    }

    template<std::forward_iterator Iter, std::floating_point MeanValueType>
    static auto covariance(Iter first_begin, Iter first_end,
                           Iter second_begin, Iter second_end,
                           MeanValueType first_mean,
                           MeanValueType second_mean) {
      const size_t first_size = std::distance(first_begin, first_end);
      const size_t second_size = std::distance(second_begin, second_end);

      if (first_size != second_size) {
        throw std::logic_error("Collections to which iterators point are of different lengths");
      }

      IterValueType<Iter> sum = 0;
      for (; first_begin != first_end; ++first_begin, ++second_begin) {
        sum += (*first_begin - first_mean) * (*second_begin - second_mean);
      }

      const auto result = sum / first_size;
      return result;
    }

    template<std::forward_iterator Iter>
    static Tensor<IterValueType<Iter>>
    covariance_tensor(Iter first_begin, Iter first_end,
                      Iter second_begin, Iter second_end,
                      Iter third_begin, Iter third_end) {
      const auto first_mean = Mean::mean(first_begin, first_end);
      const auto second_mean = Mean::mean(second_begin, second_end);
      const auto third_mean = Mean::mean(third_begin, third_end);
      auto result = covariance_tensor(first_begin, first_end,
                                      second_begin, second_end,
                                      third_begin, third_end,
                                      first_mean, second_mean, third_mean);
      return result;
    }

    template<std::forward_iterator Iter>
    static Tensor<IterValueType<Iter>>
    covariance_tensor(Iter first_begin, Iter first_end,
                      Iter second_begin, Iter second_end,
                      Iter third_begin, Iter third_end,
                      IterValueType<Iter> first_mean,
                      IterValueType<Iter> second_mean,
                      IterValueType<Iter> third_mean) {
      const auto c_11 = covariance(first_begin, first_end, first_begin, first_end, first_mean, first_mean);
      const auto c_12 = covariance(first_begin, first_end, second_begin, second_end, first_mean, second_mean);
      const auto c_13 = covariance(first_begin, first_end, third_begin, third_end, first_mean, third_mean);
      const auto c_22 = covariance(second_begin, second_end, second_begin, second_end, second_mean, second_mean);
      const auto c_23 = covariance(second_begin, second_end, third_begin, third_end, second_mean, third_mean);
      const auto c_33 = covariance(third_begin, third_end, third_begin, third_end, third_mean, third_mean);

      return { c_11, c_12, c_13, c_12, c_22, c_23, c_13, c_23, c_33 };
    }
  };
}

#endif